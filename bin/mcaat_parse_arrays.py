#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mcaat_parse_arrays.py — parse MCAAT ``CRISPR_Arrays_*.txt`` into machine-readable tables.

This is the ONLY place in the pipeline that knows the MCAAT array output grammar.
Nothing else -- no workflow, no config, no module -- may hard-code a column name or a
format assumption about ``CRISPR_Arrays_*.txt``.

Authoritative grammar (verified against the writer in
``include/crispr_array/post_processing.h`` -- ``PostProcessor``'s "9: Output" step, the
``open_new_file()`` lambda and the ``consensus_arrays`` loop that follows it -- and
against MCAAT's own reference reader in ``include/core/array_to_nodes.h``)::

    # MCAAT — CRISPR Array Output            <- literal U+2014 EM DASH, UTF-8
    # Generated : 2026-07-21 11:04:33
    # Arrays    : 412                        <- GLOBAL run total, repeated verbatim
    # Spacers   : 3117                       <- GLOBAL run total, repeated verbatim
                                             <- blank line (5-line preamble)
    >Array_7  spacers=8                      <- '>' Array_ N, two spaces, spacers=<count>
    GTTTCAATCCACGCACCCGTGGGGGTGATCC          <- consensus repeat, column 0
            ------------------------------\tAGGCTTAAGCCTGAGGCTAAGCTAG
            ------------------------------\tTTGACCTGAAGCTTAGCCTGAAGCT
                                             <- blank line closes the block
    >Array_8  spacers=3
    ...

Spacer lines are EXACTLY eight leading spaces, then the repeat variant, then a single
TAB, then the spacer.  The writer emits::

    const std::string& rv = (repeat_var == arr.consensus)
        ? std::string(arr.consensus.size(), '-')
        : repeat_var;

i.e. a run of '-' of exactly the consensus length is the on-disk encoding of
"this repeat is identical to the consensus".  It MUST be expanded back to the real
consensus sequence in every downstream table -- emitting the literal dashes as a
sequence is meaningless.

Rules implemented (R1-R10):

  R1  Input files are sorted NUMERICALLY on the ``_<N>`` suffix.  MCAAT's own reader
      sorts lexicographically, so ``_10.txt`` precedes ``_2.txt``; array numbering is
      global and continuous, so numeric order reproduces the original sequence.
  R2  The 5-line ``#`` preamble is stripped from EVERY file.  ``# Arrays`` / ``# Spacers``
      are GLOBAL run totals repeated verbatim in each split file -- they are NEVER
      summed, only used as a cross-check.
  R3  Spacer lines are split on the FIRST TAB, never on whitespace runs.  (readme.md
      shows spaces; ``post_processing.h`` writes ``\\t``.)
  R4  A repeat variant that is entirely '-' is expanded to the consensus verbatim.
  R5  A non-dashed variant is NOT assumed to be the same length as the consensus
      (the consensus is ``tail_consensus + front_consensus``; a variant is the
      untailed ``rep`` and is typically shorter).
  R6  The variant column is block-constant in practice; a mixed block is reported as
      a warning, never as an error.
  R7  Files are opened with ``encoding='utf-8'`` (the header carries an em-dash).
  R8  Blank lines and ``#`` lines are skipped; ``>`` opens a block; the next
      non-indented line is the consensus; indented lines are spacer lines; an
      indented line with no TAB is malformed and is skipped with a warning.
  R9  ``spacers=N`` must equal the number of following indented lines; a mismatch is
      a warning.
  R10 All three empty states (missing file / header-only / populated-but-filtered)
      normalise to zero data rows plus a summary row -- never to an error.  Zero
      arrays is a common, legitimate metagenomic result.

Additionally performs the ``parameters.json`` readback (``--parameters`` +
``--expected``).  MCAAT's argument parse loop has no unknown-argument branch, so a
typo'd or dropped flag is silently ignored; comparing the values MCAAT actually
recorded against the values this pipeline intended is the only detector for that.

Standard library only.  Python >= 3.8.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from statistics import median
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------

#: Column order of ``<sample>.arrays.tsv``.  Downstream code imports this name rather
#: than re-declaring the strings.
ARRAYS_TSV_COLUMNS: Tuple[str, ...] = (
    "sample",
    "array_id",
    "source_file",
    "consensus_repeat",
    "consensus_repeat_len",
    "n_spacers",
    "spacer_index",
    "repeat_variant",
    "repeat_variant_is_consensus",
    "spacer",
    "spacer_len",
)

#: Column order of ``<sample>.provenance.tsv``.
PROVENANCE_TSV_COLUMNS: Tuple[str, ...] = (
    "sample",
    "mcaat_threads",
    "mcaat_ram_gb",
    "mcaat_autoclean",
    "mcaat_threshold_multiplicity",
    "mcaat_cycle_min_length",
    "mcaat_cycle_max_length",
    "mcaat_low_abundance",
    "mcaat_spacer_min_length",
    "mcaat_spacer_max_length",
    "mcaat_repeat_min_length",
    "mcaat_repeat_max_length",
)

#: Mapping from the ``--expected`` JSON key to the dotted path inside ``parameters.json``.
#: ``output_file`` and ``cycles_folder`` are deliberately absent: MCAAT advertises paths
#: there that it never creates.
EXPECTED_TO_JSON_PATH: Dict[str, str] = {
    "threads": "threads",
    "autoclean": "autoclean",
    "threshold_multiplicity": "cycle_finder.threshold_multiplicity",
    "cycle_min_length": "cycle_finder.cycle_min_length",
    "cycle_max_length": "cycle_finder.cycle_max_length",
    "low_abundance": "cycle_finder.low_abundance",
    "spacer_min_length": "dna_sequence.spacer_min_length",
    "spacer_max_length": "dna_sequence.spacer_max_length",
    "repeat_min_length": "dna_sequence.repeat_min_length",
    "repeat_max_length": "dna_sequence.repeat_max_length",
}

_ARRAY_HEADER_RE = re.compile(r"^>\s*Array[_ ]?(\d+)\s*(?:spacers\s*=\s*(\d+))?\s*$")
_FILE_INDEX_RE = re.compile(r"_(\d+)\s*$")
_PREAMBLE_ARRAYS_RE = re.compile(r"^#\s*Arrays\s*:\s*(\d+)\s*$")
_PREAMBLE_SPACERS_RE = re.compile(r"^#\s*Spacers\s*:\s*(\d+)\s*$")
_PREAMBLE_GENERATED_RE = re.compile(r"^#\s*Generated\s*:\s*(.*?)\s*$")


def warn(msg: str) -> None:
    """Emit a namespaced warning on stderr.  Never fatal."""
    sys.stderr.write("MCAAT_PARSE_ARRAYS: WARNING: {0}\n".format(msg))


def info(msg: str) -> None:
    sys.stderr.write("MCAAT_PARSE_ARRAYS: {0}\n".format(msg))


def die(msg: str, code: int = 1) -> "NoReturn":  # type: ignore[valid-type]
    sys.stderr.write("MCAAT_PARSE_ARRAYS: ERROR: {0}\n".format(msg))
    sys.exit(code)


# --------------------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------------------


class SpacerEntry(object):
    """One (repeat variant, spacer) pair inside an array block."""

    __slots__ = ("index", "repeat_variant", "variant_is_consensus", "spacer")

    def __init__(self, index: int, repeat_variant: str, variant_is_consensus: bool, spacer: str):
        self.index = index
        self.repeat_variant = repeat_variant
        self.variant_is_consensus = variant_is_consensus
        self.spacer = spacer


class ArrayBlock(object):
    """One ``>Array_<N>`` block."""

    __slots__ = ("array_id", "declared_n_spacers", "consensus", "source_file", "entries")

    def __init__(self, array_id: int, declared_n_spacers: Optional[int], source_file: str):
        self.array_id = array_id
        self.declared_n_spacers = declared_n_spacers
        self.consensus = ""
        self.source_file = source_file
        self.entries: List[SpacerEntry] = []


class FileHeader(object):
    """The 5-line ``#`` preamble of one split file."""

    __slots__ = ("generated", "declared_arrays", "declared_spacers")

    def __init__(self) -> None:
        self.generated: Optional[str] = None
        self.declared_arrays: Optional[int] = None
        self.declared_spacers: Optional[int] = None


# --------------------------------------------------------------------------------------
# File discovery and ordering (R1)
# --------------------------------------------------------------------------------------


def numeric_file_key(path: str) -> Tuple[int, str]:
    """Sort key extracting the ``_<N>`` suffix of ``CRISPR_Arrays_<N>.txt``.

    Files without a numeric suffix sort last, then alphabetically, so the ordering is
    always total and deterministic.
    """
    stem = os.path.splitext(os.path.basename(path))[0]
    match = _FILE_INDEX_RE.search(stem)
    if match is None:
        return (1 << 30, os.path.basename(path))
    return (int(match.group(1)), os.path.basename(path))


def discover_array_files(inputs: Sequence[str]) -> List[str]:
    """Expand directories to their ``*CRISPR_Arrays_*.txt`` members, then sort (R1)."""
    found: List[str] = []
    for item in inputs:
        if os.path.isdir(item):
            for name in os.listdir(item):
                full = os.path.join(item, name)
                if not os.path.isfile(full):
                    continue
                if "CRISPR_Arrays_" in name and name.endswith(".txt"):
                    found.append(full)
        elif os.path.isfile(item):
            found.append(item)
        else:
            warn("input path does not exist and is ignored: {0}".format(item))
    # De-duplicate while preserving the numeric ordering.
    unique = sorted(set(found), key=numeric_file_key)
    return unique


# --------------------------------------------------------------------------------------
# Grammar (R2-R9)
# --------------------------------------------------------------------------------------


def is_dash_run(text: str) -> bool:
    """True when ``text`` is a non-empty run consisting only of '-' (R4)."""
    return len(text) > 0 and text.count("-") == len(text)


def parse_array_file(path: str) -> Tuple[List[ArrayBlock], FileHeader]:
    """Parse one ``CRISPR_Arrays_<N>.txt``.

    Returns the array blocks in file order plus the parsed ``#`` preamble.  Never
    raises on malformed content: unparseable lines are skipped with a warning so that
    one corrupt block cannot destroy an otherwise usable sample (R8, R10).
    """
    blocks: List[ArrayBlock] = []
    header = FileHeader()
    basename = os.path.basename(path)

    current: Optional[ArrayBlock] = None
    consensus_set = False

    # R7: the preamble carries a literal em-dash; errors='replace' guarantees that a
    # byte-level corruption downgrades to a mangled character rather than an exception.
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")

            if line.strip() == "":
                # A blank line closes the block.  MCAAT's own reader merely skips it and
                # relies on the next '>' -- we accept both, which is strictly safer.
                if current is not None:
                    blocks.append(current)
                    current = None
                    consensus_set = False
                continue

            if line.startswith("#"):
                # R2: preamble only.  NEVER summed -- these are global run totals.
                match = _PREAMBLE_ARRAYS_RE.match(line)
                if match is not None:
                    header.declared_arrays = int(match.group(1))
                    continue
                match = _PREAMBLE_SPACERS_RE.match(line)
                if match is not None:
                    header.declared_spacers = int(match.group(1))
                    continue
                match = _PREAMBLE_GENERATED_RE.match(line)
                if match is not None:
                    header.generated = match.group(1)
                continue

            if line.startswith(">"):
                if current is not None:
                    blocks.append(current)
                match = _ARRAY_HEADER_RE.match(line)
                if match is None:
                    warn(
                        "{0}:{1}: unparseable array header {2!r}; block skipped".format(
                            basename, lineno, line
                        )
                    )
                    current = None
                    consensus_set = False
                    continue
                declared = int(match.group(2)) if match.group(2) is not None else None
                current = ArrayBlock(int(match.group(1)), declared, basename)
                consensus_set = False
                continue

            if current is None:
                warn(
                    "{0}:{1}: content outside any array block; skipped: {2!r}".format(
                        basename, lineno, line[:60]
                    )
                )
                continue

            indented = line[0] == " " or line[0] == "\t"

            if not consensus_set and not indented:
                # First non-indented, non-comment line of the block is the consensus.
                current.consensus = line.rstrip(" \t")
                consensus_set = True
                continue

            if not indented:
                # Non-indented line after the consensus.  MCAAT's reference reader
                # silently ignores this; we warn so a format change is visible.
                warn(
                    "{0}:{1}: unexpected non-indented line inside array {2}; skipped: "
                    "{3!r}".format(basename, lineno, current.array_id, line[:60])
                )
                continue

            # R3: split on the FIRST tab, never on whitespace runs.
            tab = line.find("\t")
            if tab < 0:
                warn(
                    "{0}:{1}: malformed spacer line (no TAB) in array {2}; skipped".format(
                        basename, lineno, current.array_id
                    )
                )
                continue

            variant_raw = line[:tab].strip(" \t")
            spacer = line[tab + 1 :].strip(" \t\r\n")
            if spacer == "":
                warn(
                    "{0}:{1}: empty spacer in array {2}; skipped".format(
                        basename, lineno, current.array_id
                    )
                )
                continue

            # R4: a pure dash run means "identical to the consensus"; expand it back.
            if is_dash_run(variant_raw):
                if current.consensus and len(variant_raw) != len(current.consensus):
                    warn(
                        "{0}:{1}: dash run of length {2} does not match consensus length "
                        "{3} in array {4}; expanded to the consensus anyway".format(
                            basename,
                            lineno,
                            len(variant_raw),
                            len(current.consensus),
                            current.array_id,
                        )
                    )
                variant = current.consensus
                variant_is_consensus = True
            else:
                # R5: do NOT assume a non-dashed variant matches the consensus length.
                variant = variant_raw
                variant_is_consensus = variant_raw == current.consensus

            current.entries.append(
                SpacerEntry(
                    index=len(current.entries) + 1,
                    repeat_variant=variant,
                    variant_is_consensus=variant_is_consensus,
                    spacer=spacer,
                )
            )

    if current is not None:
        blocks.append(current)

    # R9 and R6 checks.
    for block in blocks:
        if block.declared_n_spacers is not None and block.declared_n_spacers != len(block.entries):
            warn(
                "{0}: array {1} declares spacers={2} but {3} spacer lines were parsed "
                "(format may have changed)".format(
                    basename, block.array_id, block.declared_n_spacers, len(block.entries)
                )
            )
        if block.consensus == "":
            warn(
                "{0}: array {1} has no consensus repeat line; rows are still emitted with "
                "an empty consensus".format(basename, block.array_id)
            )
        flags = set(entry.variant_is_consensus for entry in block.entries)
        if len(flags) > 1:
            warn(
                "{0}: array {1} mixes dashed and explicit repeat variants; the variant "
                "column is normally block-constant".format(basename, block.array_id)
            )

    return blocks, header


# --------------------------------------------------------------------------------------
# parameters.json readback
# --------------------------------------------------------------------------------------


def load_json_defensively(path: str) -> Optional[dict]:
    """Load a JSON document, downgrading any failure to a warning.

    ``parameters.json`` is hand-serialised by MCAAT with only ``"`` and ``\\`` escaped,
    so a control character or an invalid UTF-8 byte in a path would emit invalid JSON.
    A decode failure must not destroy an otherwise successful run.
    """
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return json.load(handle)
    except (OSError, ValueError) as exc:
        warn("could not parse {0}: {1}".format(path, exc))
        return None


def dotted_get(document: dict, path: str) -> Tuple[bool, object]:
    """Fetch ``a.b.c`` from a nested dict.  Returns ``(found, value)``."""
    node: object = document
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            return (False, None)
        node = node[part]
    return (True, node)


def values_equal(expected: object, observed: object) -> bool:
    """Compare an expected value against what MCAAT recorded, tolerating type drift."""
    if isinstance(expected, bool) or isinstance(observed, bool):
        def as_bool(value: object) -> Optional[bool]:
            if isinstance(value, bool):
                return value
            if isinstance(value, str):
                lowered = value.strip().lower()
                if lowered in ("true", "yes", "1"):
                    return True
                if lowered in ("false", "no", "0"):
                    return False
                return None
            if isinstance(value, (int, float)):
                return bool(value)
            return None

        left = as_bool(expected)
        right = as_bool(observed)
        return left is not None and left == right

    if isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        return abs(float(expected) - float(observed)) < 1e-9

    try:
        return abs(float(str(expected)) - float(str(observed))) < 1e-9
    except (TypeError, ValueError):
        return str(expected) == str(observed)


def readback_parameters(parameters: dict, expected: dict) -> List[str]:
    """Return a list of human-readable mismatch descriptions (empty when all agree)."""
    mismatches: List[str] = []
    for key, value in sorted(expected.items()):
        json_path = EXPECTED_TO_JSON_PATH.get(key)
        if json_path is None:
            warn(
                "expected-parameter key {0!r} is not part of the parameters.json contract "
                "and is ignored".format(key)
            )
            continue
        found, observed = dotted_get(parameters, json_path)
        if not found:
            mismatches.append(
                "{0}: expected {1!r} but parameters.json has no key {2!r}".format(
                    key, value, json_path
                )
            )
            continue
        if not values_equal(value, observed):
            mismatches.append(
                "{0} ({1}): expected {2!r} but MCAAT recorded {3!r}".format(
                    key, json_path, value, observed
                )
            )
    return mismatches


def provenance_row(sample: str, parameters: Optional[dict]) -> Dict[str, str]:
    """Flatten ``parameters.json`` into the provenance echo consumed by the aggregator."""
    row = dict((column, "NA") for column in PROVENANCE_TSV_COLUMNS)
    row["sample"] = sample
    if not parameters:
        return row
    mapping = (
        ("mcaat_threads", "threads"),
        ("mcaat_ram_gb", "ram_gb"),
        ("mcaat_autoclean", "autoclean"),
        ("mcaat_threshold_multiplicity", "cycle_finder.threshold_multiplicity"),
        ("mcaat_cycle_min_length", "cycle_finder.cycle_min_length"),
        ("mcaat_cycle_max_length", "cycle_finder.cycle_max_length"),
        ("mcaat_low_abundance", "cycle_finder.low_abundance"),
        ("mcaat_spacer_min_length", "dna_sequence.spacer_min_length"),
        ("mcaat_spacer_max_length", "dna_sequence.spacer_max_length"),
        ("mcaat_repeat_min_length", "dna_sequence.repeat_min_length"),
        ("mcaat_repeat_max_length", "dna_sequence.repeat_max_length"),
    )
    for column, json_path in mapping:
        found, value = dotted_get(parameters, json_path)
        if not found:
            continue
        if isinstance(value, bool):
            row[column] = "true" if value else "false"
        else:
            row[column] = str(value)
    return row


# --------------------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------------------


def sanitise_field(value: str) -> str:
    """Strip characters that would corrupt a TSV cell."""
    return value.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def write_arrays_tsv(path: str, sample: str, blocks: Sequence[ArrayBlock]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(ARRAYS_TSV_COLUMNS) + "\n")
        for block in blocks:
            consensus = block.consensus
            n_spacers = (
                block.declared_n_spacers
                if block.declared_n_spacers is not None
                else len(block.entries)
            )
            for entry in block.entries:
                handle.write(
                    "\t".join(
                        (
                            sanitise_field(sample),
                            str(block.array_id),
                            sanitise_field(block.source_file),
                            sanitise_field(consensus),
                            str(len(consensus)),
                            str(n_spacers),
                            str(entry.index),
                            sanitise_field(entry.repeat_variant),
                            "true" if entry.variant_is_consensus else "false",
                            sanitise_field(entry.spacer),
                            str(len(entry.spacer)),
                        )
                    )
                    + "\n"
                )


def write_spacers_fasta(path: str, sample: str, blocks: Sequence[ArrayBlock]) -> None:
    """Write one unwrapped record per spacer.

    The identifier is ``<sample>|array<array_id>|sp<spacer_index>``.  It is stable across
    reruns for a fixed input because ``array_id`` comes from MCAAT's own global numbering
    and ``spacer_index`` is the 1-based position inside the block.

    The identifier is followed by DESCRIPTION ATTRIBUTES -- ``repeat_len=`` and
    ``spacer_len=``, separated from the identifier by a single space.  That is deliberate
    and part of the output contract: everything after the first whitespace is the FASTA
    description, so tools that key on the identifier see only the ``|``-delimited part
    while a human reading the file still gets the dimensions.  A test that anchors a
    regular expression on the identifier must therefore not anchor it at end-of-line.
    """
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for block in blocks:
            for entry in block.entries:
                handle.write(
                    ">{sample}|array{array}|sp{index} repeat_len={rlen} spacer_len={slen}\n".format(
                        sample=sample,
                        array=block.array_id,
                        index=entry.index,
                        rlen=len(entry.repeat_variant),
                        slen=len(entry.spacer),
                    )
                )
                handle.write(entry.spacer + "\n")


def write_genstats_mqc(
    path: str,
    sample: str,
    n_arrays: int,
    n_spacers: int,
    n_spacers_unique: int,
    median_spacer_len: float,
    median_repeat_len: float,
) -> None:
    """Write the MultiQC general-statistics custom-content file.

    Uses the zero-config magic-filename mechanism: any file matching
    ``*_mqc.tsv`` is discovered automatically, and files sharing an ``id`` merge across
    samples into one section.  ``headers`` is a LIST OF SINGLE-KEY MAPS, not a plain map.

    This file is ALWAYS written, with zeroes when nothing was found -- a sample that
    legitimately detects zero arrays must appear as an explicit zero row, not vanish.

    The section id is ``mcaat_genstats``.  It is a General Statistics contribution, not a
    report section of its own, so it is deliberately absent from ``custom_content.order``
    in ``assets/multiqc_config.yml`` and is configured there by the
    ``table_columns_visible`` / ``table_columns_placement`` blocks instead.  The five
    cohort section ids live in ``bin/mcaat_multiqc_sections.py``.
    """
    header_block = """\
# id: 'mcaat_genstats'
# section_name: 'MCAAT CRISPR arrays'
# description: 'CRISPR arrays and spacers detected directly in unassembled reads by MCAAT.'
# plot_type: 'generalstats'
# headers:
#     - CRISPR_arrays:
#         title: 'CRISPR arrays'
#         description: 'CRISPR arrays detected in unassembled reads by MCAAT'
#         min: 0
#         format: '{:,.0f}'
#         scale: 'Blues'
#     - CRISPR_spacers:
#         title: 'Spacers'
#         description: 'Total spacers across all detected arrays'
#         min: 0
#         format: '{:,.0f}'
#         scale: 'Greens'
#     - CRISPR_spacers_unique:
#         title: 'Unique spacers'
#         description: 'Distinct spacer sequences within this sample'
#         min: 0
#         format: '{:,.0f}'
#         scale: 'Purples'
#         hidden: True
#     - CRISPR_median_spacer_len:
#         title: 'Median spacer bp'
#         description: 'Median spacer length in base pairs'
#         min: 0
#         format: '{:,.1f}'
#         scale: 'Oranges'
#         hidden: True
#     - CRISPR_median_repeat_len:
#         title: 'Median repeat bp'
#         description: 'Median consensus repeat length in base pairs'
#         min: 0
#         format: '{:,.1f}'
#         scale: 'Reds'
#         hidden: True
"""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(header_block)
        handle.write(
            "Sample\tCRISPR_arrays\tCRISPR_spacers\tCRISPR_spacers_unique\t"
            "CRISPR_median_spacer_len\tCRISPR_median_repeat_len\n"
        )
        handle.write(
            "{0}\t{1}\t{2}\t{3}\t{4:.1f}\t{5:.1f}\n".format(
                sanitise_field(sample),
                n_arrays,
                n_spacers,
                n_spacers_unique,
                median_spacer_len,
                median_repeat_len,
            )
        )


def write_provenance_tsv(path: str, row: Dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(PROVENANCE_TSV_COLUMNS) + "\n")
        handle.write("\t".join(sanitise_field(row.get(c, "NA")) for c in PROVENANCE_TSV_COLUMNS) + "\n")


def safe_median(values: Iterable[int]) -> float:
    materialised = list(values)
    if not materialised:
        return 0.0
    return float(median(materialised))


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mcaat_parse_arrays.py",
        description=(
            "Parse MCAAT CRISPR_Arrays_*.txt into a per-sample TSV, a per-sample spacer "
            "FASTA and a MultiQC general-statistics file, expanding dash-encoded repeat "
            "variants back to the consensus repeat sequence."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Exit codes:\n"
            "  0  success (including the legitimate zero-array result)\n"
            "  1  parameters.json readback mismatch, or an unwritable output path\n"
        ),
    )
    parser.add_argument("--sample", required=True, help="Sample identifier (meta.id).")
    parser.add_argument(
        "--arrays",
        required=True,
        nargs="+",
        metavar="PATH",
        help="CRISPR_Arrays_*.txt files, and/or directories to scan for them.",
    )
    parser.add_argument(
        "--parameters",
        metavar="FILE",
        help="MCAAT parameters.json for the provenance echo and the readback assertion.",
    )
    parser.add_argument(
        "--expected",
        metavar="FILE",
        help=(
            "JSON file of the parameter values this pipeline intended. Compared against "
            "parameters.json; a mismatch is fatal unless --lenient-params is given. This "
            "is the only detector for a flag MCAAT silently accepted and dropped."
        ),
    )
    parser.add_argument(
        "--lenient-params",
        action="store_true",
        help="Downgrade a parameters.json readback mismatch from fatal to a warning.",
    )
    parser.add_argument("--out-tsv", required=True, metavar="FILE", help="Per-spacer TSV output.")
    parser.add_argument(
        "--out-fasta", required=True, metavar="FILE", help="Per-sample spacer FASTA output."
    )
    parser.add_argument(
        "--out-mqc", required=True, metavar="FILE", help="MultiQC general-stats *_mqc.tsv output."
    )
    parser.add_argument(
        "--out-provenance",
        metavar="FILE",
        help="Optional per-sample provenance TSV echoing the resolved MCAAT parameters.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    sample = args.sample

    files = discover_array_files(args.arrays)
    if not files:
        # R10: a missing file is a legitimate zero-result, not an error.
        warn(
            "no CRISPR_Arrays_*.txt found for sample {0}; emitting an explicit "
            "zero-result".format(sample)
        )

    blocks: List[ArrayBlock] = []
    declared_arrays_seen: List[int] = []
    declared_spacers_seen: List[int] = []
    for path in files:
        info("reading {0}".format(path))
        file_blocks, header = parse_array_file(path)
        blocks.extend(file_blocks)
        if header.declared_arrays is not None:
            declared_arrays_seen.append(header.declared_arrays)
        if header.declared_spacers is not None:
            declared_spacers_seen.append(header.declared_spacers)

    n_arrays = len(blocks)
    spacer_lengths = [len(entry.spacer) for block in blocks for entry in block.entries]
    repeat_lengths = [len(block.consensus) for block in blocks if block.consensus]
    n_spacers = len(spacer_lengths)
    unique_spacers = set(entry.spacer for block in blocks for entry in block.entries)

    # R2: cross-check only.  These are global run totals repeated in every split file --
    # summing them would double-count.
    if declared_arrays_seen:
        declared_arrays = max(declared_arrays_seen)
        if len(set(declared_arrays_seen)) > 1:
            warn(
                "the '# Arrays' preamble differs between split files ({0}); MCAAT repeats a "
                "single global total, so this indicates the files come from different "
                "runs".format(sorted(set(declared_arrays_seen)))
            )
        if declared_arrays != n_arrays:
            warn(
                "preamble declares {0} arrays but {1} array blocks were parsed".format(
                    declared_arrays, n_arrays
                )
            )
    if declared_spacers_seen:
        declared_spacers = max(declared_spacers_seen)
        if declared_spacers != n_spacers:
            warn(
                "preamble declares {0} spacers but {1} spacer rows were parsed".format(
                    declared_spacers, n_spacers
                )
            )

    parameters: Optional[dict] = None
    if args.parameters:
        if os.path.isfile(args.parameters):
            parameters = load_json_defensively(args.parameters)
        else:
            warn("--parameters {0} does not exist".format(args.parameters))

    exit_code = 0
    if args.expected:
        expected = load_json_defensively(args.expected) or {}
        if not expected:
            warn("--expected produced no comparable values; readback skipped")
        elif parameters is None:
            message = (
                "cannot perform the parameters.json readback: parameters.json is missing or "
                "unparseable"
            )
            if args.lenient_params:
                warn(message)
            else:
                die(message)
        else:
            mismatches = readback_parameters(parameters, expected)
            if mismatches:
                banner = (
                    "parameters.json readback FAILED for sample {0}. MCAAT's argument parse "
                    "loop has no unknown-argument branch, so a flag it did not understand is "
                    "silently ignored and the run proceeds with defaults. Mismatches:".format(
                        sample
                    )
                )
                sys.stderr.write("MCAAT_PARSE_ARRAYS: ERROR: {0}\n".format(banner))
                for item in mismatches:
                    sys.stderr.write("MCAAT_PARSE_ARRAYS:   - {0}\n".format(item))
                if args.lenient_params:
                    warn("--lenient-params set; continuing despite the mismatch")
                else:
                    exit_code = 1
            else:
                info("parameters.json readback OK ({0} keys)".format(len(expected)))

    write_arrays_tsv(args.out_tsv, sample, blocks)
    write_spacers_fasta(args.out_fasta, sample, blocks)
    write_genstats_mqc(
        args.out_mqc,
        sample,
        n_arrays,
        n_spacers,
        len(unique_spacers),
        safe_median(spacer_lengths),
        safe_median(repeat_lengths),
    )
    if args.out_provenance:
        write_provenance_tsv(args.out_provenance, provenance_row(sample, parameters))

    info(
        "sample {0}: {1} arrays, {2} spacers ({3} unique) from {4} file(s)".format(
            sample, n_arrays, n_spacers, len(unique_spacers), len(files)
        )
    )

    if exit_code != 0:
        sys.stderr.write(
            "MCAAT_PARSE_ARRAYS: output files were written before failing so the mismatch "
            "can be inspected.\n"
        )
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
