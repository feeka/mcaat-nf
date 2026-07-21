#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mcaat_aggregate.py — cohort aggregation for the mcaat-nf Nextflow pipeline.

Consumes the per-sample tables produced by ``mcaat_parse_arrays.py`` and emits the
cross-sample deliverables:

  ``cohort_arrays.tsv``             every per-sample row, one header
  ``cohort_summary.tsv``            one row per sample
  ``cohort_spacer_sharing.tsv``     long-format pairwise sample x sample sharing
  ``cohort_spacer_redundancy.tsv``  one row per spacer sequence seen in >= N samples
  ``cohort_repeat_families.tsv``    one row per distinct consensus repeat
  ``cohort_spacers.unique.fasta``   one record per distinct spacer sequence

Design notes
------------
* Standard library only.  No pandas, no numpy.  The pipeline runs this in a plain
  ``python:3.12`` container.
* ``cohort_arrays.tsv`` is written by STREAMING the inputs, so the concatenation itself
  never holds the cohort in memory.  Only the per-sample summaries and the
  spacer -> samples / repeat -> samples indexes are retained.
* Pairwise sharing is computed from the spacer -> samples index rather than by
  intersecting every pair of sets: for a spacer present in k samples we add one to each
  of the k*(k-1)/2 pairs.  On real cohorts sharing is sparse, so this is dramatically
  cheaper than O(samples^2) set intersections and produces identical numbers.
* Every sample that reached this stage appears in ``cohort_summary.tsv``, including
  samples with zero arrays.  A legitimate negative result must be an explicit zero row,
  never a missing row.  The sample roster is therefore established from THREE independent
  sources, in decreasing order of reliability:

    1. ``--samples`` -- an explicit roster, normally derived by the caller from the staged
       filenames.  This is the only source that still works when the caller supplies one
       pre-concatenated cohort table, because a zero-array sample contributes no rows to
       it and its identity is gone by the time this script runs.
    2. the ``sample`` column of the ``*.provenance.tsv`` tables, which carry exactly one
       row per sample regardless of how many arrays were found.
    3. the filename of a per-sample ``<sample>.arrays.tsv`` that turned out to have no
       data rows.

* Inputs are classified by FILENAME SUFFIX, not by which option they arrived on.
  ``*.provenance.tsv`` is provenance, any other ``*.tsv`` is an arrays table -- wherever
  it was passed.  A mis-ordered call site therefore produces correct tables instead of
  silently-``NA`` provenance columns.
* An empty cohort is NOT an error.  Every table is written with its header row even when
  no sample contributed a single array, so downstream consumers (and the module's output
  declarations) always find the files they expect.
"""

from __future__ import annotations

import argparse
import os
import sys
from itertools import combinations
from statistics import median
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

# --------------------------------------------------------------------------------------
# Column contracts
# --------------------------------------------------------------------------------------

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

SUMMARY_COLUMNS: Tuple[str, ...] = (
    "sample",
    "n_arrays",
    "n_spacers",
    "n_spacers_unique",
    "n_consensus_repeats",
    "median_repeat_len",
    "median_spacer_len",
    "min_spacer_len",
    "max_spacer_len",
    "n_spacers_shared",
    "pct_spacers_shared",
    "n_repeats_shared",
    "mcaat_threads",
    "mcaat_ram_gb",
    "mcaat_threshold_multiplicity",
)

SHARING_COLUMNS: Tuple[str, ...] = (
    "sample_a",
    "sample_b",
    "n_shared_spacers",
    "jaccard",
    "n_shared_repeats",
)

REDUNDANCY_COLUMNS: Tuple[str, ...] = (
    "spacer",
    "spacer_len",
    "n_samples",
    "samples",
    "n_occurrences",
    "n_consensus_repeats",
    "consensus_repeats",
)

FAMILY_COLUMNS: Tuple[str, ...] = (
    "consensus_repeat",
    "repeat_len",
    "n_samples",
    "samples",
    "n_arrays_total",
    "n_spacers_total",
    "n_spacers_unique",
    "median_spacer_len",
)


def warn(msg: str) -> None:
    sys.stderr.write("MCAAT_AGGREGATE: WARNING: {0}\n".format(msg))


def info(msg: str) -> None:
    sys.stderr.write("MCAAT_AGGREGATE: {0}\n".format(msg))


def die(msg: str, code: int = 1) -> None:
    sys.stderr.write("MCAAT_AGGREGATE: ERROR: {0}\n".format(msg))
    sys.exit(code)


def sanitise_field(value: str) -> str:
    return value.replace("\t", " ").replace("\r", " ").replace("\n", " ")


# --------------------------------------------------------------------------------------
# Input discovery
# --------------------------------------------------------------------------------------


def collect_inputs(
    explicit: Sequence[str], directory: Optional[str], suffixes: Sequence[str]
) -> List[str]:
    """Collect files from an explicit list and/or by scanning a directory."""
    found: List[str] = []
    for item in explicit or []:
        if os.path.isfile(item):
            found.append(item)
        elif os.path.isdir(item):
            found.extend(_scan_dir(item, suffixes))
        else:
            warn("input path does not exist and is ignored: {0}".format(item))
    if directory:
        if os.path.isdir(directory):
            found.extend(_scan_dir(directory, suffixes))
        else:
            info("directory {0} is absent; treated as empty".format(directory))
    # Deterministic, de-duplicated ordering.
    return sorted(set(os.path.abspath(p) for p in found))


def classify_inputs(paths: Iterable[str]) -> Tuple[List[str], List[str]]:
    """Split a flat list of TSVs into (arrays, provenance) by filename suffix.

    Classification is deliberately independent of the option a file arrived on.  The
    module contract puts arrays first and provenance second, but a call site that swapped
    the two slots would otherwise produce a run that succeeds while writing `NA` into
    every provenance column.  Suffixes are unambiguous, so honouring them is strictly
    safer than honouring position.
    """
    arrays: List[str] = []
    provenance: List[str] = []
    for path in paths:
        name = os.path.basename(path)
        if name.endswith(".provenance.tsv"):
            provenance.append(path)
        elif name.endswith(".tsv"):
            arrays.append(path)
        else:
            warn("{0} has no recognised .tsv suffix and is ignored".format(name))
    return arrays, provenance


def split_sample_roster(values: Sequence[str]) -> List[str]:
    """Expand a roster given as repeated tokens and/or comma-separated strings."""
    roster: List[str] = []
    for value in values or []:
        for item in str(value).split(","):
            item = item.strip()
            if item and item not in roster:
                roster.append(item)
    return roster


def _scan_dir(directory: str, suffixes: Sequence[str]) -> List[str]:
    out: List[str] = []
    for name in sorted(os.listdir(directory)):
        full = os.path.join(directory, name)
        if os.path.isfile(full) and any(name.endswith(s) for s in suffixes):
            out.append(full)
    return out


def read_tsv(path: str) -> Iterator[Dict[str, str]]:
    """Yield rows of a TSV as dicts.  Tolerates ragged rows with a warning."""
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        header_line = handle.readline()
        if not header_line:
            warn("{0} is empty".format(path))
            return
        header = header_line.rstrip("\r\n").split("\t")
        width = len(header)
        for lineno, raw in enumerate(handle, start=2):
            line = raw.rstrip("\r\n")
            if line == "":
                continue
            fields = line.split("\t")
            if len(fields) != width:
                warn(
                    "{0}:{1}: expected {2} fields, found {3}; row skipped".format(
                        os.path.basename(path), lineno, width, len(fields)
                    )
                )
                continue
            yield dict(zip(header, fields))


def to_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------------------------------
# Accumulators
# --------------------------------------------------------------------------------------


class SampleStats(object):
    __slots__ = (
        "sample",
        "arrays",
        "n_spacer_rows",
        "unique_spacers",
        "consensus_repeats",
        "spacer_lengths",
        "repeat_lengths",
    )

    def __init__(self, sample: str) -> None:
        self.sample = sample
        self.arrays: Set[int] = set()
        self.n_spacer_rows = 0
        self.unique_spacers: Set[str] = set()
        self.consensus_repeats: Set[str] = set()
        self.spacer_lengths: List[int] = []
        self.repeat_lengths: List[int] = []


class RepeatFamily(object):
    __slots__ = ("repeat", "samples", "arrays", "n_spacers", "unique_spacers", "spacer_lengths")

    def __init__(self, repeat: str) -> None:
        self.repeat = repeat
        self.samples: Set[str] = set()
        self.arrays: Set[Tuple[str, int]] = set()
        self.n_spacers = 0
        self.unique_spacers: Set[str] = set()
        self.spacer_lengths: List[int] = []


class SpacerRecord(object):
    __slots__ = ("spacer", "samples", "occurrences", "repeats")

    def __init__(self, spacer: str) -> None:
        self.spacer = spacer
        self.samples: Set[str] = set()
        self.occurrences = 0
        self.repeats: Set[str] = set()


def safe_median(values: Sequence[int]) -> float:
    if not values:
        return 0.0
    return float(median(values))


# --------------------------------------------------------------------------------------
# Core aggregation
# --------------------------------------------------------------------------------------


def aggregate(
    array_tsvs: Sequence[str],
    cohort_tsv_out: str,
) -> Tuple[Dict[str, SampleStats], Dict[str, SpacerRecord], Dict[str, RepeatFamily]]:
    """Stream every per-sample TSV into ``cohort_tsv_out`` while building the indexes."""
    samples: Dict[str, SampleStats] = {}
    spacers: Dict[str, SpacerRecord] = {}
    families: Dict[str, RepeatFamily] = {}

    with open(cohort_tsv_out, "w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(ARRAYS_TSV_COLUMNS) + "\n")

        for path in array_tsvs:
            seen_any = False
            # A zero-array sample yields a header-only TSV, so the sample id can only be
            # recovered from the filename.  Strip the known suffix rather than splitting on
            # '.', because the samplesheet permits dots inside a sample identifier.
            #
            # This recovery only fires for PER-SAMPLE tables.  When the caller supplies a
            # single pre-concatenated cohort table instead, a zero-array sample leaves no
            # trace in it at all and can only be recovered from `--samples` or from the
            # provenance tables -- both of which are applied by main().
            basename = os.path.basename(path)
            for suffix in (".arrays.tsv", ".tsv"):
                if basename.endswith(suffix):
                    sample_hint = basename[: -len(suffix)]
                    break
            else:
                sample_hint = basename
            for row in read_tsv(path):
                missing = [c for c in ARRAYS_TSV_COLUMNS if c not in row]
                if missing:
                    warn(
                        "{0}: missing column(s) {1}; file skipped".format(
                            os.path.basename(path), ",".join(missing)
                        )
                    )
                    break
                seen_any = True

                sample = row["sample"]
                array_id = to_int(row["array_id"], -1)
                consensus = row["consensus_repeat"]
                spacer = row["spacer"]
                spacer_len = to_int(row["spacer_len"], len(spacer))
                repeat_len = to_int(row["consensus_repeat_len"], len(consensus))

                out.write("\t".join(sanitise_field(row[c]) for c in ARRAYS_TSV_COLUMNS) + "\n")

                stats = samples.get(sample)
                if stats is None:
                    stats = SampleStats(sample)
                    samples[sample] = stats
                stats.arrays.add(array_id)
                stats.n_spacer_rows += 1
                stats.unique_spacers.add(spacer)
                stats.spacer_lengths.append(spacer_len)
                if consensus:
                    stats.consensus_repeats.add(consensus)
                    stats.repeat_lengths.append(repeat_len)

                record = spacers.get(spacer)
                if record is None:
                    record = SpacerRecord(spacer)
                    spacers[spacer] = record
                record.samples.add(sample)
                record.occurrences += 1
                if consensus:
                    record.repeats.add(consensus)

                if consensus:
                    family = families.get(consensus)
                    if family is None:
                        family = RepeatFamily(consensus)
                        families[consensus] = family
                    family.samples.add(sample)
                    family.arrays.add((sample, array_id))
                    family.n_spacers += 1
                    family.unique_spacers.add(spacer)
                    family.spacer_lengths.append(spacer_len)

            if not seen_any:
                # A zero-array sample produces a header-only TSV.  It MUST still appear
                # in the cohort summary as an explicit zero row.
                if sample_hint and sample_hint not in samples:
                    samples[sample_hint] = SampleStats(sample_hint)
                    info(
                        "{0} contains no data rows; recording {1} as a zero-array "
                        "sample".format(os.path.basename(path), sample_hint)
                    )

    return samples, spacers, families


def pairwise_counts(
    index: Dict[str, Set[str]],
) -> Dict[Tuple[str, str], int]:
    """Count shared keys per unordered sample pair from a key -> samples index."""
    pairs: Dict[Tuple[str, str], int] = {}
    for members in index.values():
        if len(members) < 2:
            continue
        for a, b in combinations(sorted(members), 2):
            pairs[(a, b)] = pairs.get((a, b), 0) + 1
    return pairs


# --------------------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------------------


def write_summary(
    path: str,
    samples: Dict[str, SampleStats],
    spacers: Dict[str, SpacerRecord],
    families: Dict[str, RepeatFamily],
    provenance: Dict[str, Dict[str, str]],
) -> None:
    shared_spacers_per_sample: Dict[str, int] = dict((s, 0) for s in samples)
    for record in spacers.values():
        if len(record.samples) < 2:
            continue
        for sample in record.samples:
            if sample in shared_spacers_per_sample:
                shared_spacers_per_sample[sample] += 1

    shared_repeats_per_sample: Dict[str, int] = dict((s, 0) for s in samples)
    for family in families.values():
        if len(family.samples) < 2:
            continue
        for sample in family.samples:
            if sample in shared_repeats_per_sample:
                shared_repeats_per_sample[sample] += 1

    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(SUMMARY_COLUMNS) + "\n")
        for sample in sorted(samples):
            stats = samples[sample]
            n_unique = len(stats.unique_spacers)
            n_shared = shared_spacers_per_sample.get(sample, 0)
            pct_shared = (100.0 * n_shared / n_unique) if n_unique else 0.0
            prov = provenance.get(sample, {})
            handle.write(
                "\t".join(
                    (
                        sanitise_field(sample),
                        str(len(stats.arrays)),
                        str(stats.n_spacer_rows),
                        str(n_unique),
                        str(len(stats.consensus_repeats)),
                        "{0:.1f}".format(safe_median(stats.repeat_lengths)),
                        "{0:.1f}".format(safe_median(stats.spacer_lengths)),
                        str(min(stats.spacer_lengths) if stats.spacer_lengths else 0),
                        str(max(stats.spacer_lengths) if stats.spacer_lengths else 0),
                        str(n_shared),
                        "{0:.2f}".format(pct_shared),
                        str(shared_repeats_per_sample.get(sample, 0)),
                        sanitise_field(prov.get("mcaat_threads", "NA")),
                        sanitise_field(prov.get("mcaat_ram_gb", "NA")),
                        sanitise_field(prov.get("mcaat_threshold_multiplicity", "NA")),
                    )
                )
                + "\n"
            )


def write_sharing(
    path: str,
    samples: Dict[str, SampleStats],
    spacers: Dict[str, SpacerRecord],
    families: Dict[str, RepeatFamily],
    min_sharing: int,
) -> None:
    spacer_index = dict((k, v.samples) for k, v in spacers.items())
    repeat_index = dict((k, v.samples) for k, v in families.items())
    shared_spacer_pairs = pairwise_counts(spacer_index)
    shared_repeat_pairs = pairwise_counts(repeat_index)

    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(SHARING_COLUMNS) + "\n")
        for (a, b) in sorted(
            set(list(shared_spacer_pairs.keys()) + list(shared_repeat_pairs.keys()))
        ):
            n_shared = shared_spacer_pairs.get((a, b), 0)
            if n_shared < min_sharing:
                continue
            size_a = len(samples[a].unique_spacers) if a in samples else 0
            size_b = len(samples[b].unique_spacers) if b in samples else 0
            union = size_a + size_b - n_shared
            jaccard = (float(n_shared) / union) if union > 0 else 0.0
            handle.write(
                "\t".join(
                    (
                        sanitise_field(a),
                        sanitise_field(b),
                        str(n_shared),
                        "{0:.6f}".format(jaccard),
                        str(shared_repeat_pairs.get((a, b), 0)),
                    )
                )
                + "\n"
            )


def write_redundancy(path: str, spacers: Dict[str, SpacerRecord], min_samples: int) -> None:
    """Cross-sample spacer redundancy: which exact spacer sequences recur, and where."""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(REDUNDANCY_COLUMNS) + "\n")
        rows = [r for r in spacers.values() if len(r.samples) >= min_samples]
        rows.sort(key=lambda r: (-len(r.samples), -r.occurrences, r.spacer))
        for record in rows:
            handle.write(
                "\t".join(
                    (
                        sanitise_field(record.spacer),
                        str(len(record.spacer)),
                        str(len(record.samples)),
                        sanitise_field(";".join(sorted(record.samples))),
                        str(record.occurrences),
                        str(len(record.repeats)),
                        sanitise_field(";".join(sorted(record.repeats))),
                    )
                )
                + "\n"
            )


def write_families(path: str, families: Dict[str, RepeatFamily]) -> None:
    """Write the repeat-family table, ranked by cohort breadth then spacer count."""
    ordered = sorted(
        families.values(),
        key=lambda f: (-len(f.samples), -f.n_spacers, f.repeat),
    )
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(FAMILY_COLUMNS) + "\n")
        for family in ordered:
            handle.write(
                "\t".join(
                    (
                        sanitise_field(family.repeat),
                        str(len(family.repeat)),
                        str(len(family.samples)),
                        sanitise_field(";".join(sorted(family.samples))),
                        str(len(family.arrays)),
                        str(family.n_spacers),
                        str(len(family.unique_spacers)),
                        "{0:.1f}".format(safe_median(family.spacer_lengths)),
                    )
                )
                + "\n"
            )


def write_unique_fasta(path: str, spacers: Dict[str, SpacerRecord]) -> None:
    """One record per distinct spacer, in deterministic (sorted) order.

    Identifiers are ``cohort_spacer_<i>`` with the provenance in the description, so the
    file is byte-stable for a fixed cohort despite MCAAT's non-deterministic internal
    ordering.
    """
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for index, spacer in enumerate(sorted(spacers), start=1):
            record = spacers[spacer]
            handle.write(
                ">cohort_spacer_{index} len={length} n_samples={ns} n_occurrences={no} "
                "samples={samples}\n".format(
                    index=index,
                    length=len(spacer),
                    ns=len(record.samples),
                    no=record.occurrences,
                    samples=",".join(sorted(record.samples)),
                )
            )
            handle.write(spacer + "\n")


def load_provenance(paths: Sequence[str]) -> Dict[str, Dict[str, str]]:
    provenance: Dict[str, Dict[str, str]] = {}
    for path in paths:
        for row in read_tsv(path):
            sample = row.get("sample")
            if sample:
                provenance[sample] = row
    return provenance


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mcaat_aggregate.py",
        description=(
            "Aggregate per-sample MCAAT array tables into cohort deliverables: a merged "
            "TSV, a per-sample summary, a pairwise spacer-sharing table, a cross-sample "
            "spacer redundancy table, a repeat-family table and a unique-spacer FASTA."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--arrays-tsv",
        nargs="*",
        default=[],
        metavar="FILE",
        help="Per-sample <sample>.arrays.tsv files (or directories containing them).",
    )
    parser.add_argument(
        "--arrays-dir",
        metavar="DIR",
        help="Directory to scan for *.tsv (absent directory is treated as empty).",
    )
    parser.add_argument(
        "--provenance-tsv",
        nargs="*",
        default=[],
        metavar="FILE",
        help="Per-sample <sample>.provenance.tsv files echoing MCAAT's parameters.json.",
    )
    parser.add_argument("--provenance-dir", metavar="DIR", help="Directory to scan for *.tsv.")
    parser.add_argument(
        "--samples",
        nargs="*",
        default=[],
        metavar="SAMPLE",
        help="Explicit sample roster, repeated and/or comma-separated. Every listed sample "
        "is guaranteed a row in cohort_summary.tsv even if it contributed no data rows. "
        "Required for correct zero-array handling when the arrays input is a single "
        "pre-concatenated cohort table rather than one file per sample.",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        metavar="DIR",
        help="Directory for the cohort outputs (default: .).",
    )
    parser.add_argument(
        "--min-spacer-sharing",
        type=int,
        default=2,
        metavar="N",
        help="Minimum shared spacers for a sample pair to enter the sharing table (default: 2).",
    )
    parser.add_argument(
        "--min-spacer-redundancy",
        type=int,
        default=2,
        metavar="N",
        help="Minimum number of samples a spacer must appear in to enter the redundancy "
        "table (default: 2).",
    )
    parser.add_argument(
        "--no-unique-fasta",
        action="store_true",
        help="Do not write cohort_spacers.unique.fasta. Safe to set from ext.args: the "
        "COHORT_SUMMARISE module declares that output `optional: true` precisely so this "
        "flag cannot turn into a 'Missing output file' failure after the work is done.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    outdir = args.outdir
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    path = lambda name: os.path.join(outdir, name)  # noqa: E731

    # Gather every candidate from every option into ONE list, then classify by filename
    # suffix.  See classify_inputs(): position is the contract, but suffix is the truth,
    # so a call site that swaps the two slots still produces correct tables.
    candidates = collect_inputs(
        list(args.arrays_tsv) + list(args.provenance_tsv),
        None,
        (".tsv",),
    )
    for directory in (args.arrays_dir, args.provenance_dir):
        if directory:
            candidates.extend(collect_inputs([], directory, (".tsv",)))

    # Never read the file we are about to overwrite (possible when --outdir points at a
    # directory that is also being scanned).
    cohort_arrays_out = os.path.abspath(path("cohort_arrays.tsv"))
    candidates = sorted(set(p for p in candidates if p != cohort_arrays_out))

    array_tsvs, provenance_tsvs = classify_inputs(candidates)

    if not array_tsvs:
        # An empty cohort is a legitimate (if unusual) result, and the caller's output
        # declarations still expect every table to exist.  Warn loudly, then write the
        # header-only tables rather than aborting after the per-sample work succeeded.
        warn(
            "no arrays TSV found (checked --arrays-tsv and --arrays-dir); writing "
            "header-only cohort tables"
        )
    info(
        "aggregating {0} arrays table(s) and {1} provenance table(s)".format(
            len(array_tsvs), len(provenance_tsvs)
        )
    )

    samples, spacers, families = aggregate(array_tsvs, path("cohort_arrays.tsv"))
    provenance = load_provenance(provenance_tsvs)

    # Complete the sample roster.  A sample that legitimately found zero arrays must be an
    # explicit zero row in cohort_summary.tsv, never a missing row -- and when the arrays
    # input was pre-concatenated, these two sources are the only surviving evidence that
    # the sample existed at all.
    for sample in split_sample_roster(args.samples) + sorted(provenance):
        if sample and sample not in samples:
            samples[sample] = SampleStats(sample)
            info("{0} contributed no array rows; recorded as a zero-array sample".format(sample))

    write_summary(path("cohort_summary.tsv"), samples, spacers, families, provenance)
    write_sharing(
        path("cohort_spacer_sharing.tsv"), samples, spacers, families, args.min_spacer_sharing
    )
    write_redundancy(path("cohort_spacer_redundancy.tsv"), spacers, args.min_spacer_redundancy)
    write_families(path("cohort_repeat_families.tsv"), families)

    if not args.no_unique_fasta:
        write_unique_fasta(path("cohort_spacers.unique.fasta"), spacers)

    total_spacer_rows = sum(s.n_spacer_rows for s in samples.values())
    info(
        "cohort: {0} samples, {1} spacer rows, {2} unique spacers, {3} repeat "
        "families".format(len(samples), total_spacer_rows, len(spacers), len(families))
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
