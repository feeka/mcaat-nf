#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mcaat_multiqc_sections.py — render MCAAT cohort tables as MultiQC custom content.

Reads the cohort tables produced by ``mcaat_aggregate.py`` and writes MultiQC
custom-content files using the ZERO-CONFIG MAGIC FILENAME mechanism: any file matching
``*_mqc.(tsv|yaml|...)`` is discovered by MultiQC automatically.  No ``sp:`` search
pattern and no ``multiqc_config.yml`` edit is required, so the report still renders
correctly when a user supplies their own ``--multiqc_config``.

Files emitted (all into ``--outdir``), with the section id each one carries:

  ``mcaat_cohort_summary_mqc.yaml``   id ``mcaat_cohort_summary``   table
  ``mcaat_arrays_bargraph_mqc.tsv``   id ``mcaat_arrays_bargraph``  bargraph
  ``mcaat_spacerlen_mqc.tsv``         id ``mcaat_spacerlen``        linegraph
  ``mcaat_repeatlen_mqc.tsv``         id ``mcaat_repeatlen``        linegraph
  ``mcaat_sharing_mqc.tsv``           id ``mcaat_sharing``          heatmap

Those five ids are the complete, canonical set and they are the same five listed under
``custom_content.order`` in ``assets/multiqc_config.yml``.  The two lists must be kept in
lock-step: an id ordered but never emitted is silently ignored, and an id emitted but
never ordered floats to an arbitrary position in the report.  Adding or renaming a
section here therefore requires the same edit there.

``mcaat_genstats`` is NOT emitted here -- it is the per-sample General Statistics
contribution written by ``mcaat_parse_arrays.py`` as ``<sample>_mcaat_genstats_mqc.tsv``.

The length-distribution sections are skipped when no cohort arrays table was supplied,
and the sharing heatmap is skipped (with a note on stderr) when the cohort has fewer than
two samples, because MultiQC renders a degenerate single-cell plot otherwise.  The
summary table is ALWAYS written, even for an empty cohort, so that the caller's
``*_mqc.{tsv,yaml}`` output declaration always matches at least one file.

Standard library only.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple


def warn(msg: str) -> None:
    sys.stderr.write("MCAAT_MULTIQC: WARNING: {0}\n".format(msg))


def info(msg: str) -> None:
    sys.stderr.write("MCAAT_MULTIQC: {0}\n".format(msg))


def read_tsv(path: str) -> Iterator[Dict[str, str]]:
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        header_line = handle.readline()
        if not header_line:
            return
        header = header_line.rstrip("\r\n").split("\t")
        width = len(header)
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "":
                continue
            fields = line.split("\t")
            if len(fields) != width:
                continue
            yield dict(zip(header, fields))


def to_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def to_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def yaml_scalar(value: object) -> str:
    """Emit a YAML scalar.  Strings are always double-quoted and escaped."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return repr(round(value, 6))
    text = str(value)
    text = text.replace("\\", "\\\\").replace('"', '\\"')
    text = text.replace("\n", " ").replace("\r", " ").replace("\t", " ")
    return '"{0}"'.format(text)


def yaml_block(data: object, indent: int = 0) -> List[str]:
    """Serialise nested dicts/lists/scalars to YAML lines.  Deterministic ordering.

    Mapping keys are emitted through ``yaml_scalar`` exactly like values, i.e. always
    double-quoted and escaped.  This is not cosmetic: the keys of the ``data:`` mapping in
    the cohort summary section are sample identifiers taken verbatim from the samplesheet,
    and ``assets/schema_input.json`` only constrains ``sample`` to ``^\\S+$``.  Left bare,
    YAML 1.1 resolution would silently retype them -- ``007`` -> int 7, ``1.10`` -> float
    1.1, ``No`` -> bool false -- and a sample called ``#weird`` would open a comment and
    disappear from the document while it still parsed cleanly.  MultiQC would then label
    the row with the retyped value and fail to join it to the General Statistics row, whose
    id comes from the TSV as a literal string.  Quoting forces every key to stay a string.
    """
    pad = " " * indent
    lines: List[str] = []
    if isinstance(data, dict):
        for key, value in data.items():
            if isinstance(value, (dict, list)) and value:
                lines.append("{0}{1}:".format(pad, yaml_scalar(key)))
                lines.extend(yaml_block(value, indent + 4))
            elif isinstance(value, (dict, list)):
                lines.append(
                    "{0}{1}: {2}".format(
                        pad, yaml_scalar(key), "{}" if isinstance(value, dict) else "[]"
                    )
                )
            else:
                lines.append("{0}{1}: {2}".format(pad, yaml_scalar(key), yaml_scalar(value)))
    elif isinstance(data, list):
        for item in data:
            if isinstance(item, (dict, list)):
                nested = yaml_block(item, indent + 2)
                if nested:
                    first = nested[0].lstrip()
                    lines.append("{0}- {1}".format(pad, first))
                    lines.extend(nested[1:])
            else:
                lines.append("{0}- {1}".format(pad, yaml_scalar(item)))
    else:
        lines.append("{0}{1}".format(pad, yaml_scalar(data)))
    return lines


def write_yaml(path: str, document: Dict[str, object]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(yaml_block(document)) + "\n")


def write_mqc_tsv(
    path: str,
    meta: Sequence[Tuple[str, str]],
    pconfig: Sequence[Tuple[str, str]],
    header_row: Sequence[str],
    rows: Sequence[Sequence[str]],
) -> None:
    """Write a MultiQC custom-content TSV: ``#``-prefixed YAML front matter, then a table."""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for key, value in meta:
            handle.write("# {0}: {1}\n".format(key, value))
        if pconfig:
            handle.write("# pconfig:\n")
            for key, value in pconfig:
                handle.write("#     {0}: {1}\n".format(key, value))
        handle.write("\t".join(header_row) + "\n")
        for row in rows:
            handle.write("\t".join(str(cell) for cell in row) + "\n")


# --------------------------------------------------------------------------------------
# Sections
# --------------------------------------------------------------------------------------


def section_summary_table(path: str, summary_rows: Sequence[Dict[str, str]]) -> None:
    headers = {
        "n_arrays": {
            "title": "Arrays",
            "description": "Distinct CRISPR arrays detected by MCAAT",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Blues",
        },
        "n_spacers": {
            "title": "Spacers",
            "description": "Total spacer occurrences across all arrays",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Greens",
        },
        "n_spacers_unique": {
            "title": "Unique spacers",
            "description": "Distinct spacer sequences within the sample",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Purples",
        },
        "n_consensus_repeats": {
            "title": "Repeat families",
            "description": "Distinct consensus repeat sequences",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "YlGn",
        },
        "median_repeat_len": {
            "title": "Median repeat (bp)",
            "description": "Median consensus repeat length",
            "min": 0,
            "format": "{:,.1f}",
            "scale": "Reds",
        },
        "median_spacer_len": {
            "title": "Median spacer (bp)",
            "description": "Median spacer length",
            "min": 0,
            "format": "{:,.1f}",
            "scale": "Oranges",
        },
        "n_spacers_shared": {
            "title": "Shared spacers",
            "description": "Unique spacers also present in at least one other sample",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "BuPu",
        },
        "pct_spacers_shared": {
            "title": "% shared",
            "description": "Percentage of this sample's unique spacers seen elsewhere in the cohort",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "RdYlBu-rev",
        },
    }

    data: Dict[str, Dict[str, object]] = {}
    for row in summary_rows:
        sample = row.get("sample", "")
        if not sample:
            continue
        data[sample] = {
            "n_arrays": to_int(row.get("n_arrays", "0")),
            "n_spacers": to_int(row.get("n_spacers", "0")),
            "n_spacers_unique": to_int(row.get("n_spacers_unique", "0")),
            "n_consensus_repeats": to_int(row.get("n_consensus_repeats", "0")),
            "median_repeat_len": to_float(row.get("median_repeat_len", "0")),
            "median_spacer_len": to_float(row.get("median_spacer_len", "0")),
            "n_spacers_shared": to_int(row.get("n_spacers_shared", "0")),
            "pct_spacers_shared": to_float(row.get("pct_spacers_shared", "0")),
        }

    document: Dict[str, object] = {
        "id": "mcaat_cohort_summary",
        "section_name": "MCAAT cohort summary",
        "description": (
            "Per-sample CRISPR array statistics from MCAAT, computed directly on "
            "unassembled reads. Sharing columns are exact-sequence spacer matches across "
            "the cohort."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "mcaat_cohort_summary_table",
            "title": "MCAAT: cohort summary",
            "namespace": "MCAAT",
            "col1_header": "Sample",
        },
        "headers": headers,
        "data": data,
    }
    write_yaml(path, document)


def section_arrays_bargraph(path: str, summary_rows: Sequence[Dict[str, str]]) -> None:
    rows: List[List[str]] = []
    for row in summary_rows:
        sample = row.get("sample", "")
        if not sample:
            continue
        rows.append(
            [
                sample,
                str(to_int(row.get("n_arrays", "0"))),
                str(to_int(row.get("n_spacers_unique", "0"))),
            ]
        )
    write_mqc_tsv(
        path,
        meta=[
            ("id", "'mcaat_arrays_bargraph'"),
            ("section_name", "'MCAAT arrays and spacers'"),
            (
                "description",
                "'CRISPR arrays and distinct spacer sequences detected per sample by MCAAT.'",
            ),
            ("plot_type", "'bargraph'"),
        ],
        pconfig=[
            ("id", "'mcaat_arrays_bargraph_plot'"),
            ("title", "'MCAAT: arrays and unique spacers per sample'"),
            ("ylab", "'Count'"),
            ("cpswitch", "False"),
        ],
        header_row=["Sample", "CRISPR arrays", "Unique spacers"],
        rows=rows,
    )


def _length_histograms(
    cohort_rows: Iterator[Dict[str, str]],
) -> Tuple[Dict[str, Dict[int, int]], Dict[str, Dict[int, int]]]:
    """Build per-sample spacer- and repeat-length histograms in one pass."""
    spacer_hist: Dict[str, Dict[int, int]] = {}
    repeat_hist: Dict[str, Dict[int, int]] = {}
    seen_arrays: Set[Tuple[str, str]] = set()
    for row in cohort_rows:
        sample = row.get("sample", "")
        if not sample:
            continue
        slen = to_int(row.get("spacer_len", "0"))
        if slen > 0:
            bucket = spacer_hist.setdefault(sample, {})
            bucket[slen] = bucket.get(slen, 0) + 1
        # Count each array's consensus repeat once, not once per spacer.
        array_key = (sample, row.get("array_id", ""))
        if array_key not in seen_arrays:
            seen_arrays.add(array_key)
            rlen = to_int(row.get("consensus_repeat_len", "0"))
            if rlen > 0:
                bucket = repeat_hist.setdefault(sample, {})
                bucket[rlen] = bucket.get(rlen, 0) + 1
    return spacer_hist, repeat_hist


def section_length_linegraph(
    path: str,
    histograms: Dict[str, Dict[int, int]],
    plot_id: str,
    section_name: str,
    description: str,
    xlab: str,
    ylab: str,
) -> bool:
    if not histograms:
        return False
    x_values = sorted(set(x for hist in histograms.values() for x in hist))
    if not x_values:
        return False
    rows: List[List[str]] = []
    for sample in sorted(histograms):
        hist = histograms[sample]
        rows.append([sample] + [str(hist.get(x, 0)) for x in x_values])
    write_mqc_tsv(
        path,
        meta=[
            ("id", "'{0}'".format(plot_id)),
            ("section_name", "'{0}'".format(section_name)),
            ("description", "'{0}'".format(description)),
            ("plot_type", "'linegraph'"),
        ],
        pconfig=[
            ("id", "'{0}_plot'".format(plot_id)),
            ("title", "'{0}'".format(section_name)),
            ("xlab", "'{0}'".format(xlab)),
            ("ylab", "'{0}'".format(ylab)),
            ("categories", "False"),
        ],
        header_row=["Sample"] + [str(x) for x in x_values],
        rows=rows,
    )
    return True


def section_sharing_heatmap(
    path: str, samples: Sequence[str], sharing_rows: Sequence[Dict[str, str]]
) -> bool:
    if len(samples) < 2:
        info("fewer than two samples in the cohort; the sharing heatmap is skipped")
        return False

    matrix: Dict[Tuple[str, str], int] = {}
    for row in sharing_rows:
        a = row.get("sample_a", "")
        b = row.get("sample_b", "")
        if not a or not b:
            continue
        n = to_int(row.get("n_shared_spacers", "0"))
        matrix[(a, b)] = n
        matrix[(b, a)] = n

    ordered = sorted(samples)
    rows: List[List[str]] = []
    for a in ordered:
        cells: List[str] = [a]
        for b in ordered:
            if a == b:
                cells.append("0")
            else:
                cells.append(str(matrix.get((a, b), 0)))
        rows.append(cells)

    write_mqc_tsv(
        path,
        meta=[
            ("id", "'mcaat_sharing'"),
            ("section_name", "'MCAAT spacer sharing'"),
            (
                "description",
                "'Number of exact-match spacer sequences shared between each pair of "
                "samples. The diagonal is set to zero so it does not dominate the colour "
                "scale.'",
            ),
            ("plot_type", "'heatmap'"),
        ],
        pconfig=[
            ("id", "'mcaat_sharing_plot'"),
            ("title", "'MCAAT: shared spacers between samples'"),
            ("xlab", "'Sample'"),
            ("ylab", "'Sample'"),
            ("square", "True"),
        ],
        header_row=[""] + ordered,
        rows=rows,
    )
    return True


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mcaat_multiqc_sections.py",
        description=(
            "Render MCAAT cohort tables as MultiQC custom-content sections "
            "(*_mqc.tsv / *_mqc.yaml, discovered by MultiQC's magic-filename mechanism)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--summary",
        required=True,
        metavar="FILE",
        help="cohort_summary.tsv from mcaat_aggregate.py.",
    )
    parser.add_argument(
        "--cohort-tsv",
        metavar="FILE",
        help="cohort_arrays.tsv; enables the spacer- and repeat-length distribution plots.",
    )
    parser.add_argument(
        "--sharing", metavar="FILE", help="cohort_spacer_sharing.tsv; enables the sharing heatmap."
    )
    parser.add_argument(
        "--outdir", default=".", metavar="DIR", help="Output directory (default: .)."
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    outdir = args.outdir
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    def out(name: str) -> str:
        return os.path.join(outdir, name)

    summary_rows: List[Dict[str, str]] = []
    if os.path.isfile(args.summary):
        summary_rows = list(read_tsv(args.summary))
    else:
        warn("{0} does not exist; an empty cohort summary section is written".format(args.summary))

    samples = [r.get("sample", "") for r in summary_rows if r.get("sample")]
    if not samples:
        warn("cohort summary contains no samples; only the empty summary section is written")

    # Always written, even for an empty cohort: the caller declares `*_mqc.{tsv,yaml}` as a
    # required output, so emitting nothing at all would turn a legitimate zero-result into
    # a "Missing output file" task failure.
    section_summary_table(out("mcaat_cohort_summary_mqc.yaml"), summary_rows)
    written = 1

    if not samples:
        info("wrote {0} MultiQC custom-content file(s) to {1}".format(written, outdir or "."))
        return 0

    section_arrays_bargraph(out("mcaat_arrays_bargraph_mqc.tsv"), summary_rows)
    written += 1

    if args.cohort_tsv and os.path.isfile(args.cohort_tsv):
        spacer_hist, repeat_hist = _length_histograms(read_tsv(args.cohort_tsv))
        if section_length_linegraph(
            out("mcaat_spacerlen_mqc.tsv"),
            spacer_hist,
            "mcaat_spacerlen",
            "MCAAT spacer length distribution",
            "Distribution of CRISPR spacer lengths per sample, as detected by MCAAT.",
            "Spacer length (bp)",
            "Spacers",
        ):
            written += 1
        if section_length_linegraph(
            out("mcaat_repeatlen_mqc.tsv"),
            repeat_hist,
            "mcaat_repeatlen",
            "MCAAT repeat length distribution",
            "Distribution of consensus repeat lengths per sample, counted once per array.",
            "Consensus repeat length (bp)",
            "Arrays",
        ):
            written += 1
    elif args.cohort_tsv:
        warn("{0} does not exist; length distributions skipped".format(args.cohort_tsv))

    if args.sharing and os.path.isfile(args.sharing):
        if section_sharing_heatmap(
            out("mcaat_sharing_mqc.tsv"), samples, list(read_tsv(args.sharing))
        ):
            written += 1
    elif args.sharing:
        warn("{0} does not exist; sharing heatmap skipped".format(args.sharing))

    info("wrote {0} MultiQC custom-content file(s) to {1}".format(written, outdir or "."))
    return 0


if __name__ == "__main__":
    sys.exit(main())
