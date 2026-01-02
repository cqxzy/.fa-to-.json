#!/usr/bin/env python3
"""
Convert ProteinMPNN FASTA (.fa/.fasta) outputs (single-chain) to AlphaFold Server JSON format.

AlphaFold Server expects a JSON file that is a LIST of job dicts.
Each job contains:
  - name: str
  - modelSeeds: list[str] (can be empty)
  - sequences: list[entity dicts]
  - dialect: "alphafoldserver"
  - version: 1

We generate: one FASTA record -> one AF Server job (single proteinChain).
"""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional


AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    """Parse FASTA into list of (header, sequence). Supports wrapped sequences."""
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq_lines: List[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.replace(" ", "").replace("\t", ""))

    if header is not None:
        records.append((header, "".join(seq_lines)))

    if not records:
        raise ValueError(f"No FASTA records found in: {path}")
    return records


def normalize_and_validate_sequence(seq: str, strict_aa20: bool = True) -> str:
    seq = seq.strip().upper().replace(" ", "").replace("\t", "")
    if not seq:
        raise ValueError("Empty sequence encountered.")

    # remove common non-seq characters (just in case)
    seq = seq.replace("\r", "").replace("\n", "")

    if strict_aa20:
        bad = sorted(set([c for c in seq if c not in AA20]))
        if bad:
            raise ValueError(
                "Found non-standard AA letters not allowed by AF Server (expects 20 AA). "
                f"Bad letters: {bad}. "
                "If you are sure you want to allow them, rerun with --allow-nonstandard."
            )
    return seq


def sanitize_job_name(name: str, max_len: int = 120) -> str:
    # Replace whitespace and punctuation with underscores, keep a safe subset
    name = name.strip()
    name = re.sub(r"\s+", "_", name)
    name = name.replace(",", "_").replace(";", "_").replace(":", "_").replace("/", "_")
    name = re.sub(r"[^A-Za-z0-9_.\-]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    if not name:
        name = "job"
    return name[:max_len]


def derive_suffix_from_mpnn_header(header: str) -> str:
    """
    Try to extract a nice suffix from common ProteinMPNN headers, e.g.
      "T=0.2, sample=12, score=..., ..."
    or first token before comma like "design_cut_2, score=..., ..."
    """
    t = re.search(r"T\s*=\s*([0-9.]+)", header)
    s = re.search(r"sample\s*=\s*([0-9]+)", header)
    if t and s:
        return f"T{t.group(1)}_sample{s.group(1)}"
    if s:
        return f"sample{s.group(1)}"
    # fallback: take text before first comma
    return header.split(",")[0].strip()


def build_jobs(
    records: List[Tuple[str, str]],
    job_prefix: str,
    skip_first: bool,
    max_jobs: Optional[int],
    count: int,
    max_template_date: Optional[str],
    use_structure_template: Optional[bool],
    seeds: List[str],
    strict_aa20: bool,
    name_mode: str,
) -> List[dict]:
    jobs: List[dict] = []
    start_idx = 1 if skip_first else 0

    selected = records[start_idx:]
    if max_jobs is not None:
        selected = selected[:max_jobs]

    for i, (hdr, seq) in enumerate(selected, start=1):
        seq_norm = normalize_and_validate_sequence(seq, strict_aa20=strict_aa20)

        if name_mode == "index":
            job_name = f"{job_prefix}_{i:04d}"
        else:
            suffix = derive_suffix_from_mpnn_header(hdr)
            job_name = f"{job_prefix}_{suffix}"

        job_name = sanitize_job_name(job_name)

        protein_chain = {
            "sequence": seq_norm,
            "count": int(count),
        }
        if max_template_date:
            protein_chain["maxTemplateDate"] = max_template_date
        if use_structure_template is not None:
            protein_chain["useStructureTemplate"] = bool(use_structure_template)

        job = {
            "name": job_name,
            "modelSeeds": seeds,  # keep empty list by default (recommended)
            "sequences": [{"proteinChain": protein_chain}],
            "dialect": "alphafoldserver",
            "version": 1,
        }
        jobs.append(job)

    if not jobs:
        raise ValueError("After filtering (skip/max_jobs), no jobs remain to write.")
    return jobs


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert ProteinMPNN FASTA to AlphaFold Server JSON (alphafoldserver v1).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, type=Path, help="Input FASTA (.fa/.fasta) from ProteinMPNN.")
    p.add_argument("-o", "--output", required=True, type=Path, help="Output JSON file for AlphaFold Server.")

    p.add_argument("--job-prefix", default=None, help="Prefix for each job name. Default: input filename stem.")
    p.add_argument(
        "--name-mode",
        choices=["header", "index"],
        default="header",
        help="How to name jobs: derive from MPNN header (header) or just sequential index (index).",
    )

    p.add_argument("--skip-first", action="store_true",
                   help="Skip the first FASTA record (often the original/native sequence in ProteinMPNN output).")
    p.add_argument("--max-jobs", type=int, default=None, help="Only output the first N jobs (after skip).")

    p.add_argument("--count", type=int, default=1, help="proteinChain.count in AF Server JSON.")
    p.add_argument("--max-template-date", default=None,
                   help="Optional: proteinChain.maxTemplateDate (YYYY-MM-DD).")
    p.add_argument("--use-structure-template", default=None, choices=["true", "false"],
                   help="Optional: proteinChain.useStructureTemplate (true/false).")

    p.add_argument("--seed", action="append", default=[],
                   help="Add a model seed (uint32 as string). You can pass multiple --seed. Default empty list.")
    p.add_argument("--allow-nonstandard", action="store_true",
                   help="Allow non-20AA letters (NOT recommended; AF Server may reject).")

    return p.parse_args()


def main() -> None:
    args = parse_args()

    inp: Path = args.input
    out: Path = args.output
    if not inp.exists():
        raise FileNotFoundError(f"Input FASTA not found: {inp}")

    job_prefix = args.job_prefix if args.job_prefix else inp.stem

    use_structure_template = None
    if args.use_structure_template is not None:
        use_structure_template = (args.use_structure_template.lower() == "true")

    records = parse_fasta(inp)

    jobs = build_jobs(
        records=records,
        job_prefix=job_prefix,
        skip_first=bool(args.skip_first),
        max_jobs=args.max_jobs,
        count=args.count,
        max_template_date=args.max_template_date,
        use_structure_template=use_structure_template,
        seeds=[str(s) for s in args.seed],
        strict_aa20=(not args.allow_nonstandard),
        name_mode=args.name_mode,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        json.dump(jobs, f, indent=2)

    print(f"[OK] Wrote {len(jobs)} AlphaFold Server jobs to: {out}")


if __name__ == "__main__":
    main()
