# .fa-to-.json
ProteinMPNN FASTA → AlphaFold Server JSON (alphafoldserver v1)
Batch-convert single-chain FASTA (.fa/.fasta) sequences into AlphaFold Server input JSON jobs (dialect: alphafoldserver, v1)
````markdown
# .fa to .json

Convert **ProteinMPNN** FASTA outputs (`.fa/.fasta`) into **AlphaFold Server** input JSON format (dialect: `alphafoldserver`, version `1`).

-  **Single-chain** FASTA supported (each FASTA record → one AlphaFold Server “job”)
-  Batch processing: **many sequences in one FASTA → one JSON file** containing a list of jobs
-  No third-party dependencies (pure Python 3)

> Not affiliated with AlphaFold / DeepMind / Google. This is a small utility script for preparing inputs.

---

## Why this repo?

ProteinMPNN typically outputs many designed sequences in a FASTA file.  
AlphaFold Server expects a JSON file where each “job” contains a protein chain sequence.

This script bridges that gap:
- Reads `>header` + `SEQUENCE` records from FASTA
- Validates sequences (20 standard amino acids by default)
- Writes a single JSON file containing a list of AlphaFold Server jobs

---

## Requirements

- Python **3.8+** (Linux recommended; macOS/Windows also works)
- No extra packages required

Check:
```bash
python3 --version
````



## One-line, copy-paste command (most common)

> This is the most common usage for ProteinMPNN outputs:
> **skip the first record**, generate **clean sequential job names**, and write everything into one JSON.

```bash
python mpnn_fa_to_afserver_json.py --input design_cut_2.fa --output afserver_jobs.json --job-prefix design_cut_2 --name-mode index --skip-first
```

If your FASTA is in a directory:

```bash
python mpnn_fa_to_afserver_json.py --input dir/design_cut_2.fa --output dir/afserver_jobs.json --job-prefix design_cut_2 --name-mode index --skip-first
```

---

## Quick start

### Minimal (only required args)

```bash
python mpnn_fa_to_afserver_json.py \
  --input  path/to/design.fa \
  --output path/to/afserver_jobs.json
```

### Recommended for ProteinMPNN outputs (skip the first record)

ProteinMPNN outputs often place the “original/native” sequence as the first record.

```bash
python mpnn_fa_to_afserver_json.py \
  --input  dir/design_cut_2.fa \
  --output dir/afserver_jobs.json \
  --job-prefix design_cut_2 \
  --skip-first
```

> **Shell tip:** the backslash `\` must be the **very last character** on the line (no trailing spaces).

---

## Input FASTA format

Any standard FASTA is fine:

```text
>design_0001
ACDEFGHIKLMNPQRSTVWY
>design_0002
MKT...
```

* Sequences can be wrapped across multiple lines (supported)
* By default, only the 20 standard amino acids are allowed:
  `ACDEFGHIKLMNPQRSTVWY`

---

## Output JSON format

The output JSON is a **list** of job objects. Each FASTA record becomes one job:

* `dialect`: `"alphafoldserver"`
* `version`: `1`
* `sequences`: `[{"proteinChain": {"sequence": "...", "count": 1}}]`

Example (shortened):

```json
[
  {
    "name": "design_cut_2_0001",
    "modelSeeds": [],
    "sequences": [{"proteinChain": {"sequence": "ACDE...", "count": 1}}],
    "dialect": "alphafoldserver",
    "version": 1
  }
]
```

---

## Naming: `--job-prefix` and `--name-mode`

AlphaFold Server displays each submitted item by its `name`. This script constructs names as:

* `--job-prefix` (default: input filename stem)
* plus either:

  * `--name-mode header` (default): derive suffix from FASTA header (more informative)
  * `--name-mode index`: purely sequential (`PREFIX_0001`, `PREFIX_0002`, ...) — **cleanest**

Recommended if you want clean, predictable names:

```bash
python mpnn_fa_to_afserver_json.py \
  --input design.fa \
  --output afserver_jobs.json \
  --job-prefix my_designs \
  --name-mode index \
  --skip-first
```

---

## Full CLI usage

Print help:

```bash
python mpnn_fa_to_afserver_json.py -h
```

### All options (overview)

| Option                                | What it does                                     | Typical use                         |
| ------------------------------------- | ------------------------------------------------ | ----------------------------------- |
| `--input`, `-i`                       | Input FASTA path                                 | required                            |
| `--output`, `-o`                      | Output JSON path                                 | required                            |
| `--job-prefix`                        | Prefix for job `name`                            | optional (recommended)              |
| `--name-mode {header,index}`          | Name strategy                                    | optional                            |
| `--skip-first`                        | Skip the first FASTA record                      | often useful for ProteinMPNN        |
| `--max-jobs N`                        | Only export first N jobs (after skip)            | batching                            |
| `--count N`                           | `proteinChain.count` (default 1)                 | usually keep default                |
| `--max-template-date YYYY-MM-DD`      | Set `maxTemplateDate` (optional AF Server field) | advanced                            |
| `--use-structure-template true/false` | Set `useStructureTemplate` (optional field)      | advanced                            |
| `--seed S` (repeatable)               | Add model seeds to `modelSeeds`                  | usually leave empty                 |
| `--allow-nonstandard`                 | Allow non-20AA letters (X/B/Z/...)               | not recommended (server may reject) |

---

## Practical examples

### Export only 100 designs (batching)

```bash
python mpnn_fa_to_afserver_json.py \
  --input design.fa \
  --output afserver_jobs_100.json \
  --skip-first \
  --max-jobs 100
```

### Use absolute paths

```bash
python mpnn_fa_to_afserver_json.py \
  --input  /data/proteinmpnn/design_cut_2.fa \
  --output /data/afserver_inputs/afserver_jobs.json \
  --job-prefix design_cut_2 \
  --skip-first
```

### Allow nonstandard AA letters (try only if you must)

```bash
python mpnn_fa_to_afserver_json.py \
  --input design.fa \
  --output afserver_jobs.json \
  --allow-nonstandard
```

---

## Notes & best practices

* **Server limits:** AlphaFold Server may limit number/size of jobs per submission.
  If your FASTA has hundreds/thousands of sequences, export in batches using `--max-jobs`, or split the FASTA first.
* **Sequence validation:** Default is strict 20AA validation to prevent server rejections.
* **Single-chain only:** This version assumes each FASTA record is one single chain. (Multi-chain support can be added later.)

---

## Troubleshooting

### “Found non-standard AA letters…”

Your sequence contains letters outside the standard 20 amino acids (e.g., `X`).

* Fix upstream (recommended), or
* rerun with `--allow-nonstandard` (may still be rejected by the server)

### Output JSON is empty

You probably used `--skip-first` on a FASTA with only one record, or filtered everything with `--max-jobs`.

* Remove `--skip-first` or increase `--max-jobs`


## Contributing

Issues and PRs are welcome:

* multi-chain FASTA support
* automatic chunking into multiple JSON files
* richer header parsing for ProteinMPNN naming

This work is done at [Yang Lab](https://jieyang-lab.com/) at [UVA](https://www.virginia.edu/) school of medicine, [Biochemistry-Molecular Genetics](https://med.virginia.edu/bmg/)
