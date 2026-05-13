# seqcal

seqcal is a primer-guided sequence extraction and copy-number calculation tool.

It uses `seqkit locate` to find forward and reverse primer binding sites, pairs
compatible hits, reports copy numbers, and writes extracted target sequences.

## Installation

Create the conda environment with the Python dependencies and `seqkit`:

```bash
conda env create -f environment.yml
conda activate seqcal
pip install -e .
```

## Usage

```bash
seqcal \
  -f genomes/ \
  -o results/ \
  --forward-primer GATGAAGARCGYAGYRAA \
  --reverse-primer TCCTCCGCTTATTGATATGC \
  -t 8
```

You can pass either a directory containing `.fna`/`.fasta` files or a single
`.fna`/`.fasta` file to `-f`.

## Outputs

- `copies.seqcal`: copy number table
- `seqhit.fasta`: extracted target sequences
- `target_regions.tsv`: target region coordinates

## Options

Extract sequences including the primer binding sites:

```bash
seqcal \
  -f genomes/ \
  -o results/ \
  --forward-primer GATGAAGARCGYAGYRAA \
  --reverse-primer TCCTCCGCTTATTGATATGC \
  --include-primers
```

Filter extracted sequences by length:

```bash
seqcal \
  -f genomes/ \
  -o results/ \
  --forward-primer GATGAAGARCGYAGYRAA \
  --reverse-primer TCCTCCGCTTATTGATATGC \
  --min-length 100 \
  --max-length 500
```

## Development

Run tests with:

```bash
pytest
```
