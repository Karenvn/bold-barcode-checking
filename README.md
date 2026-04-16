# BOLD Barcode Checking Workflow

This repository packages a command-line workflow for:

1. finding the mitochondrial sequence linked to a `GCA_...` assembly accession
2. extracting the COI / COX1 region with NCBI BLAST
3. querying BOLD for public barcode matches and BIN assignments
4. retrieving BIN-level distance and nearest-neighbour metadata when available


## Requirements

- Python 3.10 or newer
- an internet connection to reach NCBI and BOLD
- `ENTREZ_EMAIL` set to a valid email address before running the workflow
- optionally, `ENTREZ_API_KEY` for higher NCBI request limits

## Install

```bash
python3 -m pip install -e .
```

Or install dependencies without packaging:

```bash
python3 -m pip install biopython pandas requests
```

## Usage

Run one accession directly:

```bash
bold-coi GCA_964291005.1
```

Run several accessions and save the full table:

```bash
bold-coi GCA_964291005.1 GCA_000001405.29 --wait-seconds 8 --output-csv results.csv
```

Run from a text file:

```bash
bold-coi --accessions-file accessions.txt
```

Input files should contain one accession per line, for example:

```text
GCA_964291005.1
GCA_000001405.29
```

## Output

The workflow reports:

- mitochondrial accession and sequence length
- inferred species name from the mitochondrial record description
- extracted COI coordinates, strand, and BLAST identity
- top BOLD match, similarity, process ID, and BIN URI
- a likely self-hit flag when the top BOLD match appears to be the query itself
- optional BIN distance / nearest-neighbour metrics when BOLD exposes them

Use `--output-csv` to retain the full result table.

## Notes

- NCBI and BOLD services can change response formats or rate limits over time.
- The workflow queries public BOLD data only.
- BIN metrics are a best-effort scrape from the public BOLD cluster page and may be absent for some hits.
