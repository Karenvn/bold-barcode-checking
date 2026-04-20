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

## Citation

If you use this workflow in a methods note, report, or manuscript, cite the underlying resources:

- NCBI Datasets for mitochondrial accession discovery and organelle sequence retrieval: O'Leary NA, Cox E, Holmes JB, Anderson WR, Falk R, Hem V, Tsuchiya MTN, Schuler GD, Zhang X, Torcivia J, Ketter A, Breen L, Cothran J, Bajwa H, Tinne J, Meric PA, Hlavina W, Schneider VA. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data. 2024;11(1):732. doi:10.1038/s41597-024-03571-y.
- BLAST for COI region identification: Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990;215(3):403-410. doi:10.1016/S0022-2836(05)80360-2.
- BOLD Systems for barcode matching: Ratnasingham S, Hebert PDN. BOLD: The Barcode of Life Data System. Mol Ecol Notes. 2007;7(3):355-364. doi:10.1111/j.1471-8286.2006.01678.x.
- If you discuss BIN assignments explicitly, also cite: Ratnasingham S, Hebert PDN. A DNA-Based Registry for All Animal Species: The Barcode Index Number (BIN) System. PLoS One. 2013;8(7):e66213. doi:10.1371/journal.pone.0066213.

## Notes

- NCBI and BOLD services can change response formats or rate limits over time.
- The workflow queries public BOLD data only.
- BOLD relies on correct species identification of comparator records, so this workflow is a useful check rather than a definitive identification.
- BIN metrics are a best-effort scrape from the public BOLD cluster page and may be absent for some hits.
