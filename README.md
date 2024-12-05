
# EVOtRec
![Python Version](https://img.shields.io/pypi/pyversions/hammingdist)
![License](https://img.shields.io/github/license/user/repo)

## Overview
EVOTREC is a Python module designed for topological recurrence analysis of genome alignments using persistent homology. It processes sequence alignment data and returns the topological recurrence index (tRI) of the single nucleotide variations (SNV) recorded in the alignment.

## Dependencies
The module requires the following dependencies:
- Python 3.x
- Biopython
- hammingdist
- Ripser CLI (with representative cycles)

You can install the *python dependencies* using `pip install -r requirements.txt`.

To install Ripser with tight representative cycles:
```sh
git clone --branch tight-representative-cycles https://github.com/Ripser/ripser.git
cd ripser
make
./ripser examples/sphere_3_192.lower_distance_matrix
```
and add the `ripser` executable to your PATH, e.g. by adding the following line to your `.bashrc` or `.bash_profile`:
```sh
export PATH=$PATH:/path/to/ripser
```
For more information see the instructions provided in the [Ripser GitHub repository](https://github.com/Ripser/ripser/tree/tight-representative-cycles?tab=readme-ov-file#building).

## Usage
To use the EVOTREC module, follow these steps:

1. **Install Dependencies**:

2. **Prepare Input Files**:
    - **FASTA File**: Prepare a FASTA file containing the sequence alignments.
    - **Reference Sequence ID**: Identify the reference sequence ID within the FASTA file.
    - **Time Series Analysis**: To perform time series analysis, the sequence headers must contain the date information in the format `|YYYY-MM-DD|`.

3. **Run the EVOTREC Script**:
    ```sh
    python evotrec.py <input_fasta> <refseq_id> [--timeseries]
    ```
    - `<input_fasta>`: Path to the input FASTA file.
    - `<refseq_id>`: Reference sequence ID in the FASTA file.
    - `--timeseries`: Optional flag to enable timeseries analysis.

4. **Output Files**:
    The script generates several output files:
    - `<input_fasta>.dist`: Hamming distance file.
    - `<input_fasta>.timedist`: Rips transformed distance file (if `--timeseries` is used).
    - `<input_fasta>.ripser`: Ripser output file.
    - `<input_fasta>.csv`: tRI analysis results in CSV format.


## Example
You can run evotrec on the provided test data as follows
```sh
python evotrec.py tests/example.fasta "hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30|China" --timeseries
```

This produces the following console output
```sh
===============================================================
EVOtRec -- Topological recurrence analysis of genome alignments
(c) 2024 Andreas Ott
===============================================================

Preparing tests/example.fasta...

Start date: 2019-12-30
End date: 2021-09-30
Time range: 641 days

Number of sequences: 3594

Hammingdist...
# hammingdist status
MuRiT...
Ripser...
Analyzing cycles...
Computing tRI...

Results written to tests/example.csv.
```
and output files
```sh
└── tests
    ├── example.csv
    ├── example.dist
    ├── example.ripser
    └── example.timedist
```

### Citation
TODO: (Max) Update arxiv link, bibtex citation

If you use EVOTREC in your research, please cite the following references:

```bibtex
@article{bleher2023topological,
  title={Topological data analysis identifies emerging adaptive mutations in SARS-CoV-2},
  author={Bleher, Michael and Hahn, Lukas and Neumann, Maximilian and Patino-Galindo, Juan Angel and Carriere, Mathieu and Bauer, Ulrich and Rabadan, Raul and Ott, Andreas},
  journal={arXiv preprint arXiv:2106.07292},
  year={2023},
  url={https://arxiv.org/abs/2106.07292}
}

@article{neumann2022murit,
  title={MuRiT: Efficient Computation of Pathwise Persistence Barcodes in Multi-Filtered Flag Complexes via Vietoris-Rips Transformations},
  author={Neumann, Maximilian and Bleher, Michael and Hahn, Lukas and Braun, Samuel and Obermaier, Holger and Soysal, Mehmet and Caspart, Ren{\'e} and Ott, Andreas},
  journal={arXiv preprint arXiv:2207.03394},
  year={2022},
  url={https://arxiv.org/abs/2207.03394}
}
```
