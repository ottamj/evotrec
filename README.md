# EVOTREC

## Overview
EVOTREC is a Python module designed for topological recurrence analysis of genome alignments using persistent homology. It processes sequence alignment data, calculates Hamming distances, converts them to time distances, and performs persistent homology analysis to identify single nucleotide variant (SNV) cycles. The module then retrieves sequences and mutations in these cycles and performs topological recurrence index (tRI) analysis. 

## Usage
To use the EVOTREC module, follow these steps:

1. **Install Dependencies**:
    Ensure you have all the required dependencies installed. You can install them using the following command:
    ```sh
    pip install -r requirements.txt
    ```

2. **Prepare Input Files**:
    - **FASTA File**: Prepare a FASTA file containing the sequence alignments.
    - **Reference Sequence ID**: Identify the reference sequence ID within the FASTA file. Mutations are determined in the form `(Pos, Ref, Alt)` with respect to this reference sequence.

3. **Run the EVOTREC Script**:
    Execute the `evotrec.py` script with the required arguments:
    ```sh
    python evotrec.py <input_fasta> <refseq_id> [--timeseries]
    ```
    - `<input_fasta>`: Path to the input FASTA file.
    - `<refseq_id>`: Reference sequence ID in the FASTA file.
    - `--timeseries`: Optional flag to enable timeseries analysis. Requires sequence headers to contain dates in the format `YYYY-MM-DD`.

4. **Output Files**:
    The script generates several output files:
    - `<input_fasta>.dist`: Hamming distance file.
    - `<input_fasta>.timedist`: Time distance file (if `--timeseries` is used).
    - `<input_fasta>.ripser`: Ripser output file.
    - `<input_fasta>.csv`: tRI analysis results in CSV format.

## Dependencies
The module requires the following dependencies:
- Python 3.x
- Biopython
- hammingdist
- Ripser

You can install these dependencies using the provided [requirements.txt](http://_vscodecontentref_/0) file:
```sh
pip install -r requirements.txt
```

TODO: Add installation instructions and links for ripser

## Testing
TODO: This is internal info, remove before publishing.

You can test the functionality of EVOTREC using the provided unit tests.

Either run pytest
```sh
pytest
```

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
# hammingdist (instance dependent lines)
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
    ├── example.fasta
    ├── example.ripser
    └── example.timedist
```

### Citation
If you use EVOTREC in your research, please cite the following references:

```bibtex
@article{bleher2023topological,
  title={Topological data analysis identifies emerging adaptive mutations in SARS-CoV-2},
  author={Michael Bleher and Lukas Hahn and Maximilian Neumann and Juan Angel Patino-Galindo and Mathieu Carriere and Ulrich Bauer and Raul Rabadan and Andreas Ott},
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
