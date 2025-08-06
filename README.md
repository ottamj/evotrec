![Python Version](https://img.shields.io/pypi/pyversions/hammingdist)
![License](https://img.shields.io/github/license/ottamj/evotrec)

# EVOtRec

Topological recurrence analysis of genome alignments using persistent homology.

## Overview

EVOtRec is a Python module designed for topological recurrence analysis of genome alignments using persistent homology. It processes sequence alignment data and returns the topological recurrence index (tRI) of single nucleotide variations (SNVs) recorded in the alignment.

TODO: Describe the purpose of EVOtRec and how it can be used in genomic studies.


## Prerequisites

- **Python**: Version 3.8 or higher
- **C++ Compiler**: Required for building Ripser (gcc/clang)
- **Make**: For building Ripser from source

## Installation

### 1. Clone the Repository

```sh
git clone https://github.com/ottamj/evotrec.git
cd evotrec
```

### 2. Install Python Dependencies

The Python dependencies are:
- **Biopython**: For sequence analysis and FASTA file handling
- **hammingdist**: For efficient Hamming distance calculations

To install these dependencies run:
```sh
pip install -r requirements.txt
```


### 3. Install Ripser CLI (with Representative Cycles)

EVOtRec requires a special version of Ripser with representative cycles:

```sh
git clone --branch tight-representative-cycles https://github.com/Ripser/ripser.git
cd ripser
make
./ripser examples/sphere_3_192.lower_distance_matrix
```

Add the `ripser` executable to your PATH by adding this line to your `.bashrc` or `.bash_profile`:

```sh
export PATH=$PATH:/path/to/ripser
```

For detailed installation instructions, see the [Ripser GitHub repository](https://github.com/Ripser/ripser/tree/tight-representative-cycles?tab=readme-ov-file#building).

### 4. Verify Installation on Example Data

You can test EVOtRec using the provided example alignment in the `examples/` directory.

```sh
python evotrec.py examples/example.fasta "hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30|China"
```

The command produces the following console output:

```
===============================================================
EVOtRec -- Topological recurrence analysis of genome alignments
(c) 2024 Andreas Ott
===============================================================

Preparing examples/example.fasta...

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

Results written to examples/example.csv.
```

After successful execution, the following text files will be created:

```
examples/
├── example.csv      # tRI analysis results
├── example.dist     # Hamming distance matrix
├── example.ripser   # Ripser persistent homology output
└── example.timedist # Time-transformed distance matrix
```


## Usage

### Input Requirements

1. **FASTA File**: A multiple sequence alignment file in FASTA format
2. **Reference Sequence ID**: The identifier of the reference sequence within the FASTA file
3. **Time Series Data** (optional): For temporal analysis, sequence headers must contain date information in the format `|YYYY-MM-DD|`

### Basic Usage

```sh
python evotrec.py <input_fasta> <refseq_id> [--timeseries]
```

- `<input_fasta>`: Path to the input FASTA file containing sequence alignments
- `<refseq_id>`: Reference sequence identifier as it appears in the FASTA file
- `--timeseries`: Optional flag to enable time series analysis (requires date information in headers)


### Output

The script generates several text files with the same base name as the input file and the following filename extensions:

- **`<input_fasta>.dist`**: Hamming distance matrix between all sequences
- **`<input_fasta>.timedist`**: Rips-transformed distance matrix (generated only with `--timeseries` flag)
- **`<input_fasta>.ripser`**: Raw output from Ripser persistent homology computation
- **`<input_fasta>.csv`**: Final results containing topological recurrence index (tRI) analysis

The structure of `<input_fasta>.csv` depends on whether the `--timeseries` flag is used:

- `POS`: Genomic position of the variation
- `REF`: Reference nucleotide at this position
- `ALT`: Alternative nucleotide (the mutation)

**Without `--timeseries` flag:**
- `TRI`: Topological recurrence index value for this variation

**With `--timeseries` flag:**
- `1,2,3,...,N`: Time-binned columns representing each day in the time range, where values indicate the topological recurrence index at each time point (0 = not present, higher values = increasing topological significance)


## Support

For questions, issues, or feature requests:

1. Check the [Issues](https://github.com/ottamj/evotrec-prelim/issues) page
2. Search existing issues before creating a new one
3. Provide detailed information including:
   - Operating system and Python version
   - Input data characteristics
   - Error messages (if any)
   - Steps to reproduce the issue

## Contributing

We welcome contributions! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes and add tests
4. Run the test suite (`pytest tests/`)
5. Commit your changes (`git commit -m 'Add some amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

### Testing

Run the pytest suite to verify all components produce the desired outputs:

```sh
pytest tests/
```

Or run tests with verbose output:

```sh
pytest -v tests/
```

The test suite includes:
- **Unit tests**: Testing individual functions with mocked data
- **Integration tests**: End-to-end testing with real example data
- **Output validation**: Comparing results against expected outputs

### Code Style

- Follow PEP 8 for Python code style
- Use meaningful variable and function names
- Add docstrings for new functions
- Include unit tests for new functionality

## License

This project is licensed under the terms specified in the LICENSE file.

## Citation

If you use EVOtRec in your research, please cite the following references:

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
