![Python Version](https://img.shields.io/badge/python-3.10-blue.svg)
![License](https://img.shields.io/github/license/ottamj/evotrec)


# EVOtRec
Topological recurrence analysis of genome alignments using topological data analysis.<br>
© 2020-26 Andreas Ott, Michael Bleher, Maximilian Neumann

## Description

EVOtRec (topological Recurrence in EVOlution) is a Python module for the efficient and scalable topological recurrence analysis of large genome alignments. It uses persistent homology to compute the topological recurrence index (tRI) of single nucleotide variations in a nucleotide sequence alignment. EVOtRec enables the inference of the dynamics of convergently evolving genomic variants directly from topological patterns in the genomic dataset, without requiring the construction of a phylogenetic tree. For a paper describing EVOtRec and its applications to molecular evolution in more detail, see [Bleher et al, 2026](https://arxiv.org/abs/2106.07292).

## Prerequisites

- **Python**: version 3.8 or higher
- **C++ compiler**: for building Ripser (e.g. gcc>=11.4.0)
- **Make**: for building Ripser (>=4.3)

## Installation

The installation of all required software should take no more than 10 minutes, and the expected runtime for each demo example is less than 1 minute (tested on a system with an 8-core CPU and 16 GB RAM).

### 1. Clone the repository

```sh
git clone https://github.com/ottamj/evotrec.git
cd evotrec
```

### 2. Install Python dependencies

The Python dependencies are:
- **Biopython**: for sequence analysis and FASTA file handling (==1.86)
- **hammingdist**: for efficient Hamming distance calculations (==1.4.0)

To install these dependencies run:
```sh
pip install -r requirements.txt
```

Note: While we do not expect breaking changes in future versions of these dependencies, we cannot guarantee compatibility. It is recommended to test your installation with the provided example data.

### 3. Install Ripser (with representative cycles)

EVOtRec requires a special version of Ripser that also computes representative cycles:

```sh
git clone --branch tight-representative-cycles https://github.com/Ripser/ripser.git
cd ripser
make
```

This will create the `ripser-representatives` executable in the `ripser/` directory.
Add it to your PATH by adding the following to your `.bashrc` or `.bash_profile`:

```sh
export PATH=$PATH:/path/to/ripser
```

You can check that Ripser was installed successfully by running it on the example distance matrix provided in the ripser repository:

```sh
ripser-representatives examples/sphere_3_192.lower_distance_matrix
```

For more detailed installation instructions, see the [Ripser GitHub repository](https://github.com/Ripser/ripser/tree/tight-representative-cycles?tab=readme-ov-file#building).

### 4. Verify installation on example data

You can test EVOtRec using the example alignments in the `examples/` directory.

```sh
python evotrec.py examples/sars-cov-2_spike_gene.fasta "NC_045512.2|Severeacuterespiratorysyndromecoronavirus2isolateWuhan-Hu-1,completegenome|China|2019-12-30" --timeseries
```

The command produces the following output:

```
===============================================================
EVOtRec -- Topological recurrence analysis of genome alignments
(c) 2020-26 Andreas Ott, Michael Bleher, Maximilian Neumann
===============================================================

Preparing examples/sars-cov-2_spike_gene.fasta...

Start date: 2019-12-30
End date: 2021-01-03
Time range: 371 days

Number of sequences: 3274

Hammingdist...
# hammingdist <status>
MuRiT...
Ripser...
Analyzing cycles...
Computing tRI...

Results written to examples/sars-cov-2_spike_gene.csv.
```

After successful execution, the following text files will be created:

```
examples/
├── sars-cov-2_spike_gene.csv
├── sars-cov-2_spike_gene.dist
├── sars-cov-2_spike_gene.ripser
└── sars-cov-2_spike_gene.timedist
```


## Usage

### Input

1. **FASTA file**: a nucleotide sequence alignment in FASTA format
2. **Reference sequence header**: the header (excluding the leading ">") of the reference sequence in the FASTA file, enclosed in double quotes
3. **`--timeseries`** (optional): When using the `--timeseries` flag, EVOtRec requires that all sequence headers in the FASTA file contain a date in the following format: `>seq-id|field1|field2|...|YYYY-MM-DD`

**Requirements for sequence dates:**
- **Format**: must be exactly `YYYY-MM-DD` (ISO 8601 format)
- **Position**: must be the last field when splitting by `|`
- **Chronological order**: all sequence dates must be on or after the reference sequence date

If any sequence header does not follow this format or contains an invalid date, EVOtRec will raise a `ValueError` with details about the problematic sequence.

**Remarks:**
- For best computational performance, we recommend that the input FASTA file is deduplicated by sequence and sorted in reverse chronological order (newest sequences first).
- Hammingdist requires that input nucleotide sequences contain only the characters A, C, G, T or -.

### Basic usage

```sh
python evotrec.py <input_fasta> "<refseq_header>" [--timeseries]
```

- `<input_fasta>`: path to the input FASTA file containing sequence alignments
- `<refseq_header>`: reference sequence header as it appears in the FASTA file (excluding the leading ">")
- `--timeseries`: optional flag to enable time series analysis (requires a date for each sequence)


### Output

The script generates several text files with the same base name as the input file and the following filename extensions:

- **`<input_fasta>.dist`**: Hamming distance matrix between all sequences in sparse format
- **`<input_fasta>.timedist`**: Rips-transformed distance matrix (generated only when `--timeseries` flag is used)
- **`<input_fasta>.ripser`**: raw output from Ripser persistent homology computation
- **`<input_fasta>.csv`**: final results of topological recurrence analysis (tRI list)

The structure of `<input_fasta>.csv` depends on whether the `--timeseries` flag is used:

- `POS`: genomic position of the variation
- `REF`: reference nucleotide at this position
- `ALT`: alternative nucleotide (the mutation)

**Without `--timeseries` flag:**
- `tRI`: topological recurrence index of the variation

**With `--timeseries` flag:**
- `1,2,3,...,N`: time-binned columns representing days in the time range, and columns containing tRI values at each day


## Support

For questions, issues, or feature requests:

1. Check the [Issues](https://github.com/ottamj/evotrec-prelim/issues) page
2. Search existing issues before creating a new one
3. Provide detailed information including:
   - operating system and Python version
   - input data characteristics
   - error messages (if any)
   - steps to reproduce the issue

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
- **Unit tests**: testing individual functions with mocked data
- **Integration tests**: end-to-end testing with real example data
- **Output validation**: comparing results against expected outputs

### Code style

- Follow PEP 8 for Python code style
- Use meaningful variable and function names
- Add docstrings for new functions
- Include unit tests for new functionality

## License

This project is licensed under the terms specified in the LICENSE file.

## Citation

If you use EVOtRec in your research, please cite the following references:

```bibtex
@article{bleher2026topological,
  title={Ultrafast topological data analysis reveals pandemic-scale dynamics of convergent evolution},
  author={Bleher, Michael and Hahn, Lukas and Neumann, Maximilian and Ardern, Zachary and Patino-Galindo, Juan Angel and Carriere, Mathieu and Bauer, Ulrich and Rabadan, Raul and Ott, Andreas},
  journal={arXiv preprint arXiv:2106.07292},
  year={2026},
  url={https://arxiv.org/abs/2106.07292}
}

@article{neumann2022murit,
  title={MuRiT: Efficient Computation of Pathwise Persistence Barcodes in Multi-Filtered Flag Complexes via Vietoris-Rips Transformations},
  author={Neumann, Maximilian and Bleher, Michael and Hahn, Lukas and Braun, Samuel and Obermaier, Holger and Soysal, Mehmet and Caspart, Rene and Ott, Andreas},
  journal={arXiv preprint arXiv:2207.03394},
  year={2022},
  url={https://arxiv.org/abs/2207.03394}
}
```
