# EVOtRec examples

This directory contains example data for demonstrating EVOtRec functionality.

## Data files

Nucleotide sequence alignments, consisting of publicly available genomes downloaded from NCBI:

- **`h5n1_pb2_gene.fasta`**: influenza H5N1 PB2 gene alignment
- **`sars-cov-2_spike_gene.fasta`**: SARS-CoV-2 spike gene alignment

## Running the examples

From the project root directory, run:

```sh
python evotrec.py examples/h5n1_pb2_gene.fasta "OQ958041.1|InfluenzaAvirus(A/AmericanWigeon/SouthCarolina/22-000345-001/2021(H5N1))segment1polymerasePB2(PB2)gene,completecds|2021-12-30" --timeseries
```

```sh
python evotrec.py examples/sars-cov-2_spike_gene.fasta "NC_045512.2|Severeacuterespiratorysyndromecoronavirus2isolateWuhan-Hu-1,completegenome|China|2019-12-30" --timeseries
```

This will generate the following output files in the `examples/` directory:
- **`<filename>.csv`**: tRI analysis results
- **`<filename>.dist`**: Hamming distance matrix
- **`<filename>.ripser`**: Ripser persistent homology output
- **`<filename>.timedist`**: time-transformed distance matrix (if `--timeseries` is used)