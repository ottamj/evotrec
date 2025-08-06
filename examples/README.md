# EVOtRec Examples

This directory contains example data for demonstrating EVOtRec functionality.

## Files

- **`example.fasta`**: Sample SARS-CoV-2 sequence alignment data with temporal information

## Running the Example

From the project root directory, run:

```sh
python evotrec.py examples/example.fasta "hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30|China" --timeseries
```

This will generate output files in the `examples/` directory:
- `example.csv` - tRI analysis results
- `example.dist` - Hamming distance matrix
- `example.ripser` - Ripser persistent homology output
- `example.timedist` - Time-transformed distance matrix

## Data Description

The example dataset contains 3,594 SARS-CoV-2 sequences spanning from December 30, 2019 to September 30, 2021 (641 days). Each sequence header contains temporal information in the format `|YYYY-MM-DD|`, making it suitable for time series analysis.

The reference sequence is from the original Wuhan strain: `hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30|China`.
