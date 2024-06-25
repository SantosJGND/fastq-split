### Fastq-split

A simple script to split a fastq file into multiple files. The script will split the fastq file into a specified number of files. The number of reads per file is calculated on the fly.

The script is written to be applicable to pair end reads individually, retaining the order and number of sequences per subset file for both forward and reverse reads.

### Usage

```bash
python fastq_split.py [input fastq file] [output directory] [number of files] --compression_level [default 2]
```

The script is written to be as fast as possible, and relies for that purpose on the libraries:

- `xopen` for fast file reading and writing
- `dnaio` for fastq parsing
