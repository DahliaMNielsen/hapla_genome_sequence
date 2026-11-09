# hapla_genome_sequence

## extract_sequence.py

Extracts the DNA sequence around a nucleotide position in one contig of a FASTA file.
It needs only Python 3 (no extra packages), reads gzip-compressed files directly, and
streams the file line by line, so it works on whole genomes without loading them into memory.

```
./extract_sequence.py -f FASTA -c CONTIG -p POSITION -n FLANK [-w WIDTH]
```

| Option | Meaning |
|---|---|
| `-f`, `--fasta` | FASTA file (plain or `.gz`) |
| `-c`, `--contig` | Contig ID: the first word of its `>` header line |
| `-p`, `--position` | Nucleotide position, 1-based (first base = 1) |
| `-n`, `--flank` | Nucleotides to extract on each side of the position |
| `-w`, `--width` | Output line width; `0` prints the sequence on one line (default 60) |

The result is printed in FASTA format, `2 × FLANK + 1` nucleotides long:

```
$ ./extract_sequence.py -f genome.fasta -c contig_12 -p 15000 -n 10
>contig_12:14990-15010 position=15000 flank=10
ACGTTGCAAGCTTGACCGTAA
```

If the window runs past either end of the contig, it is cut short at that end and a
warning goes to stderr. A missing file, unknown contig, or position past the contig
end stops the program with an error and exit code 1.
