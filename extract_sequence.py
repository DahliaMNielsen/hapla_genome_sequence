#!/usr/bin/env python3
"""Extract a DNA sequence window around a nucleotide position from a FASTA file.

Positions are 1-based, as in most genome browsers and annotation formats.
The window covers FLANK nucleotides on each side of POSITION, so it is
2 * FLANK + 1 nucleotides long unless it is cut short at a contig end.

Example:
    extract_sequence.py -f genome.fasta -c contig_12 -p 15000 -n 50
"""

import argparse
import gzip
import sys


def open_fasta(path):
    """Open a FASTA file for reading text, gzip-compressed or not."""
    with open(path, "rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"
    if is_gzip:
        return gzip.open(path, "rt")
    return open(path, "r")


def extract_window(fasta_path, contig_id, position, flank):
    """Return (sequence, start, end, contig_length) for the requested window.

    start and end are the 1-based, inclusive coordinates actually extracted,
    which differ from the requested ones when the window runs past a contig
    end. Reads the file one line at a time, so whole genomes are never held
    in memory. Raises KeyError if the contig is not found.
    """
    want_start = max(1, position - flank)
    want_end = position + flank

    found = False
    seq_parts = []
    length = 0  # nucleotides of the contig read so far

    with open_fasta(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if found:
                    break  # reached the next record; the contig is complete
                # The contig ID is the first word of the header line.
                header_words = line[1:].split()
                found = bool(header_words) and header_words[0] == contig_id
                continue
            if not found:
                continue

            line = line.strip()
            line_start = length + 1
            length += len(line)
            # Keep only the part of this line that falls inside the window.
            if length >= want_start and line_start <= want_end:
                lo = max(want_start, line_start) - line_start
                hi = min(want_end, length) - line_start + 1
                seq_parts.append(line[lo:hi])

    if not found:
        raise KeyError(contig_id)

    sequence = "".join(seq_parts)
    end = min(want_end, length)
    return sequence, want_start, end, length


def wrap(sequence, width):
    """Split a sequence into lines of the given width (0 = no wrapping)."""
    if width <= 0:
        return sequence
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract the DNA sequence surrounding a nucleotide position "
                    "in a contig of a FASTA file.",
        epilog="Example: %(prog)s -f genome.fasta -c contig_12 -p 15000 -n 50",
    )
    parser.add_argument("-f", "--fasta", required=True,
                        help="FASTA file to read (may be gzip-compressed)")
    parser.add_argument("-c", "--contig", required=True,
                        help="contig ID: the first word of its '>' header line")
    parser.add_argument("-p", "--position", required=True, type=int,
                        help="1-based nucleotide position at the centre of the window")
    parser.add_argument("-n", "--flank", required=True, type=int,
                        help="number of nucleotides to extract on each side of POSITION")
    parser.add_argument("-w", "--width", type=int, default=60,
                        help="line width of the output sequence; 0 for one line "
                             "(default: 60)")
    args = parser.parse_args(argv)

    if args.position < 1:
        parser.error("--position must be 1 or greater")
    if args.flank < 0:
        parser.error("--flank must be 0 or greater")
    return args


def main(argv=None):
    args = parse_args(argv)

    try:
        sequence, start, end, length = extract_window(
            args.fasta, args.contig, args.position, args.flank)
    except FileNotFoundError:
        sys.exit(f"Error: file not found: {args.fasta}")
    except KeyError:
        sys.exit(f"Error: contig '{args.contig}' not found in {args.fasta}")

    if args.position > length:
        sys.exit(f"Error: position {args.position} is beyond the end of contig "
                 f"'{args.contig}' ({length} nt)")

    if start > args.position - args.flank or end < args.position + args.flank:
        print(f"Warning: window truncated at the contig end; extracted "
              f"{start}-{end} of contig '{args.contig}' ({length} nt)",
              file=sys.stderr)

    print(f">{args.contig}:{start}-{end} position={args.position} flank={args.flank}")
    print(wrap(sequence, args.width))


if __name__ == "__main__":
    main()
