import argparse
from operator import index
from Bio import Seq,SeqIO
import pandas
import pandas as pd

def mask_sequence(sequence, mask_ranges):
    """Mask the specified ranges in the given sequence."""
    for mask_range in mask_ranges:
        if '-' in mask_range:
            start, end = map(int, mask_range.split('-'))
            sequence = sequence[:start-1] + 'N'*(end-start+1) + sequence[end:]
        else:
            pos = int(mask_range) - 1
            sequence = sequence[:pos] + 'N' + sequence[pos+1:]
    return sequence

def main(args):
    # Read mask file and create a dictionary of sequence ids and their mask ranges
    mask_dict = {}
    df = pd.read_csv(args.mask_file, sep='\t', index_col=False, header=None)
    for line in df.itertuples(index=False):
        print(line)
        seq_id, mask_range = line
        if seq_id not in mask_dict:
            mask_dict[seq_id] = []
        mask_dict[seq_id].append(mask_range)

    # Read the input alignment and mask sequences as specified
    with open(args.input_alignment, 'r') as infile, open(args.output_alignment, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in mask_dict:
                record.seq = Seq.Seq(mask_sequence(str(record.seq), mask_dict[record.id]))
            SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mask specific regions in sequence alignments.")
    parser.add_argument("--input-alignment", required=True, help="Path to the input FASTA file.")
    parser.add_argument("--output-alignment", required=True, help="Path to the output (masked) FASTA file.")
    parser.add_argument("--mask-file", required=True, help="Path to the TSV file specifying which regions to mask.")
    
    args = parser.parse_args()
    main(args)