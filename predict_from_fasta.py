import argparse
import pickle
from collections import Counter
from pathlib import Path

from Bio import SeqIO
import numpy as np

from functions.other_datasets_porting import dict2tdm

BASEMAP = {'a': 0, 'c': 1, 'g': 2, 't': 3}


def numerize_sequences(sequences, length=115, pad='a'):
    """Convert sequences to numeric arrays accepted by the model."""
    lengths = [len(s) for s in sequences]
    most_common = Counter(lengths).most_common(1)[0][0]
    arr = np.zeros((len(sequences), length), dtype=np.int8)
    for idx, seq in enumerate(sequences):
        if len(seq) < most_common:
            seq = pad * (most_common - len(seq)) + seq
        arr[idx] = [BASEMAP.get(b, 0) for b in seq.lower()[:length]]
    return arr


def predict(fasta_path, model_path, treat_as='36N', length=115):
    """Return log10 promoter strengths for sequences in FASTA."""
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    seqs = [str(r.seq) for r in records]
    names = [r.id for r in records]

    with open(model_path, 'rb') as fh:
        model_dict = pickle.load(fh)
    tdm = dict2tdm(model_dict, treat_as=treat_as)

    num = numerize_sequences(seqs, length=length)
    bricks = tdm.sequences2bricks(num)
    pons = tdm.bricks2pons(bricks)
    return list(zip(names, np.log10(pons)))


def main():
    parser = argparse.ArgumentParser(description='Predict promoter strength from FASTA sequences.')
    parser.add_argument('fasta', help='Input FASTA file')
    parser.add_argument('model', help='Pickle file containing model parameters')
    parser.add_argument('--treat-as', default='36N', help='Dataset key inside the model file')
    parser.add_argument('--length', type=int, default=115, help='Sequence length for numerization')
    parser.add_argument('--output', type=Path, help='Optional output file to write predictions')
    args = parser.parse_args()

    results = predict(args.fasta, args.model, args.treat_as, args.length)

    out_lines = [f"{name}\t{strength:.4f}" for name, strength in results]
    if args.output:
        args.output.write_text("\n".join(out_lines))
    else:
        for line in out_lines:
            print(line)


if __name__ == '__main__':
    main()
