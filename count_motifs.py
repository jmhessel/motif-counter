'''
'''
import argparse
import collections
import itertools
import numpy as np
import tqdm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('file1',
                        help='first file to compare')
    parser.add_argument('file2',
                        help='first file to compare')
    parser.add_argument('--window',
                        type=int,
                        default=5,
                        help='Size of the motif window')
    parser.add_argument('--character_vocabulary',
                        type=str,
                        default='CGAT',
                        help='What is the character set?')
    return parser.parse_args()


def load_sequences(fname):
    '''Loads a single sequence from a file'''
    all_seqs = []
    cur_seq = ''
    with open(fname) as f:
        for line in tqdm.tqdm(f):
            if len(line.strip()) == 0: continue
            if line[0] == '>':
                all_seqs.append(cur_seq)
                cur_seq = ''
                continue
            cur_seq += line.strip()
    all_seqs.append(cur_seq)
    print('Loaded {} seqs'.format(len(all_seqs)))
    return all_seqs


def count_motifs(seqs):
    collections.defaultdict(int)
    for s in seqs:
        pass


def main():
    args = parse_args()
    vocab = list(itertools.product(
        args.character_vocabulary, repeat=args.window))
    vocab = set([''.join(v) for v in vocab])
    # make sure that the vocab size is what we expect
    assert(len(vocab) == len(args.character_vocabulary)**args.window)
    print('There were {} motifs of size {}'.format(len(vocab), args.window))
    seqs1, seqs2 = map(load_sequences, [args.file1, args.file2])
    
    
    
    
if __name__ == '__main__':
    main()
