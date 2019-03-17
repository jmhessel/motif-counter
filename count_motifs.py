'''
'''
import argparse
import collections
import itertools
import numpy as np
import scipy.stats
import tqdm
import prettytable


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
    parser.add_argument('--prior',
                        type=float,
                        default=100.,
                        help='How strong should the dirichlet prior be?')
    parser.add_argument('--pval',
                        type=float,
                        default=.0001,
                        help='What is the p-value for a true difference?')
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


def count_motifs(seqs, window):
    motif_counts = collections.defaultdict(int)
    print('Counting motifs...')
    for s in tqdm.tqdm(seqs):
        s = s.upper()
        for start_idx in range(len(s)-window+1):
            motif = s[start_idx:start_idx+window]
            motif_counts[motif] += 1
    return motif_counts


def z_score_to_pval(score):
    '''For a given Z-score, returns a 2-sided p value'''
    return scipy.stats.norm.sf(abs(score))*2


def main():
    args = parse_args()    
    vocab = list(itertools.product(
        args.character_vocabulary, repeat=args.window))
    vocab = set([''.join(v) for v in vocab])
    # make sure that the vocab size is what we expect
    assert(len(vocab) == len(args.character_vocabulary)**args.window)
    vocab2idx = {v:idx for idx, v in enumerate(sorted(list(vocab)))}

    if len(vocab) > 10**6:
        print('Warning: your vocab size is over 1M! Much larger and we might be slow...')
    
    print('There were {} motifs of size {}'.format(len(vocab), args.window))
    seqs1, seqs2 = map(load_sequences, [args.file1, args.file2])
    motifs1, motifs2 = map(lambda x: count_motifs(x, args.window), [seqs1, seqs2])

    print('{}/{} motifs in seqs1/seqs2.'.format(sum(motifs1.values()),
                                                sum(motifs2.values())))

    count_matrix = np.zeros([2, len(vocab)], dtype=np.float32)
    for v, idx in vocab2idx.items():
        count_matrix[0, idx] = motifs1[v]
        count_matrix[1, idx] = motifs2[v]

    n1, n2 = np.sum(count_matrix, axis=1)
        
    results = []
    prior, sum_prior = args.prior, len(vocab2idx) * args.prior
    for v, idx in vocab2idx.items():
        term1 = np.log((count_matrix[0, idx]+prior)/(n1+sum_prior-count_matrix[0, idx]-prior))
        term2 = np.log((count_matrix[1, idx]+prior)/(n2+sum_prior-count_matrix[1, idx]-prior))
        delta = term1 - term2
        var = 1./(count_matrix[0, idx] + prior) + 1./(count_matrix[1, idx] + prior)
        # motif, z-score, smoothed rate in 1, smoothed rate in 2
        results.append((v, delta/np.sqrt(var), np.exp(term1), np.exp(term2)))

    sig_results = [r for r in results if z_score_to_pval(r[1]) < args.pval]

    table1 = prettytable.PrettyTable()
    table1.field_names = ['Motif', 'Z-score', 'Rate in seqs1', 'Rate in seqs2']
    print('More in {} sequences:'.format(args.file1))
    for r in sorted(sig_results, key=lambda x: -x[1]):
        if r[1] > 0:
            row = [r[0], '{:.1f}'.format(r[1]), '{:.4f}'.format(r[2]), '{:.4f}'.format(r[3])]
            table1.add_row(row)
    print(table1)

    table2 = prettytable.PrettyTable()
    table2.field_names = ['Motif', 'Z-score', 'Rate in seqs1', 'Rate in seqs2']
    print('More in {} sequences:'.format(args.file2))
    for r in sorted(sig_results, key=lambda x: x[1]):
        if r[1] < 0:
            row = [r[0], '{:.1f}'.format(r[1]), '{:.4f}'.format(r[2]), '{:.4f}'.format(r[3])]
            table2.add_row(row)
    print(table2)
    
    
if __name__ == '__main__':
    main()
