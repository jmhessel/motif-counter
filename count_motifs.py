'''
'''
import argparse
import collections
import itertools
import numpy as np
import scipy.stats
import tqdm
import prettytable
import statsmodels.stats.multitest


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
                        default=10.,
                        help='How strong should the dirichlet prior be?')
    parser.add_argument('--pval',
                        type=float,
                        default=.001,
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
    seqs1, seqs2 = map(load_sequences, [args.file1, args.file2])
    motifs1, motifs2 = map(lambda x: count_motifs(x, args.window), [seqs1, seqs2])

    observed_vocab = set(list(motifs1.keys()) + list(motifs2.keys()))
    valid_vocab = []
    valid_chars = set(args.character_vocabulary)
    print('Checking each vocab item to make sure it only consists of {}'.format(
        args.character_vocabulary))
    for v in tqdm.tqdm(observed_vocab):
        if set(v).issubset(valid_chars):
            valid_vocab.append(v)
        else:
            continue
    
    print('Vocab size: {}'.format(len(valid_vocab)))
    print('{}/{} motifs in seqs1/seqs2.'.format(sum(motifs1.values()),
                                                sum(motifs2.values())))

    motifs1 = {v:motifs1[v] for v in valid_vocab}
    motifs2 = {v:motifs2[v] for v in valid_vocab}
    
    n1, n2 = sum(motifs1.values()), sum(motifs2.values())
        
    results = []
    prior, sum_prior = args.prior, len(valid_vocab) * args.prior
    print('Computing z-scores')
    
    for v in tqdm.tqdm(valid_vocab):
        term1 = np.log((motifs1[v] + prior)/(n1 + sum_prior - motifs1[v] - prior))
        term2 = np.log((motifs2[v] + prior)/(n2 + sum_prior - motifs2[v] - prior))
        delta = term1 - term2
        var = 1./(motifs1[v] + prior) + 1./(motifs2[v] + prior)
        # motif, z-score, smoothed rate in 1, smoothed rate in 2
        results.append((v, delta/np.sqrt(var), np.exp(term1), np.exp(term2)))

    # vectorized p_value computation, and multiple comparison correction
    z_scores = np.array([r[1] for r in results])
    p_vals = z_score_to_pval(z_scores)

    reject, p_val, _, _ = statsmodels.stats.multitest.multipletests(
        p_vals, alpha=args.pval, method='Holm')
    sig_results = [res for res, rej in zip(results, reject) if rej]
    print('{} significant results'.format(len(sig_results)))

    row_count = 0
    table1 = prettytable.PrettyTable()
    table1.field_names = ['Motif', 'Z-score', 'Rate in seqs1', 'Rate in seqs2']
    print('More in {} sequences:'.format(args.file1))
    for r in sorted(sig_results, key=lambda x: -x[1]):
        if r[1] > 0:
            row_count += 1
            row = [r[0], '{:.1f}'.format(r[1]), '{:.4f}'.format(r[2]), '{:.4f}'.format(r[3])]
            table1.add_row(row)
    
    if row_count == 0:
        print('None.')
    else:
        print(table1)
    
    row_count = 0
    table2 = prettytable.PrettyTable()
    table2.field_names = ['Motif', 'Z-score', 'Rate in seqs1', 'Rate in seqs2']
    print('More in {} sequences:'.format(args.file2))
    for r in sorted(sig_results, key=lambda x: x[1]):
        if r[1] < 0:
            row_count += 1
            row = [r[0], '{:.1f}'.format(r[1]), '{:.4f}'.format(r[2]), '{:.4f}'.format(r[3])]
            table2.add_row(row)
    if row_count == 0:
        print('None.')
    else:
        print(table2)
    
    
if __name__ == '__main__':
    main()
