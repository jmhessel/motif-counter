# motif-counter

Helping with Martik's project, counting genetic motifs (with python3!)

## How to install:

This script requires python3, and several python packages. After
installing python3 and pip, you can do this to get the required
packages:

```
pip3 install -r requirements.txt --user
```

## How to run:

to compare Jc_000226F.fa to JC_optix300kb.fa with a window size of 5:

```
python count_motifs.py Jc_000226F.fa JC_optix300kb.fa  --window 5
```

the output I get is:

```
More in Jc_000226F.fa sequences:
+-------+---------+---------------+---------------+
| Motif | Z-score | Rate in seqs1 | Rate in seqs2 |
+-------+---------+---------------+---------------+
| AAACA |   4.8   |     0.0027    |     0.0022    |
| ACAAA |   4.2   |     0.0028    |     0.0024    |
+-------+---------+---------------+---------------+
More in JC_optix300kb.fa sequences:
+-------+---------+---------------+---------------+
| Motif | Z-score | Rate in seqs1 | Rate in seqs2 |
+-------+---------+---------------+---------------+
| CTGTC |   -4.4  |     0.0008    |     0.0011    |
| GCCTG |   -4.3  |     0.0005    |     0.0007    |
| CCTGT |   -4.2  |     0.0006    |     0.0009    |
| TGTCT |   -3.9  |     0.0011    |     0.0013    |
+-------+---------+---------------+---------------+
```


to see all command line options (e.g., for changing window size)

```
python count_motifs.py --help
```
