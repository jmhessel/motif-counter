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
1 significant results
More in Jc_000226F.fa sequences:
None.
More in JC_optix300kb.fa sequences:
+-------+---------+---------------+---------------+
| Motif | Z-score | Rate in seqs1 | Rate in seqs2 |
+-------+---------+---------------+---------------+
| TTTTT |   -4.9  |     0.0086    |     0.0096    |
+-------+---------+---------------+---------------+
```


to see all command line options (e.g., for changing window size)

```
python count_motifs.py --help
```
