# logicali

Usage:

```
    python scripts/generate_random_alignment.py 1 2000 500 aa
    python scripts/trim.py -i random_ali_seed1_l2000_n500_AA.fasta -n 100 --chi2 -o examples/outplots_aa --plot
```

```
    python scripts/generate_random_alignment.py 1 2000 500 nt
    python scripts/trim.py -i random_ali_seed1_l2000_n500_NT.fasta -n 100 --chi2 -o examples/outplots_aa --plot
```

benshmark:
 - 10k columns times 30k sequences, trimming and plots

```
    # ~4.5 min
    python scripts/generate_random_alignment.py 1 10000 30000 nt
    # ~9 min
    python scripts/trim.py -i random_ali_seed1_l10000_n30000_NT.fasta -n 10000 --chi2 -o examples/outplots_big_nt --plot
```


 - 10k columns times 50k sequences, trimming, no plots (plotting requires too much memory)
```
    # ~7 min
    python scripts/generate_random_alignment.py 1 10000 50000 nt
    # ~11 min
    python scripts/trim.py -i random_ali_seed1_l10000_n50000_NT.fasta -n 25000 --chi2
```
