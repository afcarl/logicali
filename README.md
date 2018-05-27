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
# ~4.5 min
python scripts/trim.py -i random_ali_seed1_l10000_n30000_NT.fasta -n 10000 --chi2 -o outplots_big_nt_chi2.fasta --plot
```
*Note for testing: size of output FASTA file is 233891093*.

 - 10k columns times 50k sequences, trimming, no plots (plotting requires too much memory)
```
# ~7 min
python scripts/generate_random_alignment.py 1 10000 50000 nt

# ~2 min 1 CPU
python scripts/trim.py -i random_ali_seed1_l10000_n50000_NT.fasta -n 25000 -C 1 -o outplots_huge_nt.fasta

# ~4 min 1 CPU
python scripts/trim.py -i random_ali_seed1_l10000_n50000_NT.fasta -n 25000 --chi2 -C 1 -o outplots_huge_nt_chi2.fasta

# ~40 sec 8 CPUs
python scripts/trim.py -i random_ali_seed1_l10000_n50000_NT.fasta -n 25000 --chi2 -C 8 -o outplots_huge_nt_chi2.fasta

# ~20 sec 8 CPUs
python scripts/trim.py -i random_ali_seed1_l10000_n50000_NT.fasta -n 25000 -C 8 -o outplots_huge_nt.fasta
```

*Note for testing: size of output FASTA files are 474429401 without `--chi2` paramenter and 375931371 with it*.
