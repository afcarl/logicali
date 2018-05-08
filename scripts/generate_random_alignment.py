#!/usr/bin/env python

from random import random, seed
import sys

seed_num = int(sys.argv[1])
lseqs = int(sys.argv[2])
nseqs = int(sys.argv[3])
stype = sys.argv[4].upper()

seed(seed_num)

if stype == 'NT':
    chars = 'ATGCATGCGCATGC--'
else:
    chars = 'GGGGGGGGGPPPPPAAAAAAAVVVVVVVVVVLLLLLIIMCCFFFFFYYYYWHKKKKKKKKRRQQNNNNNEEDDDSSSTT----------'

seqs = []
for i in xrange(lseqs):
    c = chars[int(random() * len(chars))]
    if i == 0 or i % 6:
        p = int(random() * 5) / 5.
    col = []
    for _ in xrange(nseqs):
        if random() < p:
            col.append(c)
        else:
            col.append(chars[int(random() * len(chars))])
    seqs.append(col)

out = open('random_ali_seed%d_l%d_n%d_%s.fasta' % (seed_num, lseqs, nseqs, stype), 'w')
out.write(''.join('>%d\n%s\n' % (i, ''.join(l)) for i, l in enumerate(zip(*seqs))))
out.close()
