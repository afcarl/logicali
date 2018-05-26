#!/usr/bin/env python

from itertools import izip
from argparse  import ArgumentParser
from time      import time
from datetime  import datetime as dt
from sys import stdout

# generated using scipy.stats.chi2
CHI2_CUT = {'0.00001': {3: 23.0259, 4: 25.9017, 5: 28.4733, 6: 30.8562, 7: 33.1071, 8: 35.2585, 9: 37.3316, 10: 39.3407, 11: 41.2962, 12: 43.2060, 13: 45.0761, 14: 46.9116, 15: 48.7161, 16: 50.4930, 17: 52.2450, 18: 53.9743, 19: 55.6829, 20: 57.3725, 21: 59.0446, 22: 60.7003, 23: 62.3410, 24: 63.9675, 25: 65.5808, 26: 67.1818, 27: 68.7710, 28: 70.3492, 29: 71.9170, 30: 73.4749},
            '0.00005': {3: 19.8070, 4: 22.5547, 5: 25.0133, 6: 27.2937, 7: 29.4497, 8: 31.5124, 9: 33.5017, 10: 35.4312, 11: 37.3107, 12: 39.1476, 13: 40.9476, 14: 42.7153, 15: 44.4544, 16: 46.1678, 17: 47.8581, 18: 49.5274, 19: 51.1775, 20: 52.8100, 21: 54.4263, 22: 56.0276, 23: 57.6149, 24: 59.1891, 25: 60.7512, 26: 62.3017, 27: 63.8415, 28: 65.3712, 29: 66.8912, 30: 68.4022},
            '0.0001' : {3: 18.4207, 4: 21.1075, 5: 23.5127, 6: 25.7448, 7: 27.8563, 8: 29.8775, 9: 31.8276, 10: 33.7199, 11: 35.5640, 12: 37.3670, 13: 39.1344, 14: 40.8707, 15: 42.5793, 16: 44.2632, 17: 45.9249, 18: 47.5664, 19: 49.1894, 20: 50.7955, 21: 52.3860, 22: 53.9620, 23: 55.5246, 24: 57.0746, 25: 58.6130, 26: 60.1403, 27: 61.6573, 28: 63.1645, 29: 64.6624, 30: 66.1517},
            '0.0005' : {3: 15.2018, 4: 17.7300, 5: 19.9974, 6: 22.1053, 7: 24.1028, 8: 26.0178, 9: 27.8680, 10: 29.6658, 11: 31.4198, 12: 33.1366, 13: 34.8213, 14: 36.4778, 15: 38.1094, 16: 39.7188, 17: 41.3081, 18: 42.8792, 19: 44.4338, 20: 45.9731, 21: 47.4985, 22: 49.0108, 23: 50.5111, 24: 52.0002, 25: 53.4788, 26: 54.9475, 27: 56.4069, 28: 57.8576, 29: 59.3000, 30: 60.7346},
            '0.001'  : {3: 13.8155, 4: 16.2662, 5: 18.4668, 6: 20.5150, 7: 22.4577, 8: 24.3219, 9: 26.1245, 10: 27.8772, 11: 29.5883, 12: 31.2641, 13: 32.9095, 14: 34.5282, 15: 36.1233, 16: 37.6973, 17: 39.2524, 18: 40.7902, 19: 42.3124, 20: 43.8202, 21: 45.3147, 22: 46.7970, 23: 48.2679, 24: 49.7282, 25: 51.1786, 26: 52.6197, 27: 54.0520, 28: 55.4760, 29: 56.8923, 30: 58.3012},
            '0.005'  : {3: 10.5966, 4: 12.8382, 5: 14.8603, 6: 16.7496, 7: 18.5476, 8: 20.2777, 9: 21.9550, 10: 23.5894, 11: 25.1882, 12: 26.7568, 13: 28.2995, 14: 29.8195, 15: 31.3193, 16: 32.8013, 17: 34.2672, 18: 35.7185, 19: 37.1565, 20: 38.5823, 21: 39.9968, 22: 41.4011, 23: 42.7957, 24: 44.1813, 25: 45.5585, 26: 46.9279, 27: 48.2899, 28: 49.6449, 29: 50.9934, 30: 52.3356},
            '0.01'   : {3:  9.2103, 4: 11.3449, 5: 13.2767, 6: 15.0863, 7: 16.8119, 8: 18.4753, 9: 20.0902, 10: 21.6660, 11: 23.2093, 12: 24.7250, 13: 26.2170, 14: 27.6882, 15: 29.1412, 16: 30.5779, 17: 31.9999, 18: 33.4087, 19: 34.8053, 20: 36.1909, 21: 37.5662, 22: 38.9322, 23: 40.2894, 24: 41.6384, 25: 42.9798, 26: 44.3141, 27: 45.6417, 28: 46.9629, 29: 48.2782, 30: 49.5879},
            '0.025'  : {3:  7.3778, 4:  9.3484, 5: 11.1433, 6: 12.8325, 7: 14.4494, 8: 16.0128, 9: 17.5345, 10: 19.0228, 11: 20.4832, 12: 21.9200, 13: 23.3367, 14: 24.7356, 15: 26.1189, 16: 27.4884, 17: 28.8454, 18: 30.1910, 19: 31.5264, 20: 32.8523, 21: 34.1696, 22: 35.4789, 23: 36.7807, 24: 38.0756, 25: 39.3641, 26: 40.6465, 27: 41.9232, 28: 43.1945, 29: 44.4608, 30: 45.7223},
            '0.05'   : {3:  5.9915, 4:  7.8147, 5:  9.4877, 6: 11.0705, 7: 12.5916, 8: 14.0671, 9: 15.5073, 10: 16.9190, 11: 18.3070, 12: 19.6751, 13: 21.0261, 14: 22.3620, 15: 23.6848, 16: 24.9958, 17: 26.2962, 18: 27.5871, 19: 28.8693, 20: 30.1435, 21: 31.4104, 22: 32.6706, 23: 33.9244, 24: 35.1725, 25: 36.4150, 26: 37.6525, 27: 38.8851, 28: 40.1133, 29: 41.3371, 30: 42.5570},
            '0.1'    : {3:  4.6052, 4:  6.2514, 5:  7.7794, 6:  9.2364, 7: 10.6446, 8: 12.0170, 9: 13.3616, 10: 14.6837, 11: 15.9872, 12: 17.2750, 13: 18.5493, 14: 19.8119, 15: 21.0641, 16: 22.3071, 17: 23.5418, 18: 24.7690, 19: 25.9894, 20: 27.2036, 21: 28.4120, 22: 29.6151, 23: 30.8133, 24: 32.0069, 25: 33.1962, 26: 34.3816, 27: 35.5632, 28: 36.7412, 29: 37.9159, 30: 39.0875},
            '0.15'   : {3:  3.7942, 4:  5.3170, 5:  6.7449, 6:  8.1152, 7:  9.4461, 8: 10.7479, 9: 12.0271, 10: 13.2880, 11: 14.5339, 12: 15.7671, 13: 16.9893, 14: 18.2020, 15: 19.4062, 16: 20.6030, 17: 21.7931, 18: 22.9770, 19: 24.1555, 20: 25.3289, 21: 26.4976, 22: 27.6620, 23: 28.8225, 24: 29.9792, 25: 31.1325, 26: 32.2825, 27: 33.4295, 28: 34.5736, 29: 35.7150, 30: 36.8538},
            '0.2'    : {3:  3.2189, 4:  4.6416, 5:  5.9886, 6:  7.2893, 7:  8.5581, 8:  9.8032, 9: 11.0301, 10: 12.2421, 11: 13.4420, 12: 14.6314, 13: 15.8120, 14: 16.9848, 15: 18.1508, 16: 19.3107, 17: 20.4651, 18: 21.6146, 19: 22.7595, 20: 23.9004, 21: 25.0375, 22: 26.1711, 23: 27.3015, 24: 28.4288, 25: 29.5533, 26: 30.6752, 27: 31.7946, 28: 32.9117, 29: 34.0266, 30: 35.1394}}


def printime(msg):
    print (msg +
           (' ' * (79 - len(msg.split('\n')[-1]))) +
           '[' +
           str(dt.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']' + ('\n' * msg.endswith('\n')))


def get_molecules(seqs, molecule):
    chars = [s for l in seqs[:max(1, 10000 / len(seqs[0]))] for s in l[:10000]]
    if molecule == 'dna' or (molecule == 'auto' and sum(chars.count(l) for l in 'ACGTN-') == len(chars)):
        printime('\n - %s nucleotides' % ('define' if molecule != 'auto' else 'detected'))
        chars = ['A', 'C', 'G', 'T']
    else:
        printime('\n - %s amino acids' % ('define' if molecule != 'auto' else 'detected'))
        chars = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'] # 'U', 'O'
    return chars


def parse_fasta(fname):
    fh = open(fname)
    seqs = []
    names = [fh.next().strip()]
    seq = ''
    for l in fh:
        if l.startswith('>'):
            if not len(names) % 10000:
                stdout.write('       loaded {:,} sequences\r'.format(len(names)))
                stdout.flush()
            names.append(l.strip())
            seqs.append(seq)
            seq = ''
        else:
            seq += l.strip().upper()
    stdout.write('')
    return seqs, names


def plot_alignment(seqs, chars, sums, chi2, stds, prop, chi2_cut,
                   savefig=None, nox=False):
    if nox:
        import matplotlib
        matplotlib.use('Agg')

    from matplotlib import pyplot as plt, colors

    if len(chars) == 4:
        cmap = colors.LinearSegmentedColormap.from_list("",
                                                        ["#d62728", "#1f77b4",
                                                         "#ff7f0e", "#2ca02c"],
                                                        N=len(chars))
    else:
        cmap = plt.cm.get_cmap('tab20b')

    nseqs, lseqs = len(seqs), len(seqs[0])

    values = dict((c, i) for i, c in enumerate(chars, 1))
    values['-'] = float('nan')
    values['N'] = float('nan')
    values['X'] = float('nan')

    ncols = min(1 + lseqs / 100, 12)
    nrows = min(1 + nseqs / 20 ,12)

    _ = plt.figure(figsize=(1.5 * ncols + 3, 1.5 * nrows + 5))

    ax1 = plt.subplot2grid((nrows + 4, ncols + 1), (0, 0), rowspan=nrows, colspan=ncols)
    im = ax1.imshow([[values[s] for s in l] for l in seqs], aspect='auto', cmap=cmap)

    ax2 = plt.subplot2grid((nrows + 4, ncols + 1), (nrows, 0), sharex=ax1, colspan=ncols)
    ax2.plot(sums)
    ax2.set_title('count no gaps')

    ax3 = plt.subplot2grid((nrows + 4, ncols + 1), (nrows + 1, 0), sharex=ax1, colspan=ncols)
    ax3.plot(stds)
    ax3.set_yscale("log")
    ax3.axhline(0.05, color='k', lw=1, ls='--')
    ax3.grid()
    ax3.set_title('log std')

    colors = [cmap(i) for i in xrange(len(chars))]
    ax4 = plt.subplot2grid((nrows + 4, ncols + 1), (nrows + 2, 0), sharex=ax1, colspan=ncols)
    prop = zip(*prop)
    for i in xrange(len(prop)):
        ax4.plot(prop[i], color=colors[i])
    ax4.set_yticks([0.1, 0.25, 0.40])
    ax4.grid()
    ax4.set_title('proportion of sites')

    ax5 = plt.subplot2grid((nrows + 4, ncols + 1), (nrows + 3, 0), sharex=ax1, colspan=ncols)
    ax5.plot(chi2)
    ax5.grid()
    ax5.set_title('Chi2')
    ax5.axhline(chi2_cut, color='k', lw=1, ls='--')
    ax5.set_yscale("log")

    ax1.set_xlim(-0.5, lseqs + 0.5)

    axcolor = plt.subplot2grid((nrows + 4, ncols + 1), (max(0, nrows - 3), ncols), rowspan=min(3, max(1, nrows)))
    fact = len(chars) / (len(chars) + 1.)
    plt.colorbar(im, cax=axcolor, ticks=[ i * fact + 0.5 * fact
                                          for i in range(len(chars) + 1)])
    axcolor.set_yticklabels(chars)

    plt.tight_layout()

    if savefig:
        plt.savefig(savefig + '.png', format='png')
    else:
        plt.show()


def main():
    opts = get_options()

    fname  = opts.ali
    cutoff = opts.ndata
    plot = opts.plot
    output = opts.output
    filt_chi2 = opts.chi2

    print '\n Logicali trimming sequences...'

    ################################################################################
    # load data
    printime('\n - loading data...')
    seqs, names = parse_fasta(fname)

    nseqs, lseqs = len(seqs), len(seqs[0])
    printime('   * {:,} sequences with {:,} sites'.format(nseqs, lseqs))

    # find out if nucleotides or AA, checking first ~1000 chars
    chars = get_molecules(seqs, opts.molecule)

    if filt_chi2 is not None:
        try:
            chi2_cut = CHI2_CUT[filt_chi2][len(chars)]
        except KeyError:
            try:
                from scipy.stats import chi2
            except ImportError:
                raise Exception('ERROR: scipy not installed and p-value for '
                                'chi2 test not in pre-computed table (%s)' % (
                                    ', '.join(CHI2_CUT.keys())
                                ))
            printime(' - getting chi2 value with scipy')
            chi2_cut = chi2.isf(float(filt_chi2), len(chars) - 1)

    ################################################################################
    # keep only columns with data in at least a given number of sites
    printime('\n - removing sites with data in less than {} rows'.format(cutoff))
    col_cutoff = nseqs - cutoff

    good_cols = [i for i, col in enumerate(izip(*seqs)) if col.count('-') < col_cutoff]
    # TODO: TRY: instead of storing good_cols directly store cols, use them for chi2,and then zip them back
    # OR/AND: parallelize
    seqs = [[seqs[i][j] for j in good_cols] for i in xrange(nseqs)]

    printime('   * kept {:,} of {:,} columns'.format(len(good_cols), lseqs))

    lseqs = len(seqs[0])

    ################################################################################
    # keep only columns with low complexity
    if filt_chi2 is not None:
        printime(('\n - removing sites low complexity\n        '
                  '(random distribution: chi2 test p-value < %s)') % filt_chi2)
        col_cutoff = nseqs - cutoff

        expected = {}
        for c in chars[:]:
            expected[c] = sum(s==c for l in seqs for s in l)
            if not expected[c]:
                chars.remove(c)
                del(expected[c])

        total = float(sum(expected.values()))
        for c in expected:
            try:
                expected[c] /= total
            except ZeroDivisionError:
                expected[c] = 0

        printime('   * proportion of each site type: ' +
                 ', '.join(['%s: %.4f' % (c.upper(), expected[c]) for c in chars]))

        good_cols = []
        mean = 1. / len(chars)
        sums = []
        prop = []
        stds = []
        chi2 = []
        for i, col in enumerate(izip(*seqs)):
            count = [col.count(c) for c in chars]
            total = float(sum(count)) or 1  # in case total is equal to 0
            std = sum((c / total - mean)**2 for c in count)
            stds.append(std**0.5)
            sums.append(total)
            av = [total  * expected[c] for c in chars]
            chi2.append(sum((c - av[i])**2 / av[i] for i, c in enumerate(count)))
            if chi2[-1] > chi2_cut:
                good_cols.append(i)
            prop.append([c / total for c in count])

        if plot:
            printime('   * plotting')
            plot_alignment(seqs, chars, sums, chi2, stds, prop, chi2_cut,
                           savefig=output + '_filt1', nox=opts.nox)
        seqs = [[seqs[i][j] for j in good_cols] for i in xrange(nseqs)]

        printime('   * kept {:,} of {:,} columns'.format(len(good_cols), lseqs))

        lseqs = len(seqs[0])

    ################################################################################
    # mask lonely sites (surrounded by gaps)
    if opts.mask_lonely:
        printime('\n - removing single sites surrounded by gaps')
        for l in seqs:
            # rule for first and last
            if l[1] == '-':
                l[0] = '-'
            if l[-2] == '-':
                l[-1] = '-'
            for i in xrange(lseqs - 2):
                if l[i] == '-':
                    if l[i+2] == '-':
                        l[i+1] = '-'
                    elif i+3 < lseqs and l[i+3] == '-':
                        l[i+1] = '-'
                        l[i+2] = '-'

        ################################################################################
        # keep only columns with data in at least a given number of sites (second round)
        printime('\n - removing sites with data in less than {} rows'.format(cutoff))
        good_cols = [i for i, col in enumerate(izip(*seqs)) if col.count('-') < col_cutoff]

        printime('   * kept {:,} of {:,} columns'.format(len(good_cols), lseqs))

        seqs = [[seqs[i][j] for j in good_cols] for i in xrange(nseqs)]

        lseqs = len(seqs[0])

    if plot:
        printime('   * plotting')
        mean = 1. / len(chars)
        sums = []
        prop = []
        stds = []
        chi2 = []
        for i, col in enumerate(izip(*seqs)):
            count = [col.count(c) for c in chars]
            total = float(sum(count))
            std = sum((c / total - mean)**2 for c in count)
            stds.append(std**0.5)
            sums.append(total)
            if filt_chi2:
                av = [total  * expected[c] for c in chars]
                chi2.append(sum((c - av[i])**2 / av[i] for i, c in enumerate(count)))
            prop.append([c / total for c in count])
        plot_alignment(seqs, chars, sums, chi2, stds, prop, chi2_cut,
                       savefig=output + '_filt2')

    ################################################################################
    # write result
    print('\n - saving data')
    out = open(output, 'w')
    out.write(''.join('%s\n%s\n' % (names[i],
                                    ''.join(l))
                      for i, l in enumerate(seqs)))
    out.close()

    printime('\nDone.\n')


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', '--input', dest='ali', metavar="PATH", action='store',
                        default=None, type=str, required=True,
                        help='path to alignment file (FASTA format)')
    parser.add_argument('-o', '--output', dest='output', metavar="PATH", action='store',
                        default=None, type=str, required=True,
                        help='output file name, figures will use this name as prefix.')
    parser.add_argument('--plot', dest='plot',
                        action='store_true', default=False,
                        help=('generate interactive alignment plots with descriptive stats.'))
    parser.add_argument('--molecule', dest='molecule',
                        default='auto', choices=['auto', 'dna',  'protein'],
                        help=('[%(default)s] sequence type, by default inferred from input.'))
    parser.add_argument('--nox', dest='nox',
                        action='store_true', default=False,
                        help=('no X available for plotting'))
    parser.add_argument('--chi2', dest='chi2', nargs='?', default=None, const='0.05',
                        help=(('[%(default)s] filter out columns with random '
                              'distribution of sites with respect to background '
                              'proportions computed from input alignment. '
                              'Alternative p-value can be passed, the lower the more '
                               'columns filtered out (pre-computed: {}).').format(
                                  str(sorted(CHI2_CUT.keys(), key=lambda x: float(x)))
                              )))
    parser.add_argument('--mask_lonely', dest='mask_lonely',
                        action='store_true', default=False,
                        help=('treat as gaps sites surrounded by gaps '))
    parser.add_argument('-n', '--ndata', dest='ndata', metavar="INT", default=100,
                        type=int,
                        help='[%(default)s] Minimum number of sites with data per column.')

    opts = parser.parse_args()
    return opts


if __name__ == '__main__':
    exit(main())


# #    size1     size2                                    pH    Hydrophobicity
# A    89.10     71.08     2.34     9.69     -          6.00   47
# R   174.20    156.19     2.17     9.04     12.48     10.76  -26
# N   132.12    114.11     2.02     8.80     -          5.41  -41
# D   133.11    115.09     1.88     9.60     3.65       2.77  -18
# C   121.16    103.15     1.96     10.28     8.18      5.07   52
# E   147.13    129.12     2.19     9.67     4.25       3.22    8
# Q   146.15    128.13     2.17     9.13     -          5.65  -18
# G    75.07     57.05     2.34     9.60     -          5.97    0
# H   155.16    137.14     1.82     9.17     6.00       7.59  -42
# O   131.13    113.11     1.82     9.65     -          -       -
# I   131.18    113.16     2.36     9.60     -          6.02  100
# L   131.18    113.16     2.36     9.60     -          5.98  100
# K   146.19    128.18     2.18     8.95     10.53      9.74  -37
# M   149.21    131.20     2.28     9.21     -          5.74   74
# F   165.19    147.18     1.83     9.13     -          5.48   92
# P   115.13     97.12     1.99     10.6     -          6.30  -46
# U   139.11    121.09     -        -        -          5.68    -
# S   105.09     87.08     2.21     9.15     -          5.68   -7
# T   119.12    101.11     2.09     9.10     -          5.60   13
# W   204.23    186.22     2.83     9.39     -          5.89   84
# Y   181.19    163.18     2.20     9.11     10.07      5.66   49
# V   117.15     99.13     2.32     9.62     -          5.96   79
