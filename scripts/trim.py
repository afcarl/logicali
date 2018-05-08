#!/usr/bin/env python

from os import path, mkdir
import errno
from itertools import izip
from argparse import ArgumentParser


def plot_alignment(seqs, chars, sums, chi2, stds, prop, chi2_cut,
                   outdir=None, fname=''):
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

    ncols = min(1 + lseqs / 100, 12)
    nrows = min(1 + nseqs / 20 ,12)

    fig = plt.figure(figsize=(1.5 * ncols + 3, 1.5 * nrows + 5))

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

    if outdir:
        plt.savefig(path.join(outdir, 'alignment%s.png' % (('_' + fname) if fname else '')),
                    format='png')
    else:
        plt.show()


def main():
    opts = get_options()

    fname  = opts.ali
    cutoff = opts.ndata
    plot = opts.plot
    outdir = opts.outdir
    filt_chi2 = opts.chi2

    if outdir:
        try:
            mkdir(outdir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and path.isdir(outdir):
                pass
            else:
                raise

    fh = open(fname)

    ################################################################################
    # load data
    print ' - loading data...'
    seqs = []
    names = [fh.next().strip()]
    seq = ''
    for l in fh:
        if l.startswith('>'):
            if not len(names) % 10000:
                print '       loaded {:,} sequences'.format(len(names))
            names.append(l.strip())
            seqs.append(seq)
            seq = ''
        else:
            seq += l.strip()

    nseqs, lseqs = len(seqs), len(seqs[0])
    print '   * {:,} sequences with {:,} sites'.format(nseqs, lseqs)

    chars = set(s for l in seqs for s in l)
    chars = sorted([c for c in chars if c != '-'])

    chi2_cut = { 3:  5.9915,  4:  7.8147,  5:  9.4877,  6: 11.0705,  7: 12.5916,
                 8: 14.0671,  9: 15.5073, 10: 16.9190, 11: 18.3070, 12: 19.6751,
                13: 21.0261, 14: 22.3620, 15: 23.6848, 16: 24.9958, 17: 26.2962,
                18: 27.5871, 19: 28.8693, 20: 30.1435, 21: 31.4104, 22: 32.6706,
                23: 33.9244, 24: 35.1725, 25: 36.4150, 26: 37.6525, 27: 38.8851,
                28: 40.1133, 29: 41.3371, 30: 42.5570, 31: 43.7730, 32: 44.9853,
                33: 46.1943, 34: 47.3999, 35: 48.6024, 36: 49.8018, 37: 50.9985,
                38: 52.1923, 39: 53.3835, 40: 54.5722}[len(chars)]

    ################################################################################
    # keep only columns with data in at least a given number of sites
    print ' - removing sites with data in less than {} rows'.format(cutoff)
    col_cutoff = nseqs - cutoff

    good_cols = [i for i, col in enumerate(izip(*seqs)) if col.count('-') < col_cutoff]

    seqs = [[seqs[i][j] for j in good_cols] for i in xrange(nseqs)]

    print '   * kept {:,} of {:,} columns'.format(len(good_cols), lseqs)

    lseqs = len(seqs[0])

    ################################################################################
    # keep only columns with low complexity
    if filt_chi2:
        print ' - removing sites low complexity'
        col_cutoff = nseqs - cutoff

        expected = {}
        for c in chars:
            expected[c] = sum(s==c for l in seqs for s in l)
        total = float(sum(expected.values()))
        for c in expected:
            try:
                expected[c] /= total
            except ZeroDivisionError:
                expected[c] = 0

        print ('   * proportion of each site type: ' +
               ', '.join(['%s: %.4f' % (c.upper(), expected[c]) for c in chars]))

        good_cols = []
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
            av = [total  * expected[c] for c in chars]
            chi2.append(sum((c - av[i])**2 / av[i] for i, c in enumerate(count)))
            if chi2[-1] > chi2_cut:
                good_cols.append(i)
            prop.append([c / total for c in count])

        if plot:
            plot_alignment(seqs, chars, sums, chi2, stds, prop, chi2_cut, outdir=outdir, fname='filt1')
        seqs = [[seqs[i][j] for j in good_cols] for i in xrange(nseqs)]

        print '   * kept {:,} of {:,} columns'.format(len(good_cols), lseqs)


        lseqs = len(seqs[0])

    ################################################################################
    # mask lonely sites (surrounded by gaps)
    print ' - removing single sites surrounded by gaps'
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
    print ' - removing sites with data in less than {} rows'.format(cutoff)
    good_cols = [i for i, col in enumerate(izip(*seqs)) if col.count('-') < col_cutoff]

    print '   * kept {:,} of {:,} columns'.format(len(good_cols), lseqs)

    seqs = [[seqs[i][j] for j in good_cols] for i in xrange(nseqs)]

    lseqs = len(seqs[0])

    if plot:
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
        plot_alignment(seqs, chars, sums, chi2, stds, prop, chi2_cut, outdir=outdir, fname='final')

    ################################################################################
    # write result
    print ' - saving data'
    out = open(fname + '_trim.fasta', 'w')
    out.write(''.join('%s\n%s\n' % (names[i],
                                    ''.join(l))
                      for i, l in enumerate(seqs)))
    out.close()


    print '\nDone.'


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', '--input', dest='ali', metavar="PATH", action='store',
                        default=None, type=str,
                        help='path to alignment file (FASTA format)')
    parser.add_argument('-o', '--outdir', dest='outdir', metavar="PATH", action='store',
                        default=None, type=str,
                        help='path where to save figures')
    parser.add_argument('--plot', dest='plot',
                        action='store_true', default=False,
                        help=('generate interactive alignment plots with descriptive stats'))
    parser.add_argument('--chi2', dest='chi2',
                        action='store_true', default=False,
                        help=('filter out columns with random distribution '
                              'of sites with respect to background '
                              'proportions computed from input alignment'))
    parser.add_argument('-n', '--ndata', dest='ndata', metavar="INT", default=100,
                        type=int,
                        help='[%(default)s] Minimum number of sites with data per column')

    opts = parser.parse_args()
    return opts


if __name__ == '__main__':
    exit(main())
