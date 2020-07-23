# coding: utf-8

from __future__ import absolute_import

import argparse
import pysam
import os
import re
import sys

'''
command line example:
amplisep.py
    -f (--file)
    -o (--outdir)
    -d (--design)
'''

def parse_options():
    # Create parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file',
                        help='input BAM file',
                        type=str)
    parser.add_argument('-o', '--outdir',
                        help='output directory',
                        default='',
                        type=str)
    parser.add_argument('-d', '--design',
                        help='design BED file',
                        type=str)
    options = parser.parse_args()
    return options


def err_hand_cmd(ex):
    sys.stderr.write(ex.message)
    raise SystemExit(1)


def err_hand_pipe(ex):
    raise ex


def check_file(filename, eh):
    """
    Check that file exist
    """
    if not os.path.exists(filename):
        msg = "{0} not found\n".format(filename)
        eh(IOError(msg))


def check_path(pathname, eh):
    """
    Check that path exist
    """
    if not os.path.exists(pathname) or not os.path.isdir(pathname):
        msg = "{0} not found\n".format(pathname)
        eh(IOError(msg))


def check_extension(filename, ext, eh):
    """
    Check file extension
    """
    if not filename.endswith(ext):
        msg = "{0} have invalid extension\n".format(filename)
        eh(IOError(msg))


class Amplicon(object):
    def __init__(self, ampl):
        self.name = ampl[3]
        self.chr = ampl[0]
        self.start = int(ampl[1])
        self.end = int(ampl[2])
        self.pool = int(re.findall('Pool=(\d);', ampl[7])[0])
        self.length = int(ampl[2]) - int(ampl[1])
        self.length_fp = len(re.findall('ForwardPrimer=([ACGTacgt]+);', ampl[7])[0])
        self.length_rp = len(re.findall('ReversePrimer=([ACGTacgt]+)', ampl[7])[0])


def split_pools(design):
    f_dsgn = open(design, 'r')
    txt_regions = f_dsgn.readlines()[1:]
    f_dsgn.close()
    target_regions = [Amplicon(x.split('\t')) for x in txt_regions]
    pools = {}
    for region in target_regions:
        pool = region.pool
        if pool not in pools:
            pools[pool] = [region]
        else:
            pools[pool].append(region)
    return pools


def get_range(strand, ampl):
    portion = 0.3
    if strand: #reverse read
        return [range(ampl.start - ampl.length_fp, ampl.end - int(round(ampl.length*portion))), range(ampl.end, ampl.end + ampl.length_rp)]
    else: #nonreverse read
        return [range(ampl.start - ampl.length_fp, ampl.start), range(ampl.start + int(round(ampl.length*portion)), ampl.end + ampl.length_rp)]


def check_read(read, ampl):
    global outbam
    read_range = get_range(read.is_reverse, ampl)
    if read.reference_name == ampl.chr and \
                    read.reference_start in read_range[0] and \
                    read.reference_end in read_range[1]:
        #print 'write it!!!', ampl.chr, ampl.start
        outbam.write(read)


def main(fname, dsgn, outdir, errorhandler=err_hand_pipe):
    global outbam
    # Check that both files exist
    check_file(fname, errorhandler)
    check_file(dsgn, errorhandler)

    # Check that both files have right extension
    check_extension(fname, '.bam', errorhandler)
    check_extension(dsgn, '.bed', errorhandler)

    # Create output directory name
    # If --outdir not specified - create bam files in input directory
    # If --outdir has specified - check that path exists
    if outdir == '':
        outdir = os.path.dirname(fname)
    else:
        check_path(outdir, errorhandler)
    if not outdir.endswith('/'):
        outdir = outdir + '/'

    # Split AmpliSeq design per pool
    pools = split_pools(dsgn)

    inputbam = pysam.AlignmentFile(fname, "rb")
    nameshot = os.path.basename(os.path.splitext(fname)[0])
    pool_iter = 0

    for pool in pools:
        print pool
        outname = outdir + nameshot + '.pool' + str(pools.keys()[pool_iter]) + '.bam'
        pool_iter += 1
        outbam = pysam.AlignmentFile(outname, "wb", template=inputbam)
        for amp in pools[pool]:
            pileup = inputbam.fetch(region=amp.chr + ':' + str(amp.start) + '-' + str(amp.end))
            for read in pileup:
                check_read(read, amp)
        outbam.close()
    inputbam.close()

if __name__ == '__main__':
    options = parse_options()
    main(fname=options.file,
         dsgn=options.design,
         outdir=options.outdir,
         errorhandler=err_hand_cmd)
    pass
