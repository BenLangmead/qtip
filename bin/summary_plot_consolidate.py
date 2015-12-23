#!/usr/bin/env python

import os

for fn in filter(os.path.isdir, os.listdir('.')):
    toks = fn.split('_')
    roff = 3
    if toks[2] == 'r0' or toks[2] == 'r12':
        roff = 2
    assert toks[roff] == 'r0' or toks[roff] == 'r12'
    nm = '_'.join(toks[:roff])
    paired = 'F' if toks[roff] == 'r0' else 'T'
    aligner = 'bt2'
    genome = 'hg'
    sim = 'mason'
    local = False
    if toks[roff+1].endswith(toks[-2]):
        toks[roff+1] = toks[roff+1][:-len(toks[-2])]

    # Parsing aligner info
    if 'bwamem' in toks[roff+1]:
        aligner = 'bwamem'
        local = True
    if 'snap' in toks[roff+1]:
        aligner = 'snap'
        local = True
    if aligner == 'bt2' and toks[roff+1][-1] == 'l':
        toks[roff+1] = toks[roff+1][:-1]
        local = True

    # Parsing reference species info
    if toks[roff+2] == 'mm':
        genome = 'mm'
    if toks[roff+2] == 'zm':
        genome = 'zm'

    # Parsing simulator info:
    if toks[roff+2] == 'wgsim':
        sim = 'wgsim'
    if toks[roff+2] == 'art':
        sim = 'art'

    # Parsing sensitivity info:
    sensitivity = 's'
    if aligner == 'bt2':
        sensitivity = toks[roff+1][3:]
    local = 'T' if local else 'F'
    readlen = toks[-2]
    if readlen == '50to500':
        readlen = '500'

    for tt, training in [('training', 'T'), ('test', 'F')]:
        for cd, cid in [('cid', 'T'), ('csed', 'F')]:
            with open(os.path.join(fn, '_'.join([tt, cd, 'round.csv']))) as fh:
                for i, ln in enumerate(fh):
                    print ','.join([nm, paired, aligner, local, sensitivity, genome, sim, readlen, training, cid, str(i), ln.rstrip()])
