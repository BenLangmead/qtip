#!/usr/bin/env python

"""
Dump all the CSED and CID info into a giant csv.
"""

import os


def parse_aligner_local(target):
    # Parsing aligner info
    toks = target.split('_')
    aligner, local = 'bt2', False
    if 'bwamem' in toks[1]:
        aligner = 'bwamem'
        local = True
    if 'snap' in toks[1]:
        aligner = 'snap'
        local = True
    if aligner == 'bt2' and toks[1][-1] == 'l':
        toks[1] = toks[1][:-1]
        local = True
    return aligner, local


def parse_species(target):
    genome = 'hg'
    toks = target.split('_')
    # Parsing reference species info
    if toks[2] == 'mm':
        genome = 'mm'
    if toks[2] == 'zm':
        genome = 'zm'
    return genome


def parse_sim(target):
    sim = 'mason'
    toks = target.split('_')
    # Parsing simulator info:
    if toks[2] == 'wgsim':
        sim = 'wgsim'
    if toks[2] == 'art':
        sim = 'art'
    return sim


def parse_sensitivity(target, aligner):
    toks = target.split('_')
    sensitivity = 's'
    if aligner == 'bt2':
        sensitivity = toks[1][3:]
    return sensitivity


def parse_paired(target):
    return target.startswith('r12')


def parse_readlen(target):
    toks = target.split('_')
    readlen = int(toks[-2])
    if readlen == '50to500':
        readlen = 500
    return readlen


def parse_name_and_target(combined):
    toks = combined.split('_')
    roff = 3
    if toks[2] == 'r0' or toks[2] == 'r12':
        roff = 2
    assert toks[roff] == 'r0' or toks[roff] == 'r12'
    return '_'.join(toks[:roff]), '_'.join(toks[roff:])


for fn in filter(os.path.isdir, os.listdir('.')):

    nm, target = parse_name_and_target(fn)

    paired = parse_paired(target)
    aligner, local = parse_aligner_local(target)
    genome = parse_species(target)
    sim = parse_sim(target)
    sensitivity = parse_sensitivity(target, aligner)
    readlen = parse_readlen(target)

    local = 'T' if local else 'F'

    for tt, training in [('training', 'T'), ('test', 'F')]:
        for cd, cid in [('cid', 'T'), ('csed', 'F')]:
            with open(os.path.join(fn, '_'.join([tt, cd, 'round.csv']))) as fh:
                for i, ln in enumerate(fh):
                    print ','.join([nm, 'T' if paired else 'F', aligner, local,
                                    sensitivity, genome, sim, str(readlen),
                                    training, cid, str(i), ln.rstrip()])
