"""
sam.py
"""

def cigarToList(cigar):
    """ Parse CIGAR string into a list of CIGAR operations. """
    ret = [];
    i = 0
    # CIGAR operations: MIDNSHP
    op_map = {'M':0, '=':0, 'X':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6}
    while i < len(cigar):
        run = 0
        while i < len(cigar) and cigar[i].isdigit():
            run *= 10
            run += int(cigar[i])
            i += 1
        assert i < len(cigar)
        op = cigar[i]
        i += 1
        assert op in op_map
        ret.append([op_map[op], run])
    return ret

def mdzToList(md):
    """ Parse MD:Z string into a list of operations, where 0=match,
        1=read gap, 2=mismatch. """
    i = 0;
    ret = [] # list of (op, run, str) tuples
    while i < len(md):
        if md[i].isdigit(): # stretch of matches
            run = 0
            while i < len(md) and md[i].isdigit():
                run *= 10
                run += int(md[i])
                i += 1 # skip over digit
            if run > 0:
                ret.append([0, run, ""])
        elif md[i].isalpha(): # stretch of mismatches
            mmstr = ""
            while i < len(md) and md[i].isalpha():
                mmstr += md[i]
                i += 1
            assert len(mmstr) > 0
            ret.append([1, len(mmstr), mmstr])
        elif md[i] == "^": # read gap
            i += 1 # skip over ^
            refstr = ""
            while i < len(md) and md[i].isalpha():
                refstr += md[i]
                i += 1 # skip over inserted character
            assert len(refstr) > 0
            ret.append([2, len(refstr), refstr])
        else: assert False
    return ret

def cigarMdzToStacked(seq, cgp, mdp_orig):
    """ Takes parsed cigar and parsed MD:Z and generates a stacked alignment:
        a pair of strings with gap characters inserted (possibly) and where
        characters at at the same offsets are opposite each other in the
        alignment.  Returns tuple of 2 parallel strings: read string, ref
        string.
        
        TODO: Fix so it works in local mode. """
    mdp = mdp_orig[:]
    rds, rfs = "", ""
    mdo, rdoff = 0, 0
    for c in cgp:
        op, run = c
        skipping = (op == 4)
        assert skipping or mdo < len(mdp), "%d, %d" % (mdo, len(mdp))
        if op == 0:   # M
            # Look for block matches and mismatches in MD:Z string
            mdrun = 0
            runleft = run
            while runleft > 0 and mdo < len(mdp):
                op_m, run_m, st_m = mdp[mdo]
                run_comb = min(runleft, run_m)
                runleft -= run_comb
                assert op_m == 0 or op_m == 1
                rds += seq[rdoff:rdoff + run_comb]
                if op_m == 0:
                    # match?
                    rfs += seq[rdoff:rdoff + run_comb]
                else:
                    #mismatch?
                    assert len(st_m) == run_comb
                    rfs += st_m
                    for i in xrange(0, run_comb):
                        assert st_m[i] != seq[rdoff+i:rdoff+i+1] or st_m[i] == 'N'
                mdrun += run_comb
                rdoff += run_comb
                # A stretch of matches in MD:Z could span M and I sections of
                # CIGAR
                if run_comb < run_m:
                    assert op_m == 0
                    mdp[mdo][1] -= run_comb
                else:
                    mdo += 1
        elif op == 1: # I
            rds += seq[rdoff:rdoff + run]
            rfs += "-" * run
            rdoff += run
        elif op == 2: # D
            op_m, run_m, st_m = mdp[mdo]
            assert op_m == 2
            assert run == run_m
            assert len(st_m) == run
            mdo += 1
            rds += "-" * run
            rfs += st_m
        elif op == 3: # N
            # If this is really long, this probably isn't the best way to do
            # this
            rds += "-" * run
            rfs += "-" * run
        elif op == 4: # S
            rdoff += run
        elif op == 5: # H
            pass
        elif op == 6: # P
            assert False
        elif op == 7: # =
            assert False
        elif op == 8: # X
            assert False
        else: assert False
    assert mdo == len(mdp)
    return rds, rfs
