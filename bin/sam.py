"""
sam.py
"""


class Cigar(object):
    """ Encapsulates a CIGAR string and the information it encodes about the
        shape of an alignment. """

    op_map = {'M': 0, '=': 0, 'X': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6}

    @staticmethod
    def cigar_to_list(cigar_string):
        ret = []
        i = 0
        while i < len(cigar_string):
            run = 0
            while i < len(cigar_string) and cigar_string[i].isdigit():
                run *= 10
                run += int(cigar_string[i])
                i += 1
            assert i < len(cigar_string)
            op = cigar_string[i]
            i += 1
            assert op in Cigar.op_map
            ret.append([Cigar.op_map[op], run])
        return ret

    def reference_length(self):
        """ Return the number of characters spanned by the reference side of
            the alignment. """
        tot = 0
        for op, run in self.cigar_list:
            if op == 0 or op == 2 or op == 3:
                tot += run
        return tot

    def read_length(self):
        """ Return the number of characters spanned by the read side of the
            alignment. """
        tot = 0
        for op, run in self.cigar_list:
            if op == 0 or op == 1 or op == 3:
                tot += run
        return tot

    def __init__(self, cigar_string):
        self.cigar_string = cigar_string
        self.cigar_list = self.cigar_to_list(cigar_string)


class Mdz(object):

    @staticmethod
    def mdz_to_list(mdz_string):
        i = 0
        ret = []  # list of (op, run, str) tuples
        while i < len(mdz_string):
            if mdz_string[i].isdigit():  # stretch of matches
                run = 0
                while i < len(mdz_string) and mdz_string[i].isdigit():
                    run *= 10
                    run += int(mdz_string[i])
                    i += 1  # skip over digit
                if run > 0:
                    ret.append([0, run, ""])
            elif mdz_string[i].isalpha():  # stretch of mismatches
                mmstr = ""
                while i < len(mdz_string) and mdz_string[i].isalpha():
                    mmstr += mdz_string[i]
                    i += 1
                assert len(mmstr) > 0
                ret.append([1, len(mmstr), mmstr])
            elif mdz_string[i] == "^":  # read gap
                i += 1  # skip over ^
                refstr = ""
                while i < len(mdz_string) and mdz_string[i].isalpha():
                    refstr += mdz_string[i]
                    i += 1  # skip over inserted character
                assert len(refstr) > 0
                ret.append([2, len(refstr), refstr])
            else:
                assert False
        return ret

    def __init__(self, mdz_string):
        self.mdz_string = mdz_string
        self.mdz_list = self.mdz_to_list(mdz_string)


def cigar_mdz_to_stacked(seq, cigar, mdz):
    """ Takes parsed cigar and parsed MD:Z and generates a stacked alignment:
        a pair of strings with gap characters inserted (possibly) and where
        characters at at the same offsets are opposite each other in the
        alignment.  Returns tuple of 2 parallel strings: read string, ref
        string.
        
        TODO: Fix so it works in local mode. """
    mdp = mdz.mdz_list[:]
    rds, rfs = "", ""
    mdo, rdoff = 0, 0
    for op, run in cigar.cigar_list:
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
        elif op == 1:  # I
            rds += seq[rdoff:rdoff + run]
            rfs += "-" * run
            rdoff += run
        elif op == 2:  # D
            op_m, run_m, st_m = mdp[mdo]
            assert op_m == 2
            assert run == run_m
            assert len(st_m) == run
            mdo += 1
            rds += "-" * run
            rfs += st_m
        elif op == 3:  # N
            # If this is really long, this probably isn't the best way to do
            # this
            rds += "-" * run
            rfs += "-" * run
        elif op == 4:  # S
            rdoff += run
        elif op == 5:  # H
            pass
        elif op == 6:  # P
            assert False
        elif op == 7:  # =
            assert False
        elif op == 8:  # X
            assert False
        else:
            assert False
    assert mdo == len(mdp)
    return rds, rfs


def cigar_ref_to_stacked(seq, cigar, ref, ref_id, ref_off):
    """ seq is w/r/t forward reference strand
        ref_off is the offset of the 1st reference character involved in the
        non-clipped portion of the alignment
    """
    rds, rfs = [], []
    ref_substr = ref.get(ref_id, ref_off, cigar.reference_length())
    seq_off, ref_off = 0, 0
    for op, run in cigar.cigar_list:
        assert op < 6
        if op == 4 or op == 5:  # S/H
            continue
        if op == 0:   # M/X/=
            # Step through corresponding reference characters one by one
            rds.append(seq[seq_off:seq_off+run])
            rfs.append(ref_substr[ref_off:ref_off+run])
            seq_off += run
            ref_off += run
        elif op == 1:  # I
            rds.append(seq[seq_off:seq_off+run])
            rfs.append('-' * run)
            seq_off += run
        elif op == 2:  # D
            rds.append('-' * run)
            rfs.append(ref_substr[ref_off:ref_off+run])
            ref_off += run
        elif op == 3:  # N
            rds.append('-' * run)
            rfs.append('-' * run)
    return ''.join(rds), ''.join(rfs)
