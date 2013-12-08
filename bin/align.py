'''
align.py

Functions for simple pairwise alignment tasks.
'''

try:
    from numpypy import zeros
except ImportError: pass

from numpy import zeros

def hammingDistance(x, y):
    ''' Return Hamming distance between 2 same-length strings '''
    assert len(x) == len(y)
    nmm = 0
    for i in xrange(0, len(x)):
        if x[i] != y[i]:
            nmm += 1
    return nmm

def boundEditDistance(x, y):
    ''' Return lower and upper bounds on the edit distance between two
        strings with sub-O(mn) work. '''
    if x == y: return 0, 0
    absDiff = abs(len(x) - len(y))
    # if lengths differ, need at least some deletions/insertions are
    # required to make lengths the same
    lower = max(1, absDiff)
    # an upper bound is the hamming distance between the shorter string
    # and a same-length prefix of the longer string, plus the number of
    # deletions/insertions needed to make them the same length
    minlen = min(len(x), len(y))
    upper = hammingDistance(y[:minlen], x[:minlen]) + absDiff
    if absDiff > 0:
        upper = min(upper, hammingDistance(y[-minlen:], x[-minlen:]) + absDiff)
    return lower, upper

def traceTranscript(D, x, y):
    ''' Trace back from bottom-right cell and report edit
        transcript. '''
    i, j, hi = len(x), len(y), max(len(x), len(y)) + 1
    xscript = []
    while i > 0 or j > 0:
        diag, vert, horz = hi, hi, hi
        delt = None
        if i > 0 and j > 0:
            delt = 0 if x[i-1] == y[j-1] else 1
            diag = D[i-1, j-1] + delt
        if i > 0:
            vert = D[i-1, j] + 1
        if j > 0:
            horz = D[i, j-1] + 1
        if diag <= vert and diag <= horz:
            xscript.append('R' if delt == 1 else 'M')
            i -= 1; j -= 1
        elif vert <= horz:
            xscript.append('I')
            i -= 1
        else:
            xscript.append('D')
            j -= 1
    return (''.join(xscript))[::-1]

def traceStacked(D, x, y):
    ''' Trace back from bottom-right cell and report stacked
        alignment. '''
    i, j, hi = len(x), len(y), max(len(x), len(y)) + 1
    xst, yst = [], []
    while i > 0 or j > 0:
        diag, vert, horz = hi, hi, hi
        delt = None
        if i > 0 and j > 0:
            delt = 0 if x[i-1] == y[j-1] else 1
            diag = D[i-1, j-1] + delt
        if i > 0:
            vert = D[i-1, j] + 1
        if j > 0:
            horz = D[i, j-1] + 1
        if diag <= vert and diag <= horz:
            xst.append(x[i-1])
            yst.append(y[j-1])
            i -= 1; j -= 1
        elif vert <= horz:
            xst.append(x[i-1])
            yst.append('-')
            i -= 1
        else:
            xst.append('-')
            yst.append(y[j-1])
            j -= 1
    return (''.join(xst))[::-1], (''.join(yst))[::-1]

def editDistance(x, y, stacked=False):
    ''' Calculate edit distance between sequences x and y using
        matrix dynamic programming.  Return distance. '''
    D = zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    if stacked:
        return D[len(x), len(y)], traceStacked(D, x, y)
    else:
        return D[len(x), len(y)], traceTranscript(D, x, y)
