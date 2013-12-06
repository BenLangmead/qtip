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

def editDistance(x, y):
    ''' Calculate edit distance between sequences x and y using
        matrix dynamic programming.  Return distance. '''
    D = zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D[len(x), len(y)]
