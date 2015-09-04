//
//  align.cpp
//  qsim
//
//  Created by Ben Langmead on 9/4/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#include "align.h"

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
