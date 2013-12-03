import random
import bisect

class WeightedRandomGenerator(object):
    ''' Given a dictionary associating keys with weights, with each
        call to next() return a key with probability equal to that
        key's fraction of the total weight. '''
    
    def __init__(self, weights):
        self.keys = sorted(weights.keys())
        self.totals = []
        cum = 0
        for k in self.keys:
            cum += weights[k]
            self.totals.append(cum)
    
    def next(self):
        rnd = random.random() * self.totals[-1]
        return self.keys[bisect.bisect_right(self.totals, rnd)]

class ReservoirSampler(object):
    ''' Simple reservoir sampler '''
    
    def __init__(self, k):
        ''' Initialize given k, the size of the reservoir '''
        self.k = k
        self.r = []
        self.n = 0
    
    def add(self, obj):
        ''' Add object to sampling domain '''
        if self.n < self.k:
            self.r.append(obj)
        else:
            j = random.randint(0, self.n)
            if j < self.k:
                self.r[j] = obj
        self.n += 1
    
    def draw(self):
        ''' Make a draw from the reservoir '''
        return random.choice(self.r)
    
    def __len__(self):
        ''' Return number of items added to the domain '''
        return self.n
    
    def empty(self):
        ''' Return true iff no items have been added '''
        return len(self) == 0
