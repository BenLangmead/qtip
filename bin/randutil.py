import random
import bisect


class WeightedRandomGenerator(object):
    """ Given a dictionary associating keys with weights, with each
        call to next() return a key with probability equal to that
        key's fraction of the total weight. """
    
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
    """ Simple reservoir sampler """
    
    def __init__(self, k):
        """ Initialize given k, the size of the reservoir """
        self.k = k
        self.r = []
        self.n = 0

    def add_step_1(self):
        """ Step 1 of stateful add: check whether this will actually get added """
        if self.n < self.k:
            return self.n
        else:
            j = random.randint(0, self.n)
            return j if j < self.k else None

    def add_step_2(self, j, obj):
        """ Step 2 of stateful add: actually add it """
        if j >= len(self.r):
            assert j == len(self.r)
            self.r.append(obj)
        else:
            self.r[j] = obj
        self.n += 1
    
    def add(self, obj):
        """ Add object to sampling domain """
        if self.n < self.k:
            self.r.append(obj)
        else:
            j = random.randint(0, self.n)
            if j < self.k:
                self.r[j] = obj
        self.n += 1
    
    def draw(self):
        """ Make a draw from the reservoir """
        return random.choice(self.r)
    
    def empty(self):
        """ Return true iff no items have been added """
        return len(self) == 0

    def __iter__(self):
        """ Return iterator over the sample """
        return iter(self.r)

    def num_added(self):
        return self.n

    def num_sampled(self):
        return len(self.r)
