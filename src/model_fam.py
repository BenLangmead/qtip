import copy
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor

__author__ = 'langmead'


class ModelFamily(object):
    """ Encapsulates a model family and a simple interface for search
        hyperparameter space. """

    def __init__(self, new_predictor, params, round_to, min_separation, start_in_middle=True):
        """
        new_predictor: function that takes set of params and returns
                       new predictor with those params
        params: list of lists
        """
        self.new_predictor = new_predictor  # returns new predictor given parameters
        self.params = params  # space of possible parameter choices
        self.last_params = None  # remember last set of params used for predictor
        self.round_to = round_to
        self.min_separation = min_separation  # have to improve by at least this much to declare new best
        self.best = self.best_base = self.best_rounded = float('-inf')
        self.best_translated_params = None
        if start_in_middle:
            center = tuple([int(round(len(x) / 2)) for x in params])
        else:
            center = tuple([0] * len(params))
        self.workset = {center}
        self.added_to_workset = copy.copy(self.workset)
        self._add_neighbors_to_workset(center)

    def _add_neighbors_to_workset(self, center):
        for i in range(len(self.params)):
            if center[i] > 0:
                neighbor = list(center[:])
                neighbor[i] -= 1  # add next-lowest neighbor
                neighbor = tuple(neighbor)
                if neighbor not in self.added_to_workset:
                    self.added_to_workset.add(neighbor)
                    self.workset.add(neighbor)
            if center[i] < len(self.params[i])-1:
                neighbor = list(center[:])
                neighbor[i] += 1  # add next-highest neighbor
                neighbor = tuple(neighbor)
                if neighbor not in self.added_to_workset:
                    self.added_to_workset.add(neighbor)
                    self.workset.add(neighbor)

    def _idxs_to_params(self, idxs):
        return [self.params[i][j] for i, j in enumerate(idxs)]

    def next_predictor(self):
        if len(self.workset) > 0:
            self.last_params = self.workset.pop()
            translated_params = self._idxs_to_params(self.last_params)
            return translated_params, self.new_predictor(translated_params)
        self.last_params = None
        return None, None

    def set_score(self, score):
        assert self.last_params is not None
        assert self.last_params in self.added_to_workset
        score_rounded = int(score / self.round_to) * self.round_to
        if self.best == float('-inf') or score > self.best_base * (1.0 + self.min_separation):
            self.best, self.best_base, self.best_rounded = score, score, score_rounded
            self.best_translated_params = self._idxs_to_params(self.last_params)
            self._add_neighbors_to_workset(self.last_params)
            return True, True
        elif score > self.best:
            self.best, self.best_rounded = score, score_rounded
            self.best_translated_params = self._idxs_to_params(self.last_params)
            return True, False
        return False, False

    def best_predictor(self):
        return self.best_translated_params, self.new_predictor(self.best_translated_params)


def random_forest_models(random_seed=33, round_to=1e-5, min_separation=0.01):
    # These perform OK but not as well as the extremely random trees
    def _gen(params):
        return RandomForestRegressor(n_estimators=params[0], max_depth=params[1],
                                     random_state=random_seed,
                                     max_features=2,
                                     oob_score=True, bootstrap=True)
    return lambda: ModelFamily(_gen, [range(5, 105, 10), range(2, 10)],
                               round_to, min_separation=min_separation)


def extra_trees_models(random_seed=33, round_to=1e-5, min_separation=0.002):
    # These perform quite well
    def _gen(params):
        return ExtraTreesRegressor(n_estimators=params[0], max_depth=params[1],
                                   random_state=random_seed,
                                   max_features=0.5,
                                   oob_score=True, bootstrap=True)
    return lambda: ModelFamily(_gen, [range(5, 85, 2), range(3, 16, 1)],
                               round_to, min_separation=min_separation)


def model_family(family, seed, optimization_tolerance):
    """ Given command-line arguments, return appropriate model family """
    if family == 'RandomForest':
        return random_forest_models(seed, optimization_tolerance)
    elif family == 'ExtraTrees':
        return extra_trees_models(seed, optimization_tolerance)
    else:
        raise RuntimeError('Bad value for --model-family: "%s"' % family)
