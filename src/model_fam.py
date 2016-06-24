import copy
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, GradientBoostingRegressor

__author__ = 'langmead'


class ModelFamily(object):
    """ Model family and simple hill-climbing hyperparameter search. """

    def __init__(self, name, new_predictor, params, round_to, min_separation, start_in_middle=True, calculates_oob=True):
        """
        name: name indicating the type of model, e.g. "RandomForestRegressor"
        new_predictor: function taking hyperparameters and returning new predictor with those params
        params: hyperparameters
        round_to: round scores off to this granularity when hill climbing
        min_separation: require an improvement of at least this much when optimizing
        start_in_middle: if true, hyperparameter search starts in middle of hyperparameter ranges.
                         if false, starts at index 0
        """
        self.name = name
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
        self._calculates_oob = calculates_oob

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

    def predictor_from_params(self, params):
        return self.new_predictor(params)

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

    def calculates_oob(self):
        return self._calculates_oob

    def best_predictor(self):
        return self.best_translated_params, self.new_predictor(self.best_translated_params)


def random_forest_models(random_seed=33, round_to=1e-5, min_separation=0.01):
    def _gen(params):
        return RandomForestRegressor(n_estimators=int(round(params[0])),
                                     max_depth=int(round(params[1])),
                                     random_state=random_seed,
                                     max_features=2,
                                     oob_score=True, bootstrap=True)
    return lambda: ModelFamily('RandomForestRegressor',
                               _gen, [range(5, 40, 2), range(4, 11, 1)],
                               round_to, min_separation=min_separation)


def extra_trees_models(random_seed=33, round_to=1e-5, min_separation=0.002):
    def _gen(params):
        return ExtraTreesRegressor(n_estimators=int(round(params[0])),
                                   max_depth=int(round(params[1])),
                                   random_state=random_seed,
                                   max_features=0.5,
                                   oob_score=True, bootstrap=True)
    return lambda: ModelFamily('ExtraTreesRegressor',
                               _gen, [range(5, 40, 2), range(4, 11, 1)],
                               round_to, min_separation=min_separation)


def gradient_boosting_models(random_seed=33, round_to=1e-5, min_separation=0.002):
    def _gen(params):
        return GradientBoostingRegressor(
            n_estimators=int(round(params[0])),
            max_depth=int(round(params[1])),
            learning_rate=params[2],
            random_state=random_seed,
            loss='ls')
    return lambda: ModelFamily('GradientBoostingRegressor',
                               _gen, [range(20, 60, 2), range(1, 11, 1), [x/10.0 for x in range(5, 11, 1)]],
                               round_to, min_separation=min_separation, calculates_oob=False)


def model_family(family, seed, optimization_tolerance):
    """ Given command-line arguments, return appropriate model family """
    if family == 'RandomForest':
        return random_forest_models(seed, optimization_tolerance)
    elif family == 'ExtraTrees':
        return extra_trees_models(seed, optimization_tolerance)
    elif family == 'GradientBoostingRegressor':
        return gradient_boosting_models(seed, optimization_tolerance)
    else:
        raise RuntimeError('Bad value for --model-family: "%s"' % family)
