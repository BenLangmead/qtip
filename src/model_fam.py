import copy
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, GradientBoostingRegressor

__author__ = 'langmead'


class ModelFamily(object):
    """ Model family and simple hill-climbing hyperparameter search. """

    def __init__(self, name, new_predictor, params, min_separation, start_in_middle=True, calculates_oob=True):
        """
        name: name indicating the type of model, e.g. "RandomForestRegressor"
        new_predictor: function taking hyperparameters and returning new predictor with those params
        params: hyperparameters
        min_separation: require an improvement of at least this much when optimizing
        start_in_middle: if true, hyperparameter search starts in middle of hyperparameter ranges.
                         if false, starts at index 0
        """
        self.name = name
        self.new_predictor = new_predictor  # returns new predictor given parameters
        self.params = params  # space of possible parameter choices
        self.last_params = None  # remember last set of params used for predictor
        self.min_separation = min_separation  # have to improve by at least this much to declare new best
        self.best = self.best_base = float('-inf')
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
        if self.best == float('-inf') or score > self.best_base * (1.0 + self.min_separation):
            self.best, self.best_base = score, score
            self.best_translated_params = self._idxs_to_params(self.last_params)
            self._add_neighbors_to_workset(self.last_params)
            return True, True
        elif score > self.best:
            self.best = score
            self.best_translated_params = self._idxs_to_params(self.last_params)
            return True, False
        return False, False

    def calculates_oob(self):
        return self._calculates_oob

    def best_predictor(self):
        return self.best_translated_params, self.new_predictor(self.best_translated_params)


def random_forest_models(random_seed, min_separation, num_trees_str, max_tree_depth_str, max_features_str):
    def _gen(params):
        return RandomForestRegressor(n_estimators=int(round(params[0])),
                                     max_depth=int(round(params[1])),
                                     random_state=random_seed,
                                     max_features=params[2],
                                     oob_score=True,
                                     bootstrap=True)
    num_trees = map(float, num_trees_str.split(','))
    max_tree_depth = map(float, max_tree_depth_str.split(','))
    max_features = map(float, max_features_str.split(','))
    return lambda: ModelFamily('RandomForestRegressor',
                               _gen, [num_trees, max_tree_depth, max_features],
                               min_separation=min_separation)


def extra_trees_models(random_seed, min_separation, num_trees_str, max_tree_depth_str, max_features_str):
    def _gen(params):
        return ExtraTreesRegressor(n_estimators=int(round(params[0])),
                                   max_depth=int(round(params[1])),
                                   random_state=random_seed,
                                   max_features=params[2],
                                   oob_score=True,
                                   bootstrap=True)
    num_trees = map(float, num_trees_str.split(','))
    max_tree_depth = map(float, max_tree_depth_str.split(','))
    max_features = map(float, max_features_str.split(','))
    return lambda: ModelFamily('ExtraTreesRegressor',
                               _gen, [num_trees, max_tree_depth, max_features],
                               min_separation=min_separation)


def gradient_boosting_models(random_seed, min_separation,
                             num_trees_str, max_tree_depth_str, learning_rate_str):
    def _gen(params):
        return GradientBoostingRegressor(
            n_estimators=int(round(params[0])),
            max_depth=int(round(params[1])),
            learning_rate=params[2],
            random_state=random_seed,
            loss='ls')
    num_trees = map(float, num_trees_str.split(','))
    max_tree_depth = map(float, max_tree_depth_str.split(','))
    learning_rate = map(float, learning_rate_str.split(','))
    return lambda: ModelFamily('GradientBoostingRegressor',
                               _gen, [num_trees, max_tree_depth, learning_rate],
                               min_separation=min_separation, calculates_oob=False)


def add_args(parser):
    parser.add_argument('--model-family', metavar='family', type=str, required=False,
                        default='RandomForest', help='{RandomForest | ExtraTrees | GradientBoosting}')
    parser.add_argument('--model-params', metavar='parameters', type=str, required=False,
                        help='hyperparameters for model; for RandomForest or ExtraTrees: (#trees):(max depth)')

    default_num_trees_range = range(5, 60, 2)
    parser.add_argument('--num-trees', metavar='int,int,...', type=str,
                        default=','.join(map(str, default_num_trees_range)),
                        help='number of decision trees to try; relevant for all model families')

    default_max_depth_range = range(3, 15, 1)
    parser.add_argument('--max-tree-depth', metavar='int,int,...', type=str,
                        default=','.join(map(str, default_max_depth_range)),
                        help='maximum decision tree depth to try; relevant for all model families')

    default_max_features_range = [0.05, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    parser.add_argument('--max-features', metavar='float,float,...', type=str,
                        default=','.join(map(str, default_max_features_range)),
                        help='maximum number of features to consider at each decision tree node; '
                             'relevant for RandomForest and ExtraTrees')

    default_learning_rate_range = [0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    parser.add_argument('--learning-rate', metavar='float,float,...', type=str,
                        default=','.join(map(str, default_learning_rate_range)),
                        help='learning rate to use when fitting; only relevant for GradientBoosting')

    parser.add_argument('--optimization-tolerance', metavar='fraction', type=float, default=0.001,
                        help='When using hill climbing procedure to optimize hyperparamters,'
                             'stop when OOB score can\'t be improved by this relative factor')


def model_family(args, random_seed):
    """ Given command-line arguments, return appropriate model family """
    if args['model_family'] == 'RandomForest':
        if args['model_params'] is not None and args['model_params'].count(':') != 2:
            raise RuntimeError('--model-params for RandomForest must have 3 fields separated by :')
        return random_forest_models(random_seed, args['optimization_tolerance'],
                                    args['num_trees'], args['max_tree_depth'], args['max_features'])
    elif args['model_family'] == 'ExtraTrees':
        if args['model_params'] is not None and args['model_params'].count(':') != 2:
            raise RuntimeError('--model-params for ExtraTrees must have 3 fields separated by :')
        return extra_trees_models(random_seed, args['optimization_tolerance'],
                                  args['num_trees'], args['max_tree_depth'], args['max_features'])
    elif args['model_family'] == 'GradientBoosting':
        if args['model_params'] is not None and args['model_params'].count(':') != 2:
            raise RuntimeError('--model-params for GradientBoosting must have 3 fields separated by :')
        return gradient_boosting_models(random_seed, args['optimization_tolerance'],
                                        args['num_trees'], args['max_tree_depth'], args['learning_rate'])
    else:
        raise RuntimeError('Bad value for --model-family: "%s"' % args['model_family'])
