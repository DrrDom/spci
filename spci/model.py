#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2014
# license         : LGPLv3
#==============================================================================

import os
import argparse
import numpy as np
import sklearn.model_selection as ms
import warnings

from time import strftime
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, RandomForestClassifier, \
    GradientBoostingClassifier
from sklearn.cross_decomposition import PLSRegression
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
import joblib
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from multiprocessing import cpu_count
from scipy.sparse import dok_matrix

from .predict import get_major_vote


def load_sirms_txt(fname, names=None):
    if names:
        names = set(names)
    with open(fname) as f:
        descr_names = f.readline().strip().split("\t")[1:]
        case_names = []
        x = []
        for line in f:
            tmp = line.strip().split("\t")
            if not names or tmp[0] in names:
                case_names.append(tmp[0])
                x.append(tuple(map(float, tmp[1:])))

    return descr_names, case_names, np.asarray(x)


def load_sirms_svm(fname, names=None):
    descr_names = [v.strip() for v in open(os.path.splitext(fname)[0] + '.colnames').readlines()]
    case_names = [v.strip() for v in open(os.path.splitext(fname)[0] + '.rownames').readlines()]
    if names:
        case_names_ids = set(i for i, n in enumerate(case_names) if n in names)
        case_names = [case_names[i] for i in sorted(case_names_ids)]
    else:
        case_names_ids = set(range(len(case_names)))
    x = dok_matrix((len(case_names_ids), len(descr_names)), dtype=np.float32)
    with open(fname) as f:
        row = 0
        for i, line in enumerate(f):
            if i in case_names_ids:
                tmp = line.strip().split(' ')
                for v in tmp:
                    col, value = v.split(':')
                    x[row, int(col)] = value
                row += 1
    return descr_names, case_names, x.toarray()


def load_y(fname):
    with open(fname) as f:
        d = dict()
        f.readline()
        for line in f:
            tmp = line.strip().split('\t')
            d[tmp[0]] = float(tmp[1])
    return d


def save_object(model, fname):
    joblib.dump(model, fname)


def mean(l):
    return sum(l) / len(l)


def get_sensitivity(obs, pred, classes=(0, 1)):
    tp = sum(np.logical_and((pred == classes[1]), (obs == pred)))
    # tp = sum(obs + pred == 2 * classes[1])
    t = sum(obs == classes[1])
    return tp / t


def get_specificity(obs, pred, classes=(0, 1)):
    tn = sum(np.logical_and((pred == classes[0]), (obs == pred)))
    # tn = sum(obs + pred == 2 * classes[0])
    n = sum(obs == classes[0])
    return tn / n


def get_kappa(obs, pred, classes=(0, 1)):
    obs_0 = sum(obs == classes[0])
    pred_0 = sum(pred == classes[0])
    obs_1 = sum(obs == classes[1])
    pred_1 = sum(pred == classes[1])
    baseline = (obs_0 * pred_0 + obs_1 * pred_1) / (len(obs) ** 2)
    acc = sum(obs == pred) / len(obs)
    return (acc - baseline) / (1 - baseline)


def get_r2_test(obs, pred, mean_training):
    return 1 - (sum((obs - pred) ** 2) / sum((obs - mean_training) ** 2))


def get_mse(obs, pred):
    return sum((obs - pred) ** 2) / (len(pred) - 1)


def calc_model_stat_2(y, pred, model_type):

    if model_type == 'reg':

        d = {"r2": get_r2_test(y, pred, mean(y)),
             "mse": get_mse(y, pred),
             "rmse": get_mse(y, pred) ** 0.5}

    elif model_type == 'class':

        d = {"balanced accuracy": (get_sensitivity(y, pred) + get_specificity(y, pred)) / 2,
             "sensitivity": get_sensitivity(y, pred),
             "specificity": get_specificity(y, pred),
             "kappa": get_kappa(y, pred)}

    return d


def save_model_stat_2(model_name, file_name, model_params, y, pred, model_type, verbose):

    d = calc_model_stat_2(y, pred, model_type)

    if os.path.isfile(file_name):
        lines = open(file_name).readlines()
        lines = lines[0:2] + [line for line in lines[2:] if line.strip().split('\t')[1] != model_name]

    else:
        if model_type == 'reg':
            lines = ["Regression\n", "Time\tModel\tR2\tRMSE\tMSE\tOptimal_parameters\n"]
        elif model_type == 'class':
            lines = ["Classification\n", "Time\tModel\tBalanced accuracy\tSensitivity\tSpecificity\tKappa\tOptimal_parameters\n"]

    if model_type == 'reg':
        lines.append(strftime("%Y-%m-%d %H:%M:%S") + "\t" + model_name + "\t" + "{0:.2f}".format(d["r2"]) + "\t" +
             "{0:.2f}".format(d["rmse"]) + "\t" + "{0:.2f}".format(d["mse"]) + "\t" + model_params + "\n")

    elif model_type == 'class':
        lines.append(strftime("%Y-%m-%d %H:%M:%S") + "\t" + model_name + "\t" + "{0:.2f}".format(d["balanced accuracy"]) + "\t" +
                     "{0:.2f}".format(d["sensitivity"]) + "\t" + "{0:.2f}".format(d["specificity"]) + "\t" +
                     "{0:.2f}".format(d["kappa"]) + "\t" + model_params + "\n")

    if verbose:
        print(lines[-1].strip())

    open(file_name, "wt").writelines(lines)


def save_bound_box_constrains(x, fname):
    d = dict()
    d['min_values'] = tuple(np.amin(x, axis=0))
    d['max_values'] = tuple(np.amax(x, axis=0))
    save_object(d, fname)


def add_obj_to_file(fname, obj):
    with open(fname, "ab") as f:
        joblib.dump(obj, f)


def make_subsets(y, seed=42):

    np.random.seed(seed)
    ids_min = []
    subsets = []
    if len(np.where(y == 0)[0]) / len(np.where(y == 1)[0]) > 1.5:
        ids_min = list(np.where(y == 1)[0])
        ids_maj = list(np.where(y == 0)[0])
    elif len(np.where(y == 1)[0]) / len(np.where(y == 0)[0]) > 1.5:
        ids_min = list(np.where(y == 0)[0])
        ids_maj = list(np.where(y == 1)[0])

    if ids_min:
        ratio = round(len(ids_maj) / len(ids_min))
        np.random.shuffle(ids_maj)
        np.random.shuffle(ids_min)
        subset_sizes = (len(ids_maj) // ratio) * np.ones(ratio, dtype=np.int)
        subset_sizes[:len(ids_maj) % ratio] += 1
        current = 0
        for subset_size in subset_sizes:
            start, stop = current, current + subset_size
            subsets.append(ids_maj[start:stop] + ids_min)
            current = stop

    return subsets


def main_params(x_fname, y_fname, model_names, models_dir, ncores, model_type, verbose, cv_predictions, imbalanced, input_format):

    seed = 42

    # create models subdir
    if models_dir is None:
        models_dir = os.path.join(os.path.dirname(x_fname), "models")
    if not os.path.exists(models_dir):
        os.makedirs(models_dir)
    for m in model_names:
        fpath = os.path.join(models_dir, m + ".pkl")
        if os.path.exists(fpath):
            os.remove(fpath)

    model_stat_fname = os.path.join(models_dir, "models_stat.txt")

    # load y
    y = load_y(y_fname)

    # load x
    if input_format == 'txt':
        descr_names, mol_names, x = load_sirms_txt(x_fname, names=y.keys())
    elif input_format == 'svm':
        descr_names, mol_names, x = load_sirms_svm(x_fname, names=y.keys())
    else:
        print("Illegal value of input format: " % input_format)
        exit()

    # process y
    y = np.asarray([y[n] for n in mol_names])

    # process x
    save_bound_box_constrains(x, os.path.join(models_dir, "bound_box.pkl"))
    save_object(descr_names, os.path.join(models_dir, "var_names.pkl"))

    # scale
    scale = StandardScaler().fit(x)
    save_object(scale, os.path.join(models_dir, "scale.pkl"))

    x = scale.transform(x)

    if model_type == "class":
        cv = ms.StratifiedKFold(n_splits=5, random_state=seed, shuffle=True)
    elif model_type == "reg":
        cv = ms.KFold(n_splits=5, random_state=seed, shuffle=True)

    if model_type == "class" and imbalanced:
        subsets = make_subsets(y, seed)
        if not subsets:
            warnings.warn("The data set is balanced (ratio majority:minority < 1.5)."
                          "No multiple undersampling will be done", Warning)
            subsets = [list(range(y.shape[0]))]
    else:
        subsets = [list(range(y.shape[0]))]

    # build models

    for current_model in model_names:

        if verbose:
            print(current_model.upper() + ' model building...')

        models_lst = []  # this lst refreshes on each model name; here we store either 1 model in balanced case or list of models=number of subsets in the case of imbalanced

        for subset in subsets:

            if current_model == "rf":

                # choosing optimal parameters
                param_grid = {"max_features": [x.shape[1] // 10, x.shape[1] // 7, x.shape[1] // 5, x.shape[1] // 3],
                              "n_estimators": [500]}

                if model_type == "reg":
                    m = ms.GridSearchCV(RandomForestRegressor(random_state=seed), param_grid, n_jobs=ncores, cv=cv, refit=False, verbose=verbose)
                elif model_type == "class":
                    m = ms.GridSearchCV(RandomForestClassifier(random_state=seed), param_grid, n_jobs=ncores, cv=cv, refit=False, verbose=verbose)
                m.fit(x[subset], y[subset])

                # final model
                if model_type == "reg":
                    m = RandomForestRegressor(n_estimators=m.best_params_["n_estimators"],
                                              max_features=m.best_params_["max_features"],
                                              bootstrap=True, random_state=seed)
                elif model_type == "class":
                    m = RandomForestClassifier(n_estimators=m.best_params_["n_estimators"],
                                               max_features=m.best_params_["max_features"],
                                               bootstrap=True, random_state=seed)
                models_lst.append(m)

            if current_model == "gbm":

                # choosing optimal parameters
                param_grid = {"n_estimators": [100, 200, 300, 400, 500]}

                if model_type == "reg":
                    m = ms.GridSearchCV(GradientBoostingRegressor(subsample=0.5, max_features=0.5, random_state=seed),
                                        param_grid, n_jobs=ncores, cv=cv, refit=False, verbose=verbose)
                elif model_type == "class":
                    m = ms.GridSearchCV(GradientBoostingClassifier(subsample=0.5, max_features=0.5, random_state=seed),
                                        param_grid, n_jobs=ncores, cv=cv, refit=False, verbose=verbose)
                m.fit(x[subset], y[subset])

                # final model
                if model_type == "reg":
                    m = GradientBoostingRegressor(n_estimators=m.best_params_["n_estimators"],
                                                  subsample=0.5, max_features=0.5, random_state=seed)
                elif model_type == "class":
                    m = GradientBoostingClassifier(n_estimators=m.best_params_["n_estimators"],
                                                   subsample=0.5, max_features=0.5, random_state=seed)
                models_lst.append(m)

            if current_model == "svm":

                # choosing optimal parameters
                if model_type == "reg":
                    param_grid = {"C": [10 ** i for i in range(0, 5)],
                                  "epsilon": [0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01]}
                    m = ms.GridSearchCV(svm.SVR(kernel='rbf'), param_grid, n_jobs=ncores, cv=cv, refit=False,
                                        verbose=verbose)
                elif model_type == "class":
                    param_grid = {"C": [10 ** i for i in range(0, 5)],
                                  "gamma": [10 ** i for i in range(-6, 0)]}
                    m = ms.GridSearchCV(svm.SVC(kernel='rbf', random_state=seed), param_grid, n_jobs=ncores, cv=cv, refit=False,
                                        verbose=verbose)
                m.fit(x[subset], y[subset])

                # final model
                if model_type == "reg":
                    m = svm.SVR(kernel='rbf', C=m.best_params_["C"], epsilon=m.best_params_["epsilon"])
                elif model_type == "class":
                    m = svm.SVC(kernel='rbf', C=m.best_params_["C"], gamma=m.best_params_["gamma"],
                                probability=True, random_state=seed)
                models_lst.append(m)

            if current_model == "pls" and model_type == "reg":

                # choosing optimal parameters
                param_grid = {"n_components": [i for i in range(1, 8)]}
                m = ms.GridSearchCV(PLSRegression(), param_grid, n_jobs=ncores, cv=cv, refit=False, verbose=verbose)
                m.fit(x[subset], y[subset])

                # final model
                m = PLSRegression(n_components=m.best_params_["n_components"])
                models_lst.append(m)

            if current_model == "knn":

                # choosing optimal parameters
                param_grid = {"n_neighbors": [i for i in range(3, 21)]}
                if model_type == "reg":
                    m = ms.GridSearchCV(KNeighborsRegressor(), param_grid, n_jobs=ncores, cv=cv, refit=False, verbose=verbose)
                elif model_type == "class":
                    m = ms.GridSearchCV(KNeighborsClassifier(), param_grid, n_jobs=ncores, cv=cv, refit=False, verbose=verbose)
                m.fit(x[subset], y[subset])

                # final model
                if model_type == "reg":
                    m = KNeighborsRegressor(n_neighbors=m.best_params_["n_neighbors"])
                elif model_type == "class":
                    m = KNeighborsClassifier(n_neighbors=m.best_params_["n_neighbors"])
                models_lst.append(m)

        # return cv predictions
        ncol = len(models_lst) + 1 if len(models_lst) > 1 else len(models_lst)   # +1 column if consensus
        cv_pred = np.column_stack((y, np.full((y.shape[0], ncol), np.nan)))

        for i, (m, subset) in enumerate(zip(models_lst, subsets)):
            pred = ms.cross_val_predict(estimator=m, X=x[subset], y=y[subset], cv=cv)
            if current_model == 'pls':   # reshape for pls because it returns 2d array and we need 1d
                pred = pred.reshape(len(subset))
            cv_pred[subset, i + 1] = pred

            # build final model, save it and its stat
            m.fit(x[subset], y[subset])
            add_obj_to_file(os.path.join(models_dir, current_model + '.pkl'), m)
            save_model_stat_2(current_model + '_%i' % i, model_stat_fname, str(m.get_params())[1:-1],
                              y[subset],
                              cv_pred[subset, i + 1],
                              model_type,
                              verbose)

        # calc cv consensus and save stat
        if model_type == "class" and len(models_lst) > 1:
            cv_pred[:, -1] = np.apply_along_axis(get_major_vote, 1, cv_pred[:, 1:])
            # cv_pred[:, -1] = np.around(np.nanmean(cv_pred[:, 1:], axis=1))
            save_model_stat_2(current_model + "_consensus", model_stat_fname, "",
                              y,
                              cv_pred[:, -1],
                              model_type,
                              verbose)

        # save cv predictions
        if cv_predictions:
            np.savetxt(os.path.join(models_dir, current_model + "_cv_pred.txt"),
                       np.column_stack([mol_names, np.round(cv_pred, 3)]),
                       fmt="%s",
                       delimiter="\t",
                       comments="",
                       header="Mol\tObs\t" +
                              "\t".join("%s_%i" % (current_model, i) for i in range(len(models_lst))) +
                              "\t" + current_model + "_consensus")

        if verbose:
            print(current_model.upper() + ' model was built\n')


def main():

    parser = argparse.ArgumentParser(description='Build QSAR/QSPR models.')
    parser.add_argument('-x', '--descriptors', metavar='descriptors.txt', required=True,
                        help='input file with descriptors in txt or svm formats.')
    parser.add_argument('-f', '--descriptors_format', metavar='txt', default='txt',
                        help='format of the input file with descriptors (txt|svm). Svm file has to have two '
                             'supplementary files with the same name and extensions rownames and colnames. '
                             'Default txt.')
    parser.add_argument('-y', '--property', metavar='property.txt', required=True,
                        help='input text file with property/activity values.')
    parser.add_argument('-m', '--models', metavar='[rf gbm svm pls knn]', default=['rf'], nargs='*',
                        help='models to build.')
    parser.add_argument('-d', '--models_dir', default=None, help='path to dir where models will be saved.')
    parser.add_argument('-c', '--ncores', metavar='1|2|3|...|all', default=str(cpu_count() - 1),
                        help='number of cores used for models building. Default: all cores - 1.')
    parser.add_argument('-t', '--model_type', metavar='reg|class', default='reg',
                        help='models type: reg for regression and class for classification.')
    parser.add_argument('-v', '--verbose', default=1,
                        help='Integer value. 0 - print no details. 1 and more - verbose output. Default: 1.')
    parser.add_argument('--imbalanced', default=False, action='store_true',
                        help='if set this flag the input data set will be considered imbalanced and '
                             'multiple undersampling will be applied. Works only with classification tasks '
                             'when the ratio between two classes is >= 1.5.')
    parser.add_argument('-p', '--cv_predictions', action='store_true', default=True,
                        help='True/False to save cross-validation predictions in text files along '
                             'with model files.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "descriptors": x_fname = v
        if o == "property": y_fname = v
        if o == "models": model_names = v
        if o == "ncores":
            if v == 'all':
                ncores = cpu_count()
            else:
                ncores = int(v)
                if ncores > cpu_count():
                    ncores = cpu_count()
        if o == "model_type": model_type = v
        if o == "verbose": verbose = int(v)
        if o == "cv_predictions": cv_predictions = v
        if o == "descriptors_format": input_format = v
        if o == "models_dir": models_dir = v
        if o == "imbalanced": imbalanced = v
    if input_format not in ['txt', 'svm']:
        print("Input file format is wrong - %s. Only txt and svm are allowed." % input_format)
        exit()

    main_params(x_fname=x_fname,
                y_fname=y_fname,
                model_names=model_names,
                models_dir=models_dir,
                ncores=ncores,
                model_type=model_type,
                verbose=verbose,
                cv_predictions=cv_predictions,
                imbalanced=imbalanced,
                input_format=input_format)


if __name__ == '__main__':
    main()
