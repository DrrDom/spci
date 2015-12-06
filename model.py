#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2014
# license         : GPL3
#==============================================================================

import os
import argparse
import pprint
import numpy as np
import sklearn.cross_validation as cv

from time import strftime
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, RandomForestClassifier, \
    GradientBoostingClassifier
from sklearn.cross_decomposition import PLSRegression
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.externals import joblib
from sklearn.grid_search import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn.metrics import make_scorer
from multiprocessing import cpu_count
from scipy.sparse import dok_matrix


def load_sirms_txt(fname):
    with open(fname) as f:
        descr_names = f.readline().strip().split("\t")[1:]
        case_names = []
        x = []
        for line in f:
            tmp = line.strip().split("\t")
            case_names.append(tmp[0])
            x.append(tuple(map(float, tmp[1:])))
    return descr_names, case_names, np.asarray(x)


def load_sirms_svm(fname):
    descr_names = [v.strip() for v in open(os.path.splitext(fname)[0] + '.colnames').readlines()]
    case_names = [v.strip() for v in open(os.path.splitext(fname)[0] + '.rownames').readlines()]
    x = dok_matrix((len(case_names), len(descr_names)), dtype=np.float32)
    with open(fname) as f:
        for row, line in enumerate(f):
            tmp = line.strip().split(' ')
            for v in tmp:
                col, value = v.split(':')
                x[row, int(col)] = value
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
    joblib.dump(model, fname, compress=1, cache_size=1e8)


def mean(l):
    return sum(l) / len(l)


def get_sensitivity(obs, pred, classes=(0, 1)):
    tp = sum(obs + pred == 2 * classes[1])
    t = sum(obs == classes[1])
    return tp/t


def get_specificity(obs, pred, classes=(0, 1)):
    tn = sum(obs + pred == 2 * classes[0])
    n = sum(obs == classes[0])
    return tn/n


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
             "mse": get_mse(y, pred)}

    elif model_type == 'class':

        d = {"accuracy": (get_sensitivity(y, pred) + get_specificity(y, pred)) / 2,
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
            lines = ["Regression\n", "Time\tModel\tR2\tMSE\tOptimal_parameters\n"]
        elif model_type == 'class':
            lines = ["Classification\n", "Time\tModel\tAccuracy\tSensitivity\tSpecificity\tKappa\tOptimal_parameters\n"]

    if model_type == 'reg':
        lines.append(strftime("%Y-%m-%d %H:%M:%S") + "\t" + model_name + "\t" + "{0:.2f}".format(d["r2"]) + "\t" +
             "{0:.2f}".format(d["mse"]) + "\t" + model_params + "\n")

    elif model_type == 'class':
        lines.append(strftime("%Y-%m-%d %H:%M:%S") + "\t" + model_name + "\t" + "{0:.2f}".format(d["accuracy"]) + "\t" +
                     "{0:.2f}".format(d["sensitivity"]) + "\t" + "{0:.2f}".format(d["specificity"]) + "\t" +
                     "{0:.2f}".format(d["kappa"]) + "\t" +model_params + "\n")

    if verbose:
        print(lines[-1].strip())

    open(file_name, "wt").writelines(lines)


def main_params(x_fname, y_fname, model_names, ncores, model_type, verbose, cv_predictions, input_format):

    seed = 42

    # create models subdir
    model_dir = os.path.join(os.path.dirname(x_fname), "models")
    if not os.path.exists(model_dir):
        os.mkdir(model_dir)

    model_stat_fname = os.path.join(model_dir, "models_stat.txt")

    # load x and scale
    if input_format == 'txt':
        descr_names, mol_names, x = load_sirms_txt(x_fname)
    elif input_format == 'svm':
        descr_names, mol_names, x = load_sirms_svm(x_fname)

    scale = StandardScaler().fit(x)
    save_object(scale, os.path.join(model_dir, "scale.pkl"))
    x = scale.transform(x)

    # load y and sort y values according to x rows
    y = load_y(y_fname)
    y = np.asarray([y[n] for n in mol_names])

    cv5 = cv.KFold(n=len(y), n_folds=5, random_state=42, shuffle=True)

    cv_pred = np.copy(y)

    # build models

    for current_model in model_names:

        if verbose:
            print(current_model.upper() + ' model building...')

        if current_model == "rf":

            # choosing optimal parameters
            param_grid = {"max_features": [x.shape[1] // 10, x.shape[1] // 7, x.shape[1] // 5, x.shape[1] // 3],
                          "n_estimators": [500]}
            if model_type == "reg":
                m = GridSearchCV(RandomForestRegressor(), param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            elif model_type == "class":
                m = GridSearchCV(RandomForestClassifier(), param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            m.fit(x, y)

            # final model
            if model_type == "reg":
                m = RandomForestRegressor(n_estimators=m.best_params_["n_estimators"],
                                          max_features=m.best_params_["max_features"],
                                          bootstrap=True, random_state=seed)
            elif model_type == "class":
                m = RandomForestClassifier(n_estimators=m.best_params_["n_estimators"],
                                           max_features=m.best_params_["max_features"],
                                           bootstrap=True, random_state=seed)

        if current_model == "gbm":

            # choosing optimal parameters
            param_grid = {"n_estimators": [100, 200, 300, 400, 500]}
            if model_type == "reg":
                m = GridSearchCV(GradientBoostingRegressor(subsample=0.5, max_features=0.5),
                                 param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            elif model_type == "class":
                m = GridSearchCV(GradientBoostingClassifier(subsample=0.5, max_features=0.5),
                                 param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            m.fit(x, y)

            # final model
            if model_type == "reg":
                m = GradientBoostingRegressor(n_estimators=m.best_params_["n_estimators"],
                                              subsample=0.5, max_features=0.5, random_state=seed)
            elif model_type == "class":
                m = GradientBoostingClassifier(n_estimators=m.best_params_["n_estimators"],
                                               subsample=0.5, max_features=0.5, random_state=seed)

        if current_model == "svm":

            # choosing optimal parameters
            param_grid = {"C": [10 ** i for i in range(0, 5)],
                          "gamma": [10 ** i for i in range(-6, 0)]}
            if model_type == "reg":
                m = GridSearchCV(svm.SVR(kernel='rbf'), param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            elif model_type == "class":
                m = GridSearchCV(svm.SVC(kernel='rbf'), param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            m.fit(x, y)

            # final model
            if model_type == "reg":
                m = svm.SVR(kernel='rbf', C=m.best_params_["C"], gamma=m.best_params_["gamma"])
            elif model_type == "class":
                m = svm.SVC(kernel='rbf', C=m.best_params_["C"], gamma=m.best_params_["gamma"],
                            probability=True, random_state=seed)

        if current_model == "pls" and model_type == "reg":

            # choosing optimal parameters
            param_grid = {"n_components": [i for i in range(1, 8)]}
            m = GridSearchCV(PLSRegression(), param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            m.fit(x, y)

            # final model
            m = PLSRegression(n_components=m.best_params_["n_components"])

        if current_model == "knn":
            # choosing optimal parameters
            param_grid = {"n_neighbors": [i for i in range(3, 21)]}
            if model_type == "reg":
                m = GridSearchCV(KNeighborsRegressor(), param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            elif model_type == "class":
                m = GridSearchCV(KNeighborsClassifier(), param_grid, n_jobs=ncores, cv=cv5, refit=False, verbose=verbose)
            m.fit(x, y)

            # final model
            if model_type == "reg":
                m = KNeighborsRegressor(n_neighbors=m.best_params_["n_neighbors"])
            elif model_type == "class":
                m = KNeighborsClassifier(n_neighbors=m.best_params_["n_neighbors"])

        # return cv predictions
        pred = cv.cross_val_predict(m, x, y, cv5)
        if current_model == 'pls':
            pred = pred.reshape(pred.shape[0] * pred.shape[1])

        # build final model, save it and its stat
        m.fit(x, y)

        save_object(m, os.path.join(model_dir, current_model + '.pkl'))
        save_model_stat_2(current_model, model_stat_fname, str(m.get_params())[1:-1], y, pred,
                          model_type, verbose)

        # save cv predictions
        if cv_predictions:
            cv_pred = np.vstack((cv_pred, pred))

        if verbose:
            print(current_model.upper() + ' model was built\n')

    # save CV predictions
    if cv_predictions:
        np.savetxt(os.path.join(model_dir, "models_cv_pred.txt"),
                   np.column_stack([mol_names, np.round(np.transpose(cv_pred), 3)]),
                   fmt="%s",
                   delimiter="\t",
                   comments="",
                   header="Mol\tObs\t" + "\t".join(model_names))


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
    parser.add_argument('-c', '--ncores', metavar='1|2|3|...|all', default=str(cpu_count() - 1),
                        help='number of cores used for models building. Default: all cores - 1.')
    parser.add_argument('-t', '--model_type', metavar='reg|class', default='reg',
                        help='models type: reg for regression and class for classification.')
    parser.add_argument('-v', '--verbose', default=1,
                        help='Integer value. 0 - print no details. 1 and more - verbose output. Default: 1.')
    parser.add_argument('-p', '--cv_predictions', action='store_true', default=True,
                        help='True/False to output cross-validation predictions.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "descriptors": x_fname = v
        if o == "property": y_fname = v
        if o == "models": model_names = v
        if o == "ncores":
            if v == 'all':
                ncores = cpu_count()
            else:
                ncores = int(str)
                if ncores > cpu_count():
                    ncores = cpu_count()
        if o == "model_type": model_type = v
        if o == "verbose": verbose = int(v)
        if o == "cv_predictions": cv_predictions = v
        if o == "descriptors_format": input_format = v
        if input_format not in ['txt', 'svm']:
            print("Input file format is wrong - %s. Only txt and svm are allowed." % input_format)
            exit()

    main_params(x_fname, y_fname, model_names, ncores, model_type, verbose, cv_predictions, input_format)


if __name__ == '__main__':
    main()
