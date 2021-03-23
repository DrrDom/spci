#!/usr/bin/env python
# author          : Pavel
# date            : 05.11.15
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel 2015
# license         : LGPLv3
#==============================================================================

import os
import argparse
import numpy as np
import joblib
from collections import defaultdict, Counter
import re
import math
import operator

from .sirms_file_class import SirmsFile


def get_var_names(fname, fformat):
    # fformat = txt or svm
    output = None
    if fformat == 'txt':
        output = open(fname).readline().strip().split('\t')[1:]
    elif fformat == 'svm':
        with open(os.path.splitext(fname)[0] + '.colnames') as f:
            output = []
            for line in f:
                output.append(line.strip())
    return output


def adjust_dataset(x, x_col_names, ref_col_names, default_value=0):

    # init output with default value
    if default_value != 0:
        output = np.empty((np.shape(x)[0], len(ref_col_names)))
        output.fill(default_value)
    else:
        output = np.zeros((np.shape(x)[0], len(ref_col_names)))

    # get indices in output (-1 means no such column in output)
    ids = [-1] * len(x_col_names)
    for i, col in enumerate(x_col_names):
        try:
            ids[i] = ref_col_names.index(col)
        except ValueError:
            pass

    keep_cols = [i != -1 for i in ids]
    ids = [i for i in ids if i != -1]

    output[::, np.array(ids)] = x[::, np.array(keep_cols)]

    return output


def predict(x, model, model_name, model_type):

    pred = None

    if model_type == 'reg':
        pred = model.predict(x)
        if re.match("pls^|pls_", model_name):  # we added _0,1,2..
            pred = pred.reshape(pred.shape[0] * pred.shape[1])
        pred = pred.tolist()

    elif model_type == 'class':
        pos_id = model.classes_.tolist().index(1)
        pred = [i[pos_id] for i in model.predict_proba(x)]

    return pred


def add_consensus(pred_array, pred_model_names, ids=None, current_model_name=None):
    # ids and current_model_name are used when we compute not overall consensus
    # but consensus of multiple models of single ML
    # method (svm_0. svm_1...)(example: ids=[0,1,2], current_model_name="svm")

    if np.shape(pred_array)[1] > 1:  # otherwise nothing to compute
        if not ids:  # all columns will be selected
            ids = range(np.shape(pred_array)[1])
        pred_array = np.append(pred_array, np.mean(pred_array[:, ids], 1).reshape((-1, 1)), 1)
        if current_model_name:  # case when we compute consensus for single ml model
            pred_model_names.append(current_model_name + "_" + "consensus")
        else:  # case when we compute overall consensus
            pred_model_names.append("consensus")
    return pred_array  # this arrary returned because it is not changed outside function.  pred_model_names ont returned they are changed


def save_predictions(fname, pred_array, pred_model_names, mol_names, ad_dict):
    # add mol_names ass the first column
    pred_array = np.column_stack((mol_names, pred_array))
    pred_model_names.insert(0, 'Compounds')
    # add AD values as additional columns
    if ad_dict is not None:
        for k, v in ad_dict.items():
            pred_array = np.column_stack((pred_array, v))
            pred_model_names.append(k)
    np.savetxt(fname, pred_array, '%s', delimiter='\t', header='\t'.join(pred_model_names))


def load_bound_box_constrains(fname):
    d = joblib.load(fname)
    return d['min_values'], d['max_values']


def check_bound_box(x, bound_box_constrains):
    # bound_box_constrains: tuple of tuples of min and max values
    a = np.greater_equal(x, bound_box_constrains[0])
    b = np.less_equal(x, bound_box_constrains[1])
    c = np.logical_and(a, b)
    d = np.all(c, axis=1)
    return d


def load_models(model_dir, model_names):
    """
    Read all model objects from each model file and store in the list.
    Returns corresponding names of models by adding _N (e.g. svm_0, etc)
    :param model_dir: path to dir where models are stored
    :param model_names: list of model names (file names without .pkl extension)
    :return: list of models and list of model names
    """
    models = []
    pred_model_names = []
    for m in model_names:
        for i, j in enumerate(load_multi_obj(os.path.join(model_dir, m + ".pkl"))):
            models.append(j)
            pred_model_names.append(m + "_" + str(i))
    return models, pred_model_names


def load_multi_obj(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield joblib.load(f)
            except EOFError:
                break


def get_major_vote(data):
    """
    :param data: iterable
    :return: item which is most often occur
    """
    data = [i for i in data if not math.isnan(i)]
    return max(Counter(data).items(), key=operator.itemgetter(1))[0]


def main_params(x_fname, input_format, out_fname, model_names, model_dir, model_type, ad, verbose):

    if ad is not None:
        ad_dict = {item: [] for item in ad}
        bound_box_constrains = load_bound_box_constrains(os.path.join(model_dir, "bound_box.pkl"))
    else:
        ad_dict = None

    # scale
    if os.path.isfile(os.path.join(model_dir, "scale.pkl")):
        scale = joblib.load(os.path.join(model_dir, "scale.pkl"))
    else:
        scale = None

    pred = defaultdict(list)

    # load models
    models, pred_model_names = load_models(model_dir, model_names)

    input_sirms_file = SirmsFile(x_fname, file_format=input_format, chunks=1000)
    ref_var_names = joblib.load(os.path.join(model_dir, "var_names.pkl"))

    input_sirms_file.reset_read()
    mol_names, var_names, x = input_sirms_file.read_next()

    while mol_names:

        x = adjust_dataset(x, var_names, ref_var_names)

        if ad is not None:
            for ad_name in ad:
                ad_dict[ad_name].extend(check_bound_box(x, bound_box_constrains))

        if scale is not None:
            x = scale.transform(x)

        for model_name, m in zip(pred_model_names, models):
            pred[model_name].extend(predict(x, m, model_name, model_type))

        mol_names, var_names, x = input_sirms_file.read_next()

    # convert to numpy array

    pred = np.array([pred[k] for k in pred_model_names]).transpose()

    # add consensuses of multiple models if we have them (svm_0, ...1, ...2..etc) and final consensus
    if len(model_names) < len(pred_model_names):  # multiple models
        for model_name in model_names:
            if model_type == "class":
                cols = [i for i, j in enumerate(pred_model_names) if re.fullmatch(model_name + "_\d*", j)]
                pred = np.append(pred, np.mean(np.around(pred[:, cols]), 1).reshape((-1, 1)), 1)
                pred_model_names.append(model_name + "_consensus")
        # add final consensus
        if model_type == "class":
            cols = [i for i, j in enumerate(pred_model_names) if j.endswith("_consensus")]
            pred = np.append(pred, np.mean(np.around(pred[:, cols]), 1).reshape((-1, 1)), 1)
            pred_model_names.append("consensus")
    else:
        if model_type == "class":
            pred = np.append(pred, np.mean(np.around(pred), 1).reshape((-1, 1)), 1)
            pred_model_names.append("consensus")
        elif model_type == "reg":
            pred = np.append(pred, np.mean(pred, 1).reshape((-1, 1)), 1)
            pred_model_names.append("consensus")
    pred = np.around(pred, 3)
    save_predictions(out_fname, pred, pred_model_names, input_sirms_file.get_mol_names(), ad_dict)


def main():

    parser = argparse.ArgumentParser(description='Predict properties of compounds based on QSAR models.')
    parser.add_argument('-i', '--input', metavar='descriptors.txt', required=True,
                        help='input text file with descriptors for compounds with removed fragments.')
    parser.add_argument('-f', '--input_format', metavar='txt|svm', required=False, default='txt',
                        help='format of the input file with descriptors (txt|svm). Default: txt.')
    parser.add_argument('-o', '--output', metavar='predictions.txt', required=True,
                        help='output text file with predicitions. Contains predictions of all models and their consensus. Special case:if' ' there are multiple models of individual'
                             'ML methods in the models_dir, then consensus for each metod and consensus of those consensuses (final consensus) will be computed.')
    parser.add_argument('-m', '--models', metavar='[rf gbm svr pls knn]', required=True, nargs='*',
                        help='file names of saved models (without extension).')
    parser.add_argument('-d', '--models_dir', metavar='path_to_models', required=True,
                        help='path to stored model files. The file scale.pkl must be in that dir also.')
    parser.add_argument('-t', '--model_type', metavar='reg|class', default='reg',
                        help='models type: reg for regression and class for classification')
    parser.add_argument('-a', '--applicability_domain', metavar='none|bound_box', required=False, nargs='*', default=None,
                        help='name(s) of applicability domain(s) to apply. If several - provide a space separated '
                             'list. Possible values: none - do not compute; bound_box - compounds with descriptor '
                             'values out of those for training set compounds are outside AD. Default: none')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='show progress on the screen.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": x_fname = v
        if o == "input_format": input_format = v
        if o == "output": out_fname = v
        if o == "models": model_names = v
        if o == "models_dir": model_dir = v
        if o == "model_type": model_type = v
        if o == "applicability_domain": ad = v
        if o == "verbose": verbose = v
    if ad is not None and 'none' in ad:
        ad.remove('none')
        if not ad:
            ad = None

    main_params(x_fname, input_format, out_fname, model_names, model_dir, model_type, ad, verbose)


if __name__ == '__main__':
    main()

