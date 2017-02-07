#!/usr/bin/env python
# author          : Pavel
# date            : 05.11.15
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel 2015
# license         : GPL3
#==============================================================================

import os
import argparse
import numpy as np
from sklearn.externals import joblib
from collections import defaultdict

from sirms_file_class import SirmsFile


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
        if model_name == 'pls':
            pred = pred.reshape(pred.shape[0] * pred.shape[1])
        pred = pred.tolist()

    elif model_type == 'class':
        pos_id = model.classes_.tolist().index(1)
        pred = [i[pos_id] for i in model.predict_proba(x)]

    return pred


def add_consensus(pred_array, pred_model_names):
    if np.shape(pred_array)[1] > 1:
        pred_array = np.append(pred_array, np.mean(pred_array, 1).reshape((-1, 1)), 1)
        pred_model_names.append('consensus')
    return pred_array, pred_model_names


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
    # return np.all(np.logical_and(np.greater_equal(x, bound_box_constrains[0]),
    #                              np.less_equal(x, bound_box_constrains[1])),
    #               axis=1).tolist()


def main_params(x_fname, input_format, out_fname, train_x_fname, train_format, model_names, model_dir, model_type, ad, verbose):

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
    models = dict()
    for m in model_names:
        models[m] = joblib.load(os.path.join(model_dir, m + ".pkl"))

    input_sirms_file = SirmsFile(x_fname, file_format=input_format, chunks=1000)
    ref_var_names = get_var_names(train_x_fname, train_format)

    input_sirms_file.reset_read()
    mol_names, var_names, x = input_sirms_file.read_next()

    while mol_names:

        x = adjust_dataset(x, var_names, ref_var_names)

        if ad is not None:
            for ad_name in ad:
                ad_dict[ad_name].extend(check_bound_box(x, bound_box_constrains))

        if scale is not None:
            x = scale.transform(x)

        for model_name, m in models.items():
            pred[model_name].extend(predict(x, m, model_name, model_type))

        mol_names, var_names, x = input_sirms_file.read_next()

    # convert to numpy array
    pred_model_names = list(sorted(pred.keys()))
    pred = np.array([pred[k] for k in pred_model_names]).transpose()

    pred, pred_model_names = add_consensus(pred, pred_model_names)
    save_predictions(out_fname, pred, pred_model_names, input_sirms_file.get_mol_names(), ad_dict)


def main():

    parser = argparse.ArgumentParser(description='Predict properties of compounds based on QSAR models.')
    parser.add_argument('-i', '--input', metavar='descriptors.txt', required=True,
                        help='input text file with descriptors for compounds with removed fragments.')
    parser.add_argument('-f', '--input_format', metavar='txt|svm', required=False, default='txt',
                        help='format of the input file with descriptors (txt|svm). Default: txt.')
    parser.add_argument('-o', '--output', metavar='predictions.txt', required=True,
                        help='output text file with fragments contributions.')
    parser.add_argument('-r', '--training_set', metavar='training_set_descriptors.txt', required=True,
                        help='text file with descriptors nad their names of training set compounds.')
    parser.add_argument('--train_format', metavar='txt|svm', required=False, default='txt',
                        help='format of the input file with descriptors (txt|svm). Default: txt.')
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
        if o == "training_set": train_x_fname = v
        if o == "train_format": train_format = v
        if o == "models": model_names = v
        if o == "models_dir": model_dir = v
        if o == "model_type": model_type = v
        if o == "applicability_domain": ad = v
        if o == "verbose": verbose = v
    if ad is not None and 'none' in ad:
        ad.remove('none')
        if not ad:
            ad = None

    main_params(x_fname, input_format, out_fname, train_x_fname, train_format, model_names, model_dir, model_type, ad, verbose)


if __name__ == '__main__':
    main()
