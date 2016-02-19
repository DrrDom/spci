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
import platform
import numpy as np
from sklearn.externals import joblib
from collections import defaultdict, OrderedDict

from sirms_file_class import SirmsFile

mol_frag_sep = "###"


def load_sirms(fname):
    with open(fname) as f:
        descr_names = f.readline().strip().split("\t")[1:]
        case_names = []
        x = []
        for line in f:
            tmp = line.strip().split("\t")
            case_names.append(tmp[0])
            x.append(tuple(map(float, tmp[1:])))
    return descr_names, case_names, np.asarray(x)


def split_mol_frag_names(names, ids):
    """
    names: list of names in format mol_name#frag_name
    ids: list of boolean
    """
    mol_names = []
    frag_names = []
    for i, el in enumerate(ids):
        if el:
            tmp = names[i].split(mol_frag_sep)
            mol_names.append(tmp[0])
            frag_names.append(tmp[1])
    return mol_names, frag_names


def prepare_dataset(x_train, x_frag, prop_name, descr_names, x_train_mol_names, x_frag_mol_names):
    """
    x_train: numpy array with descriptors for training set
    x_frag: numpy array with descriptors for dataset after removal of specified fragments
    prop_name: name of atomic property, which contribution is estimated
    descr_names: list of descriptors names in x_train and x_frag
    x_train_mol_names: list of molecule names in x_train
    x_frag_mol_names: list of molecule names in x_frag

    http://stackoverflow.com/questions/26086645/manipulation-with-arrays-in-numpy
    """
    ids = np.array(["|" + prop_name + "|" in el for el in descr_names])
    rows = np.where(np.array(x_train_mol_names) == np.array(x_frag_mol_names)[:, np.newaxis])[1]
    output = np.where(ids, x_frag, x_train[rows])
    return output


def save_contrib(fname, contrib_dict, frag_full_names, long_format):
    if not long_format:
        with open(fname, "wt") as f:
            f.write("Model_Contribution\t" + "\t".join(frag_full_names) + "\n")
            for model_name in contrib_dict.keys():
                for prop_name in contrib_dict[model_name].keys():
                    f.write(model_name + "_" + prop_name + "\t" +
                            "\t".join(["{0:.6f}".format(i) for i in contrib_dict[model_name][prop_name]]) + "\n")
    else:
        with open(fname, "wt") as f:
            f.write("Compound\tFragment\tModel\tContribution_type\tContribution_value\n")
            for model_name in contrib_dict.keys():
                for prop_name in contrib_dict[model_name].keys():
                    for frag_full_name, value in zip(frag_full_names, contrib_dict[model_name][prop_name]):
                        f.write('\t'.join(frag_full_name.split(mol_frag_sep)) + '\t' + model_name + '\t' + prop_name + '\t' + "{0:.6f}".format(value) + '\n')


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


def main_params(x_fname, out_fname, model_names, model_dir, prop_names, model_type, verbose, save_pred, input_format, long_format):

    if save_pred:
        save_pred_fname = os.path.splitext(out_fname)[0] + "_pred.txt"
        if os.path.isfile(save_pred_fname):
            print("File with prediction was erased.")
            os.remove(save_pred_fname)

    sirms_file = SirmsFile(x_fname, file_format=input_format, frag_file=True, chunks=1000)

    # scale
    if os.path.isfile(os.path.join(model_dir, "scale.pkl")):
        scale = joblib.load(os.path.join(model_dir, "scale.pkl"))
    else:
        scale = None

    # frag_contrib = dict()
    frag_contrib = OrderedDict()

    if model_type == 'reg' or model_type == 'class':

        for model_name in model_names:

            if verbose:
                print("Contributions calculation. %s model proceed" % model_name)

            frag_contrib[model_name] = defaultdict(list)

            m = joblib.load(os.path.join(model_dir, model_name + ".pkl"))

            sirms_file.reset_read()
            mol_names, var_names, x = sirms_file.read_next()

            while mol_names:

                if scale is not None:
                    x = scale.transform(x)

                x_frag_ids = [mol_frag_sep in n for n in mol_names]
                x_train_ids = [not i for i in x_frag_ids]

                x_frag_mol_names, x_frag_frag_names = split_mol_frag_names(mol_names, x_frag_ids)
                x_train_mol_names = [mol_names[i] for i, el in enumerate(x_train_ids) if el]

                train_pred = predict(x[np.asarray(x_train_ids), ], m, model_name, model_type)
                train_pred = dict(zip(x_train_mol_names, train_pred))

                for prop_name in prop_names:
                    if prop_name == "overall":
                        frag_pred = predict(x[np.asarray(x_frag_ids), ], m, model_name, model_type)
                    else:
                        x_frag_prep = prepare_dataset(x[np.asarray(x_train_ids), ],
                                                      x[np.asarray(x_frag_ids), ],
                                                      prop_name,
                                                      var_names,
                                                      x_train_mol_names,
                                                      x_frag_mol_names)
                        frag_pred = predict(x_frag_prep, m, model_name, model_type)

                    frag_contrib[model_name][prop_name].extend([train_pred[mol_name] - frag_pred[i] for i, mol_name in enumerate(x_frag_mol_names)])

                    if save_pred:
                        with open(save_pred_fname, "at") as f:
                            for k, v in train_pred.items():
                                f.write(k + "\t" + model_name + "\t" + prop_name + "\t" + str(v) + "\n")
                            for m_name, f_name, pred_value in zip(x_frag_mol_names, x_frag_frag_names, frag_pred):
                                f.write(m_name + mol_frag_sep + f_name + "\t" + model_name + "\t" + prop_name + "\t" + str(pred_value) + "\n")

                mol_names, var_names, x = sirms_file.read_next()

    # elif model_type == 'class':
    #
    #     for model_name in model_names:
    #
    #         frag_contrib[model_name] = dict()
    #         m = joblib.load(os.path.join(model_dir, model_name + ".pkl"))
    #
    #         pos_id = m.classes_.tolist().index(1)
    #         xtrain_pred = [i[pos_id] for i in m.predict_proba(x_train)]
    #         xtrain_pred = dict(zip(x_train_mol_names, xtrain_pred))
    #
    #         for prop_name in prop_names:
    #             if prop_name != "overall":
    #                 x_frag_prep = prepare_dataset(x_train, x_frag, prop_name, descr_names, x_train_mol_names,
    #                                               x_frag_mol_names)
    #                 xfrag_pred = [i[pos_id] for i in m.predict_proba(x_frag_prep)]
    #             else:
    #                 xfrag_pred = [i[pos_id] for i in m.predict_proba(x_frag)]
    #
    #             frag_contrib[model_name][prop_name] = [xtrain_pred[mol_name] - xfrag_pred[i] for i, mol_name in enumerate(x_frag_mol_names)]

    save_contrib(out_fname, frag_contrib, sirms_file.get_frag_full_names(), long_format)


def main():

    parser = argparse.ArgumentParser(description='Calculate contributions of molecular fragments based on QSAR models.')
    parser.add_argument('-i', '--input', metavar='frag_descriptors.txt', required=True,
                        help='input text file with descriptors for compounds with removed fragments.')
    parser.add_argument('-f', '--input_format', metavar='txt', default='txt',
                        help='format of input file with descriptors for compounds with removed fragments (txt|svm). '
                             'Svm sparse format file ha to have two files with the same name and extensions '
                             'rownames and colnames. Default: txt.')
    parser.add_argument('-o', '--output', metavar='contributions.txt', required=True,
                        help='output text file with fragments contributions.')
    parser.add_argument('--output_long_format', action='store_true', default=False,
                        help='store output in long instead of short format. Long format is easier to parse but takes '
                             'more disk space. Default: false.')
    parser.add_argument('-m', '--models', metavar='[rf gbm svr pls knn]', required=True, nargs='*',
                        help='file name of saved models (without extension)')
    parser.add_argument('-d', '--models_dir', metavar='path_to_models', required=True,
                        help='path to stored model files')
    parser.add_argument('-p', '--properties', metavar='[overall charge hb ...]', nargs='*', default=["overall"],
                        help='name of atomic properties which should be considered during calculation of '
                             'contribution of fragments')
    parser.add_argument('-t', '--model_type', metavar='reg|class', default='reg',
                        help='models type: reg for regression and class for classification')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='show progress on the screen.')
    parser.add_argument('-s', '--save_intermediate_prediction', action='store_true', default=False,
                        help='save intermediate prediction to text file.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": x_fname = v
        if o == "output": out_fname = v
        if o == "models": model_names = v
        if o == "models_dir": model_dir = v
        if o == "properties": prop_names = v
        if o == "model_type": model_type = v
        if o == "verbose": verbose = v
        if o == "save_intermediate_prediction": save_pred = v
        if o == "input_format": input_format = v
        if o == "output_long_format": long_format = v
    if input_format not in ['txt', 'svm']:
        print("Wrong input file format - %s. Only txt and svm file formats are allowed." % input_format)
        exit()

    main_params(x_fname, out_fname, model_names, model_dir, prop_names, model_type, verbose, save_pred, input_format, long_format)


if __name__ == '__main__':
    main()
