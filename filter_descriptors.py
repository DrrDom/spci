#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 19-01-2016
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2016
# license         : GPL3
#==============================================================================

import os
import argparse
import numpy as np
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


def save_sirms_txt(x, descr_names, case_names, fname):
    print(np.shape(x))
    print(len(descr_names))
    print(len(case_names))
    with open(fname, 'wt') as f:
        f.write('Compounds' + '\t' + '\t'.join(descr_names) + '\n')
        for i, row in enumerate(x):
            f.write(case_names[i] + '\t' + '\t'.join(map(str, row.tolist())) + '\n')


def save_sirms_svm(x, descr_names, case_names, fname):
    open(os.path.splitext(fname)[0] + '.rownames', 'wt').write('\n'.join(case_names))
    open(os.path.splitext(fname)[0] + '.colnames', 'wt').write('\n'.join(descr_names))
    with open(fname, 'wt') as f:
        for row in x:
            ids = np.nonzero(row)[0].tolist()
            row = row[row != 0].tolist()
            line = [str(i) + ':' + str(v) for i, v in zip(ids, row)]
            f.write(' '.join(line) + '\n')


def main_params(in_fname, out_fname, input_format):

    if input_format == 'txt':
        descr_names, mol_names, x = load_sirms_txt(in_fname)
    elif input_format == 'svm':
        descr_names, mol_names, x = load_sirms_svm(in_fname)

    ids = [i for i, el in enumerate(descr_names) if "|HB|I,I,I,I|" in el]
    x = np.delete(x, ids, 1)
    descr_names = [v for i, v in enumerate(descr_names) if i not in ids]

    if input_format == 'txt':
        save_sirms_txt(x, descr_names, mol_names, out_fname)
    elif input_format == 'svm':
        save_sirms_svm(x, descr_names, mol_names, out_fname)


def main():

    parser = argparse.ArgumentParser(description='Filter descriptors: discard HB descriptors with all atoms labeled I.')
    parser.add_argument('-i', '--input', metavar='input_descriptors.txt', required=True,
                        help='input file with descriptors in txt or svm formats.')
    parser.add_argument('-o', '--output', metavar='output_descriptors.txt', default=None,
                        help='output file with descriptors in txt or svm formats. If missing then output will be saved '
                             'into the input file.')
    parser.add_argument('-f', '--file_format', metavar='txt|svm', default='txt',
                        help='input file format tt or svm. Default: txt.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": out_fname = v
        if o == "file_format": input_format = v
    if input_format not in ['txt', 'svm']:
        print("Input file format is wrong - %s. Only txt and svm are allowed." % input_format)
        exit()
    if out_fname is None:
        out_fname = in_fname

    main_params(in_fname, out_fname, input_format)


if __name__ == '__main__':
    main()
