#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 19-01-2016
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2016
# license         : LGPLv3
#==============================================================================

import os
import argparse
import shutil
# import numpy as np
# from scipy.sparse import dok_matrix


# def load_sirms_txt(fname):
#     with open(fname) as f:
#         descr_names = f.readline().strip().split("\t")[1:]
#         case_names = []
#         x = []
#         for line in f:
#             tmp = line.strip().split("\t")
#             case_names.append(tmp[0])
#             x.append(tuple(map(float, tmp[1:])))
#     return descr_names, case_names, np.asarray(x)
#
#
# def load_sirms_svm(fname):
#     descr_names = [v.strip() for v in open(os.path.splitext(fname)[0] + '.colnames').readlines()]
#     case_names = [v.strip() for v in open(os.path.splitext(fname)[0] + '.rownames').readlines()]
#     x = dok_matrix((len(case_names), len(descr_names)), dtype=np.float32)
#     with open(fname) as f:
#         for row, line in enumerate(f):
#             tmp = line.strip().split(' ')
#             for v in tmp:
#                 col, value = v.split(':')
#                 x[row, int(col)] = value
#     return descr_names, case_names, x.toarray()
#
#
# def save_sirms_txt(x, descr_names, case_names, fname):
#     print(np.shape(x))
#     print(len(descr_names))
#     print(len(case_names))
#     with open(fname, 'wt') as f:
#         f.write('Compounds' + '\t' + '\t'.join(descr_names) + '\n')
#         for i, row in enumerate(x):
#             f.write(case_names[i] + '\t' + '\t'.join(map(str, row.tolist())) + '\n')
#
#
# def save_sirms_svm(x, descr_names, case_names, fname):
#     open(os.path.splitext(fname)[0] + '.rownames', 'wt').write('\n'.join(case_names))
#     open(os.path.splitext(fname)[0] + '.colnames', 'wt').write('\n'.join(descr_names))
#     with open(fname, 'wt') as f:
#         for row in x:
#             ids = np.nonzero(row)[0].tolist()
#             row = row[row != 0].tolist()
#             line = [str(i) + ':' + str(v) for i, v in zip(ids, row)]
#             f.write(' '.join(line) + '\n')


def main_params(in_fname, out_fname, file_format):

    if file_format == 'txt':

        with open(in_fname) as f_in:
            with open(in_fname + '.tmp' if out_fname == in_fname else out_fname, 'wt') as f_out:
                header = f_in.readline().strip().split('\t')
                ids = ["|HB|I,I,I,I|" not in item for item in header]
                f_out.write('\t'.join([item for item, id in zip(header, ids) if id]) + '\n')
                for line in f_in:
                    f_out.write('\t'.join([item for item, id in zip(line.strip().split('\t'), ids) if id]) + '\n')
        if out_fname == in_fname:
            shutil.move(in_fname + '.tmp', in_fname)

    elif file_format == 'svm':

        colnames = open(os.path.splitext(in_fname)[0] + '.colnames').readlines()
        new_colnames = [item for item in colnames if "|HB|I,I,I,I|" not in item]

        ids = dict()  # {old_num: new_num, ...} as strings
        for i in range(len(colnames)):
            try:
                ids[str(i)] = str(new_colnames.index(colnames[i]))
            except ValueError:
                continue

        with open(in_fname) as f_in:
            with open(in_fname + '.tmp' if out_fname == in_fname else out_fname, 'wt') as f_out:
                for line in f_in:
                    items = [tuple(item.split(':')) for item in line.strip().split(' ')]
                    output = []
                    for item in items:
                        id = ids.get(item[0], None)
                        if id is not None:
                            output.append(id + ':' + item[1])
                    f_out.write(' '.join(output) + '\n')

        if out_fname == in_fname:
            open(os.path.splitext(in_fname)[0] + '.colnames', 'wt').writelines(new_colnames)
            shutil.move(in_fname + '.tmp', in_fname)
        else:
            open(os.path.splitext(out_fname)[0] + '.colnames', 'wt').writelines(new_colnames)
            shutil.copyfile(os.path.splitext(in_fname)[0] + '.rownames', os.path.splitext(out_fname)[0] + '.rownames')


def main():

    parser = argparse.ArgumentParser(description='Filter descriptors: discard HB descriptors with all atoms labeled I.')
    parser.add_argument('-i', '--input', metavar='input_descriptors.txt', required=True,
                        help='input file with descriptors in txt or svm formats.')
    parser.add_argument('-o', '--output', metavar='output_descriptors.txt', default=None,
                        help='output file with descriptors in txt or svm formats. If missing then input file will be '
                             'rewritten through temporary file.')
    parser.add_argument('-f', '--file_format', metavar='txt|svm', default='txt',
                        help='input file format txt or svm. Default: txt.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": out_fname = v
        if o == "file_format": file_format = v
    if file_format not in ['txt', 'svm']:
        print("Input file format is wrong - %s. Only txt and svm are allowed." % file_format)
        exit()
    if out_fname is None:
        out_fname = in_fname

    main_params(in_fname, out_fname, file_format)


if __name__ == '__main__':
    main()
