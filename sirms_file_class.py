#!/usr/bin/env python
# author          : Pavel
# date            : 05.11.15
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel 2015
# license         : GPL3
#==============================================================================

import numpy as np


class SirmsFile():

    def __init__(self, fname, frag_file=False, chunks=float("Inf")):
        self.__varnames = open(fname).readline().strip().split('\t')[1:]
        self.__file = open(fname)
        self.__frag_full_names = []          # keep the names of all fragments which were read
                                             # order is important since calculated contributions are saved in
                                             # the same order
        self.__is_frag_full_names_read = False
        self.__nlines = chunks               # file is read by chunks to avoid memory problems,
                                             # to read the whole file at once set to float('Inf')
        self.__is_frag_file = frag_file      # set which file to read ordinary (False) or with fragments (True)
                                             # in the latter case each chunk will contain all fragments for each
                                             # read molecules (another reading process)

    def reset_read(self):
        self.__file.seek(0)

    def read_next(self):

        # skip header
        if self.__file.tell() == 0:
            self.__file.readline()

        lines = []

        # read first nlines
        i = 0
        while True:
            line = self.__file.readline()
            if not line:
                break
            lines.append(line.strip().split('\t'))
            i += 1
            if i >= self.__nlines:
                break

        # read next lines until meet a new molecule
        if self.__is_frag_file and lines:
            while True:
                cur_mol = lines[-1][0].split("#")[0]
                cur_pos = self.__file.tell()
                line = self.__file.readline()
                if not line:
                    break
                if line.split('\t', 1)[0].split('#')[0] == cur_mol:
                    lines.append(line.strip().split('\t'))
                else:
                    self.__file.seek(cur_pos)
                    break

        mol_names = []
        x = []
        for line in lines:
            mol_names.append(line[0])
            x.append(tuple(map(float, line[1:])))

        # add read frag names to the list (remain only mol name and frag name, e.g. mol1#frag1, not mol1#frag1#1)
        if not self.__is_frag_full_names_read:
            self.__frag_full_names.extend(['#'.join(mol_name.split('#')[0:2]) for mol_name in mol_names if mol_name.find('#') > -1])

        # if EOF then stop updating the list of fragments (it is read only once)
        if not self.__is_frag_full_names_read and len(mol_names) < self.__nlines:
            self.__is_frag_full_names_read = True

        return mol_names, self.__varnames, np.asarray(x)
