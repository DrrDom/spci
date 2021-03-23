#!/usr/bin/env python
# author          : Pavel
# date            : 05.11.15
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel 2015
# license         : LGPLv3
#==============================================================================

import os
import numpy as np


mol_frag_sep = "###"


class SirmsFile():

    def __init__(self, fname, file_format='txt', frag_file=False, chunks=float("Inf")):

        if file_format == 'txt':
            self.__varnames = open(fname).readline().strip().split('\t')[1:]
            self.__mol_full_names = []           # keep the names of all molecules which were read
            self.__frag_full_names = []          # keep the names of all fragments which were read
                                                 # order is important since calculated contributions are saved in
                                                 # the same order (mol1###frag1#1 <- frag id is present)
            self.__is_mol_full_names_read = False
            self.__is_frag_full_names_read = False

        elif file_format == 'svm':
            self.__varnames = [v.strip() for v in open(os.path.splitext(fname)[0] + '.colnames').readlines()]
            self.__mol_full_names = [v.strip() for v in open(os.path.splitext(fname)[0] + '.rownames').readlines()]
            self.__frag_full_names = [v for v in self.__mol_full_names if v.find(mol_frag_sep) > -1]
            self.__is_frag_full_names_read = True

        self.__file = open(fname)
        self.__file_format = file_format
        self.__nlines = chunks               # file is read by chunks to avoid memory problems,
                                             # to read the whole file at once set to float('Inf')
        self.__is_frag_file = frag_file      # set which file to read ordinary (False) or with fragments (True)
                                             # in the latter case each chunk will contain all fragments for each
                                             # read molecules (another reading process)
        self.__cur_mol_read = 0              # number of mols (lines) already read (or index of line to read)

    def get_frag_full_names(self):
        return self.__frag_full_names

    def get_mol_names(self):
        return self.__mol_full_names

    def reset_read(self):
        self.__cur_mol_read = 0
        self.__file.seek(0)

    def __read_svm_next(self):
        if self.__cur_mol_read == len(self.__mol_full_names):
            return [], self.__varnames, np.asarray([])
        else:
            start = self.__cur_mol_read
            end = start + self.__nlines
            if end > len(self.__mol_full_names):
                end = len(self.__mol_full_names)
            else:
                cur_mol = self.__mol_full_names[end - 1].split(mol_frag_sep)[0]
                # end += 1
                while end < len(self.__mol_full_names) and self.__mol_full_names[end].split(mol_frag_sep)[0] == cur_mol:
                    end += 1
            # end -= 1
            lines = [self.__file.readline().strip() for _ in range(start, end)]
            x = []
            for line in lines:
                out = [0] * len(self.__varnames)
                for entry in line.split():
                    index, value = entry.split(':')
                    out[int(index)] = float(value)
                x.append(out)
            x = np.asarray(x)
            self.__cur_mol_read = end
            return self.__mol_full_names[start:end], self.__varnames, x

    def __read_txt_next(self):

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
                cur_mol = lines[-1][0].split(mol_frag_sep)[0]
                cur_pos = self.__file.tell()
                line = self.__file.readline()
                if not line:
                    break
                if line.split('\t', 1)[0].split(mol_frag_sep)[0] == cur_mol:
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
            self.__frag_full_names.extend([mol_name for mol_name in mol_names if mol_name.find(mol_frag_sep) > -1])

        if not self.__is_mol_full_names_read:
            self.__mol_full_names.extend(mol_names)

        # if EOF then stop updating the list of fragments and mol names (it is read only once)
        if len(mol_names) < self.__nlines:
            self.__is_frag_full_names_read = True
            self.__is_mol_full_names_read = True

        return mol_names, self.__varnames, np.asarray(x)

    def read_next(self):
        if self.__file_format == 'txt':
            return self.__read_txt_next()
        elif self.__file_format == 'svm':
            return self.__read_svm_next()
