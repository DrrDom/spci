#!/usr/bin/env python
# author          : Pavel
# date            : 04.04.17
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2017
# license         : LGPLv3
#==============================================================================

import os
import sys
import argparse
from datetime import datetime
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from .find_frags_auto_rdkit import replace_no2


def add_atoms(atom, remove_atoms):
    for a in atom.iterateNeighbors():
        d = set([ai.index() for ai in a.iterateNeighbors()]).difference(remove_atoms)
        if len(d) == 1:
            remove_atoms.append(a.index())
            add_atoms(a, remove_atoms)


def main_params(in_sdf, out_txt, verbose, error_fname):

    with open(out_txt, "wt") as f:

        for i, mol in enumerate(Chem.SDMolSupplier(in_sdf, sanitize=False, removeHs=False)):

            if mol:

                # modify representation of NO2 groups to charged version
                mol = replace_no2(mol)
                err = Chem.SanitizeMol(mol, catchErrors=True)
                if err != 0:
                    with open(error_fname, 'at') as f_err:
                        f_err.write('%s\t%s\t%s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                      os.path.basename(__file__),
                                                      'Molecule %s failed to sanitize due to: ' % mol.GetProp("_Name") + str(err)))

                s = MurckoScaffold.GetScaffoldForMol(mol)
                ids = mol.GetSubstructMatches(s)
                if len(ids) == 1:
                    if len(ids[0]) > 0:   # to skip mol without scaffold
                        ids = [j + 1 for j in ids[0]]
                        f.write(mol.GetProp("_Name") + "\t" + Chem.MolToSmiles(s) + "\t" + "\t".join(map(str, ids)) + "\n")

                else:
                    sys.stderr.write('More than one scaffold match was found. Molecule: %s, scaffold: %s' %
                                     (Chem.MolToSmiles(mol), Chem.MolToSmiles(s)))
                    sys.stderr.flush()

            else:
                with open(error_fname, 'at') as f_err:
                    f_err.write('%s\t%s\tcould not read molecule number %i from file\n' %
                                (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                 os.path.basename(__file__), i + 1))

            if verbose and (i + 1) % 100 == 0:
                print('%i molecules passed' % (i + 1))


def main():

    parser = argparse.ArgumentParser(description='Find indices of Murcko framework atoms in molecules from input sdf-file.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf file with standardized structures, molecules should have titles.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file; each line contains tab-separated fields: name of a molecule, '
                             'name of a fragment and list of corresponding atom numbers.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='show progress on the screen.')
    parser.add_argument('-e', '--error_file', metavar='log_file_name.txt', default="rdkit_errors.txt",
                        help='save names of molecules which cause error to a text log file. Default file name '
                             'rdkit_errors.txt.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_sdf = v
        if o == "out": out_txt = v
        if o == "verbose": verbose = v
        if o == "error_file": error_fname = v

    main_params(in_sdf, out_txt, verbose, error_fname)


if __name__ == '__main__':
    main()
