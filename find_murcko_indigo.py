#!/usr/bin/env python
# author          : Pavel
# date            : 22.06.15
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2015
# license         : GPL3
#==============================================================================

import os
import argparse
from datetime import datetime
from indigo import Indigo, IndigoException


indigo = Indigo()


def add_atoms(atom, remove_atoms):
    for a in atom.iterateNeighbors():
        d = set([ai.index() for ai in a.iterateNeighbors()]).difference(remove_atoms)
        if len(d) == 1:
            remove_atoms.append(a.index())
            add_atoms(a, remove_atoms)


def main_params(in_sdf, out_txt, verbose, error_fname):

    with open(out_txt, "wt") as f:

        for mol in indigo.iterateSDFile(in_sdf):

            try:

                mol.dearomatize()  # in order to avoid different aromatic forms in different molecules
                mol_name = mol.name()
                if verbose:
                    print("Searching for Murcko frameworks in", mol_name)

                # return Murcko framework
                remove_atoms = []
                for a in mol.iterateAtoms():
                    if a.degree() == 1:
                        remove_atoms.append(a.index())
                        add_atoms(a, remove_atoms)

                murcko_ids = list(set(range(mol.countAtoms())).difference(remove_atoms))

                # if a murcko framework is available in a molecule than write it
                if len(murcko_ids) > 0:
                    murcko_smiles = mol.createSubmolecule(murcko_ids).canonicalSmiles().split(" ")[0]
                    murcko_ids = list(map(lambda x: x + 1, murcko_ids))
                    f.write(mol_name + "\t" + murcko_smiles + "\t" + "\t".join(map(str, murcko_ids)) + "\n")

            except IndigoException as e:

                print('%s was skipped due to error' % mol.name())
                print(e)
                with open(error_fname, 'at') as f_err:
                    f_err.write('%s\t%s\t%s\t%s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                      os.path.basename(__file__), mol.name(), e))




def main():

    parser = argparse.ArgumentParser(description='Find indices of Murcko framework atoms in molecules from input sdf-file.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf file with standardized structures, molecules should have titles.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file; each line contains tab-separated fields: name of a molecule, '
                             'name of a fragment and list of corresponding atom numbers.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='show progress on the screen.')
    parser.add_argument('-e', '--error_file', metavar='log_file_name.txt', default="indigo_errors.txt",
                        help='save names of molecules which cause error to a text log file. Default file name '
                             'indigo_errors.txt.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_sdf = v
        if o == "out": out_txt = v
        if o == "verbose": verbose = v
        if o == "error_file": error_fname = v

    main_params(in_sdf, out_txt, verbose, error_fname)


if __name__ == '__main__':
    main()
