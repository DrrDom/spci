#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2014
# license         : LGPLv3
#==============================================================================

import os
import sys
import argparse
from datetime import datetime
from rdkit import Chem
from .find_frags_auto_rdkit import replace_no2


def load_smarts_file(fname):

    with open(fname) as f:
        frags = []
        for line in f:
            # skip comment lines and blank lines
            if line.strip().startswith('#') or not line.strip():
                continue
            tmp = line.strip().split("\t")
            if len(tmp) == 1 or (len(tmp) == 2 and tmp[1] == ""):
                frags.append({"name": tmp[0], "query": Chem.MolFromSmarts(tmp[0])})
            elif len(tmp) == 2:
                frags.append({"name": tmp[1], "query": Chem.MolFromSmarts(tmp[0])})
            elif len(tmp) == 3:
                frags.append({"name": tmp[1], "query": Chem.MolFromSmarts(tmp[0]),
                              "attached_atoms": list(map(int, tmp[2].split(",")))})
            else:
                print("Line " + line + " contains more than 3 tab-separated fields. Check file format.")
    return frags


def load_smi_fragments(fname):

    frags = []
    cansmi = set()

    with open(fname) as f:
        for line in f:
            tmp = line.strip().split()
            mol = Chem.MolFromSmiles(tmp[0])
            if mol is not None:
                smi = Chem.MolToSmiles(mol)
                if smi not in cansmi:
                    cansmi.add(smi)
                    frags.append({"name": smi, "query": mol})

    return frags


def load_sdf_fragments(fname):

    frags = []
    cansmi = set()

    for mol in Chem.SDMolSupplier(fname):
        if mol is not None:
            smi = Chem.MolToSmiles(mol)
            if smi not in cansmi:
                cansmi.add(smi)
                frags.append({"name": smi, "query": mol})

    return frags


def load_query_fragments(fname):

    frags = None

    if fname.endswith(".smi") or fname.endswith(".smiles"):
        frags = load_smi_fragments(fname)
    if fname.endswith(".sdf"):
        frags = load_sdf_fragments(fname)
    elif fname.endswith(".smart") or fname.endswith(".smarts") or fname.endswith(".txt"):
        frags = load_smarts_file(fname)
    else:
        print("Unsupported extension of fragments file. Use only smi or smiles for SMILES; "
              "sdf for SDF; smart, smarts, txt for SMARTS")

    return frags


def extend_atoms(mol, atom_ids, atom_numbers):

    def add_atoms(atom, atom_ids, atom_numbers):
        for nbr_atom in atom.GetNeighbors():
            if nbr_atom.GetIdx() not in atom_ids and nbr_atom.GetAtomicNum() in atom_numbers:
                atom_ids.append(nbr_atom.GetIdx())
                atom_ids = add_atoms(nbr_atom, atom_ids, atom_numbers)
        return atom_ids

    for atom in mol.GetAtoms():
        if atom.GetIdx() in atom_ids:
            atom_ids = add_atoms(atom, list(atom_ids), atom_numbers)

    return atom_ids


def get_map_atom_ids_list(mol, query, attached_atoms):
    """
    Returns list of tuples with atom ids if match and None otherwise
    ids started from 1, 2, 3, ...
    """

    m = Chem.AddHs(mol)
    ids = m.GetSubstructMatches(query)

    if attached_atoms is not None:
        ids = tuple(extend_atoms(m, i, attached_atoms) for i in ids)

    ids = tuple(tuple(sorted(j+1 for j in i)) for i in ids)   # 1-based ids 1, 2, 3, ...

    return ids if ids else None


def main_params(in_sdf, out_txt, in_frags, remove_all, verbose, error_fname):

    frags = load_query_fragments(in_frags)

    if frags is None:
        sys.exit()

    # find atoms ids which belongs to the fragments which should be removed
    # output example
    # {'paliperidone': {'hydroxy': [(12, 39)]},
    #  'salbutamol': {'hydroxy': [(16, 35), (13, 31), (15, 34)], 'amino': [(4, 24)]}, ...}

    with open(out_txt, "wt") as f:

        for i, mol in enumerate(Chem.SDMolSupplier(in_sdf, sanitize=False, removeHs=False)):
            if mol is not None:

                # modify representation of NO2 groups to charged version
                mol = replace_no2(mol)
                err = Chem.SanitizeMol(mol, catchErrors=True)
                if err != 0:
                    with open(error_fname, 'at') as f_err:
                        f_err.write('%s\t%s\t%s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                      os.path.basename(__file__),
                                                      'Molecule %s failed to sanitize due to: ' % mol.GetProp("_Name") + str(err)))

                mol_name = mol.GetProp("_Name")
                if verbose:
                    print("Searching for fragments in", mol_name)

                for frag in frags:
                    ids = get_map_atom_ids_list(mol, frag["query"], frag.get("attached_atoms"))
                    if ids:
                        if remove_all:
                            tmp = set()
                            for v in ids:
                                tmp.update(v)
                            f.write(mol_name + "\t" + frag["name"] + "\t" + "\t".join(map(str, tuple(tmp))) + "\n")
                        else:
                            for v in ids:
                                f.write(mol_name + "\t" + frag["name"] + "\t" + "\t".join(map(str, v)) + "\n")

            else:
                print('Molecules #%i was skipped due to error' % i)
                with open(error_fname, 'at') as f_err:
                    f_err.write('%s\t%s\t%i\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                  os.path.basename(__file__), i + 1))


def main():

    parser = argparse.ArgumentParser(description='Find indices of fragments atoms in molecules from input sdf-file.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf file with standardized structures, molecules should have titles.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file; each line contains tab-separated fields: name of a molecule, '
                             'name of a fragment and list of corresponding atom numbers.')
    parser.add_argument('-f', '--frag', metavar='fragments.smarts', required=True,
                        help='file containing fragments in SMARTS(.txt|.smarts)/SMILES(.smi)/SDF(.sdf) format. '
                             'Tab-separated SMARTS/SMILES file contains 1) SMARTS/SMILES, 2) fragment name (optional), '
                             '3) comma-separated list of atom numbers for recursive addition (optional). '
                             'Comment lines start with #')
    parser.add_argument('-a', '--all', action='store_true', default=False,
                        help='if set this flag to true all identical fragments occurred more than once '
                             'in a molecule will be removed simultaneously.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='show progress on the screen.')
    parser.add_argument('-e', '--error_file', metavar='log_file_name.txt', default="indigo_errors.txt",
                        help='save names of molecules which cause error to a text log file. Default file name '
                             'indigo_errors.txt.')
    # parser.add_argument('-d', '--allow_duplicate_fragments', action='store_true', default=False,
    #                     help='pre-filtering of duplicate fragments. Default: false (duplicate fragments not allowed).')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_sdf = v
        if o == "out": out_txt = v
        if o == "frag": in_frags = v
        if o == "all": remove_all = v
        if o == "verbose": verbose = v
        if o == "error_file": error_fname = v
        # if o == "allow_duplicate_fragments": dupl = v

    main_params(in_sdf, out_txt, in_frags, remove_all, verbose, error_fname)


if __name__ == '__main__':
    main()
