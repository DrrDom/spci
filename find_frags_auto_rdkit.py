#!/usr/bin/env python
# author          : Pavel
# date            : 04.12.15
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2015
# license         : GPL3
#==============================================================================

import os
import argparse

from itertools import combinations
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import AllChem


def replace_no2(mol):
    query = Chem.MolFromSmarts('n(:o):o')
    repl = Chem.MolFromSmiles('[N+](=O)[O-]')
    return AllChem.ReplaceSubstructs(mol, query, repl, replaceAll=True)[0]


def fix_cansmi_attach_point(smi):
    # fix cansmi name by replacement of attachment points: Cl%91.[*:1]%91 Cl[*:1]
    p = smi.split(".")
    output = p[0]
    for i in p[1:]:
        att = i.split('%')
        output = output.replace('%' + att[1], "(" + att[0] + ")")
    return output


def frag_mol_by_cuts(mol, cut_list):

    em = Chem.EditableMol(mol)

    output = []
    ap_ids = []

    for cut in cut_list:
        em.RemoveBond(cut[0], cut[1])
        ap1 = em.AddAtom(Chem.Atom(0))
        ap_ids.append(ap1)
        em.AddBond(cut[0], ap1, Chem.BondType.SINGLE)
        ap2 = em.AddAtom(Chem.Atom(0))
        ap_ids.append(ap2)
        em.AddBond(cut[1], ap2, Chem.BondType.SINGLE)

    mol = em.GetMol()

    for ids in Chem.GetMolFrags(mol):
        cansmi = Chem.MolFragmentToSmiles(mol, atomsToUse=ids)
        # cansmi = fix_cansmi_attach_point(cansmi)
        ids = [i for i in ids if i not in ap_ids]
        output.append([cansmi, ids])

    return output


def filter_dupl(frag_list):
    # input format: [['F', [0]], ['C#N', [3, 4]], ... ]
    res = dict()
    for item in frag_list:
        res[frozenset(sorted(item[1]))] = item[0]
    return [[v, list(k)] for k, v in res.items()]


def fragment_mol(mol, query, max_cuts):
    # returns list of lists: [['F', [0]], ['C#N', [3, 4]], ... ]

    # mol = Chem.AddHs(mol)

    # modify representation of NO2 groups to charged version
    mol = replace_no2(mol)
    err = Chem.SanitizeMol(mol, catchErrors=True)
    if err != 0:
        print('Molecule %s failed to sanitize due to: ' % mol.GetProp("_Name") + str(err))
        return []

    output = []

    all_cuts = mol.GetSubstructMatches(query)

    for i in range(1, max_cuts + 1):
        for comb in combinations(all_cuts, i):
            output.extend(frag_mol_by_cuts(mol, comb))

    output = filter_dupl(output)

    return output


def main_params(in_sdf, out_txt, query, max_cuts, verbose, error_fname):

    query = Chem.MolFromSmarts(query)

    with open(out_txt, 'wt') as f:
        for i, mol in enumerate(Chem.SDMolSupplier(in_sdf, sanitize=False, removeHs=False)):
            if mol is not None:
                if verbose:
                    print(mol.GetProp("_Name") + ' is processing')
                res = fragment_mol(mol, query, max_cuts)
                for item in res:
                    ids = [i + 1 for i in item[1]]  # save as 1-based ids
                    f.write(mol.GetProp("_Name") + '\t' + item[0] + '\t' + '\t'.join(map(str, ids)) + '\n')
            else:
                print('Molecule #%i was skipped due to error' % i)
                with open(error_fname, 'at') as f_err:
                    f_err.write('%s\t%s\tMolecule #%i\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                            os.path.basename(__file__), i))


def main():

    parser = argparse.ArgumentParser(description='Generate all possible fragments and find indices of them in '
                                                 'molecules from input sdf-file.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf file with standardized structures, molecules should have titles.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file; each line contains tab-separated fields: name of a molecule, '
                             'name of a fragment and list of corresponding atom numbers.')
    parser.add_argument('-q', '--query', metavar='smarts', default='[#6+0;!$(*=,#[!#6])]!@!=!#[*]',
                        help='SMARTS string to match bonds to cleave. Default: [#6+0;!$(*=,#[!#6])]!@!=!#[*]')
    parser.add_argument('-u', '--upper_cuts_number', metavar='integer', default=3,
                        help='maximal number of bonds cleaved simultaneously. Default: 3')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='show progress on the screen.')
    parser.add_argument('-e', '--error_file', metavar='log_file_name.txt', default="indigo_errors.txt",
                        help='save names of molecules which cause error to a text log file. Default file name '
                             'indigo_errors.txt.')


    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_sdf = v
        if o == "out": out_txt = v
        if o == "query": query = v
        if o == "upper_cuts_number": max_cuts = int(v)
        if o == "verbose": verbose = v
        if o == "error_file": error_fname = v

    main_params(in_sdf, out_txt, query, max_cuts, verbose, error_fname)


if __name__ == '__main__':
    main()
