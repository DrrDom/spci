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
from indigo import Indigo, IndigoException


indigo = Indigo()


def get_bond_indigo(mol, atom1_id, atom2_id):
    for nei in mol.getAtom(atom1_id).iterateNeighbors():
        if nei.index() == atom2_id:
            return nei.bond()
    return None


def replace_no2(mol):
    query = indigo.loadSmarts('n(:o):o')
    matcher = indigo.substructureMatcher(mol)
    for match in matcher.iterateMatches(query):
        ids = []
        for atom in query.iterateAtoms():
            if match.mapAtom(atom) is not None:
                ids.append(match.mapAtom(atom).index())
        for i in ids:
            if mol.getAtom(i).atomicNumber() == 7:
                ids.remove(i)
                mol.getAtom(i).setCharge(1)
                mol.getAtom(ids[0]).setCharge(-1)
                get_bond_indigo(mol, i, ids[0]).setBondOrder(1)
                get_bond_indigo(mol, i, ids[1]).setBondOrder(2)


def fix_cansmi_attach_point(smi):
    # fix cansmi name by replacement of attachment points: Cl%91.[*:1]%91 Cl[*:1]
    p = smi.split(".")
    output = p[0]
    for i in p[1:]:
        att = i.split('%')
        output = output.replace('%' + att[1], "(" + att[0] + ")")
    return output


def frag_mol_by_cuts(mol, cut_list):

    output = []

    for cut in cut_list:
        a1 = mol.getAtom(cut[0])
        for nei in a1.iterateNeighbors():
            if nei.index() == cut[1]:
                nei.bond().remove()
                break
        # add attachment points
        mol.getAtom(cut[0]).setAttachmentPoint(1)
        mol.getAtom(cut[1]).setAttachmentPoint(1)

    for comp in mol.iterateComponents():
        cansmi = comp.clone().canonicalSmiles().split(" ")[0]
        cansmi = fix_cansmi_attach_point(cansmi)
        ids = [a.index() for a in comp.iterateAtoms()]
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

    # modify representation of NO2 groups to charged version
    replace_no2(mol)

    output = []

    matcher = indigo.substructureMatcher(mol)
    all_cuts = []
    for match in matcher.iterateMatches(query):
        ids = []
        for atom in query.iterateAtoms():
            if match.mapAtom(atom) is not None:
                ids.append(match.mapAtom(atom).index())
        all_cuts.append(ids)

    for i in range(1, max_cuts + 1):

        for comb in combinations(all_cuts, i):
            output.extend(frag_mol_by_cuts(mol.clone(), comb))

        # for cut in all_cuts:
        #     output.extend(frag_mol_by_cuts(mol.clone(), [cut]))
        #
        # for comb in combinations(all_cuts, 2):
        #     output.extend(frag_mol_by_cuts(mol.clone(), comb))
        #
        # for comb in combinations(all_cuts, 3):
        #     output.extend(frag_mol_by_cuts(mol.clone(), comb))

    output = filter_dupl(output)

    return output


def main_params(in_sdf, out_txt, query, max_cuts, verbose, error_fname):

    query = indigo.loadSmarts(query)

    with open(out_txt, 'wt') as f:
        for mol in indigo.iterateSDFile(in_sdf):
            try:
                if verbose:
                    print(mol.name() + ' is processing')
                res = fragment_mol(mol, query, max_cuts)
                for item in res:
                    ids = [i + 1 for i in item[1]]
                    f.write(mol.name() + '\t' + item[0] + '\t' + '\t'.join(map(str, ids)) + '\n')
            except IndigoException as e:
                print('%s was skipped due to error' % mol.name())
                print(e)
                with open(error_fname, 'at') as f_err:
                    f_err.write('%s\t%s\t%s\t%s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                      os.path.basename(__file__), mol.name(), e))


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
        if o =="upper_cuts_number": max_cuts = v
        if o == "verbose": verbose = v
        if o == "error_file": error_fname = v

    main_params(in_sdf, out_txt, query, max_cuts, verbose, error_fname)


if __name__ == '__main__':
    main()
