#!/usr/bin/env python3
# author          : Pavel
# date            : 04.12.15
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2015
# license         : GPL3
#==============================================================================

import os
import re
import argparse

from itertools import combinations, permutations
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import AllChem

from mol_context import get_canon_context_core, get_submol


patt = re.compile("\[\*\:[0-9]+\]")  # to change CC([*:1])O to CC([*])O


def replace_no2(mol):
    query = Chem.MolFromSmarts('n(:o):o')
    repl = Chem.MolFromSmiles('[N+](=O)[O-]')
    return AllChem.ReplaceSubstructs(mol, query, repl, replaceAll=True)[0]


def frag_mol_by_cuts(mol, cut_list, keep_stereo, radius):

    em = Chem.EditableMol(mol)

    output = []
    ap_ids = []

    for i, cut in enumerate(cut_list):   # this can be replaced with FragmentOnBonds
        em.RemoveBond(cut[0], cut[1])
        ap1 = em.AddAtom(Chem.Atom(0))
        ap_ids.append(ap1)
        em.AddBond(cut[0], ap1, Chem.BondType.SINGLE)
        ap2 = em.AddAtom(Chem.Atom(0))
        ap_ids.append(ap2)
        em.AddBond(cut[1], ap2, Chem.BondType.SINGLE)

    mol = em.GetMol()

    # label cut points
    for i, ids in enumerate(zip(ap_ids[0::2], ap_ids[1::2])):
        for id in ids:
            mol.GetAtomWithIdx(id).SetAtomMapNum(i + 1)

    # if split produce more than one fragment with several cuts it is illegible
    # ex: C* *CO* *CC* *O
    att_num = []
    frags = Chem.GetMolFrags(mol, asMols=False)
    for ids in frags:
        att_num.append(sum(mol.GetAtomWithIdx(i).GetAtomicNum() == 0 for i in ids))
    if sum(i > 1 for i in att_num) > 1:
        return output

    # remove Hs from ids
    frags = tuple(tuple(i for i in ids if mol.GetAtomWithIdx(i).GetAtomicNum() != 1) for ids in frags)

    # one cut
    if len(frags) == 2:
        for core_ids, context_ids in permutations(frags, 2):
            core = get_submol(mol, core_ids)
            if radius > 0:
                context = get_submol(mol, context_ids)
                context, core = get_canon_context_core(context, core, radius, keep_stereo)
                output.append((core + '|' + context, tuple(sorted(i for i in core_ids if i not in ap_ids))))
            else:
                # remove all atom map numbers to obtain SMILES compatible with SMILES having atom map numbers
                context, core = get_canon_context_core('', core, radius, keep_stereo)
                output.append((core, tuple(sorted(i for i in core_ids if i not in ap_ids))))
    # two amd more cuts
    else:
        core_ids = []
        context_ids = []
        for n, f in zip(att_num, frags):
            if n == 1:
                context_ids.extend(f)
            else:
                core_ids.extend(f)
        core = get_submol(mol, core_ids)
        if radius > 0:
            context = get_submol(mol, context_ids)
            context, core = get_canon_context_core(context, core, keep_stereo, radius)
            output.append((core + '|' + context, tuple(sorted(i for i in core_ids if i not in ap_ids))))
        else:
                context, core = get_canon_context_core('', core, radius, keep_stereo)
                output.append((core, tuple(sorted(i for i in core_ids if i not in ap_ids))))

    return output


def fragment_mol(mol, query, max_cuts, keep_stereo, radius):
    # returns list of lists: [['F', [0]], ['C#N', [3, 4]], ... ]

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
            output.extend(frag_mol_by_cuts(mol, comb, keep_stereo, radius))

    return output


def main_params(in_sdf, out_txt, query, max_cuts, keep_stereo, radius, verbose, error_fname):

    query = Chem.MolFromSmarts(query)

    with open(out_txt, 'wt') as f:
        for i, mol in enumerate(Chem.SDMolSupplier(in_sdf, sanitize=False, removeHs=False)):
            if mol is not None:
                if verbose:
                    print(mol.GetProp("_Name") + ' is processing')
                res = fragment_mol(mol, query, max_cuts, keep_stereo, radius)
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
    parser.add_argument('-q', '--query', metavar='smarts', default='[#6+0;!$(*=,#[!#6])]!@!=!#[!#1]',
                        help='SMARTS string to match bonds to cleave. Default: [#6+0;!$(*=,#[!#6])]!@!=!#[!#1]')
    parser.add_argument('-u', '--upper_cuts_number', metavar='integer', default=3,
                        help='maximal number of bonds cleaved simultaneously. Default: 3')
    parser.add_argument('-r', '--radius', metavar='integer', default=0,
                        help='radius of molecular context (in bonds) which will be taken into account. '
                             '0 means no context. Default: 0.')
    parser.add_argument('-s', '--keep_stereo', action='store_true', default=False,
                        help='set this flag to keep stereo in context and core parts.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='show progress on the screen.')
    parser.add_argument('-e', '--error_file', metavar='log_file_name.txt', default="rdkit_errors.txt",
                        help='save a number of molecules which cause error. Default file name: rdkit_errors.txt.')


    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_sdf = v
        if o == "out": out_txt = v
        if o == "query": query = v
        if o == "upper_cuts_number": max_cuts = int(v)
        if o == "verbose": verbose = v
        if o == "error_file": error_fname = v
        if o == "keep_stereo": keep_stereo = v
        if o == "radius": radius = int(v)

    main_params(in_sdf, out_txt, query, max_cuts, keep_stereo, radius, verbose, error_fname)


if __name__ == '__main__':
    main()
