#!/usr/bin/env python3

import argparse
from collections import OrderedDict
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs, Torsions

from sirms.files import SvmSaver, LoadFragments
from sirms.sirms import SaveSimplexes


mol_frag_sep = "###"


def CalcMolFP(m, i, opt_noH, f, frags=None, per_atom_fragments=None, id_field_name=None):

    def get_fp_as_dict(mol, opt_noH, f):
        """calc specified fingerprint for input molecule and return nonzero elements of it as dict
        """
        
        if opt_noH:
            # Chem.RemoveHs(mol)  # !it doesnt help, anyway next line gets them Hs back
            mol = Chem.RWMol(mol)
            for idx in reversed(range(mol.GetNumAtoms())):  # reverse because ids of atoms change while iter
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 1:
                    mol.RemoveAtom(idx)
            Chem.FastFindRings(mol)  # needs to calc morganfp, otherwise err "no ringinfo"
        if f in (GetMorganFingerprint_2, Pairs.GetAtomPairFingerprint, Torsions.GetTopologicalTorsionFingerprint,
                 Get_RDKFP_24):
            
            return {str(k): v for k, v in f(mol).GetNonzeroElements().items()}  # keys must be str input to saver, val-int
        elif f in [GetMorganFingerprint_2_bin, Pairs.GetAtomPairFingerprintAsBitVect, Get_RDKFP_24_bin]:
            return {str(k): int(v) for k, v in enumerate(DataStructs.BitVectToText(f(mol))) if v == '1'}

    mol_dict = OrderedDict()
    if id_field_name is not None:
        nm = m.GetProp(id_field_name)
    elif m.GetProp("_Name") == "":
        nm = 'auto_generated_id_' + str(i + 1)  # 1-based as in sirms.py
    else:
        nm = m.GetProp("_Name")
    res = get_fp_as_dict(m, opt_noH, f)
    if res:  # if all descriptors for a molecule are zero skip fragmentation (NC(=O)N gives no TT, removal of C(=O)N from it gives TT)
        mol_dict[nm] = res
        if per_atom_fragments:
            counter = 0
            for idx in range(m.GetNumAtoms()):
                if m.GetAtomWithIdx(idx).GetAtomicNum() > 1:
                    rw_m = Chem.RWMol(m)
                    rw_m.GetAtoms()[idx].SetAtomicNum(0)
                    mol_dict[nm + mol_frag_sep + str(idx + 1) + "#" + str(counter)] = get_fp_as_dict(rw_m, opt_noH, f)
                    counter += 1
        elif frags and nm in frags:
            for k, v in frags[nm].items():
                rw_m = Chem.RWMol(m)
                for idx in sorted(v, reverse=True):  # note we don't check if atom== H (is it ok?)
                    rw_m.GetAtoms()[idx-1].SetAtomicNum(0)
                mol_dict[nm + mol_frag_sep + k] = get_fp_as_dict(rw_m, opt_noH, f)
    return mol_dict


def main_params(in_fname, out_fname, output_format, get_fp, opt_verbose, opt_noH, frag_fname,
                per_atom_fragments, id_field_name):

    funcs = {'bMG2': GetMorganFingerprint_2_bin,
             'MG2': GetMorganFingerprint_2,
             'bAP': Pairs.GetAtomPairFingerprintAsBitVect,
             'AP': Pairs.GetAtomPairFingerprint,
             'TT': Torsions.GetTopologicalTorsionFingerprint,
             'bRDK': Get_RDKFP_24_bin,
             'RDK': Get_RDKFP_24}
    get_fp = funcs[get_fp]

    # load sdf and get dict of fp (like sirms dict)
    input_file_extension = in_fname.strip().split(".")[-1].lower()
    if input_file_extension == 'sdf':
        saver = None
        mols = None
        if output_format == "svm":
            saver = SvmSaver(out_fname)
        if output_format == "txt":
            mols = OrderedDict()  # key - molname, val- mol; if frags: key - molname or mol+fragname, val-mol for mol or part b
        frags = LoadFragments(frag_fname)
        for i, m in enumerate(Chem.SDMolSupplier(in_fname, removeHs=False)):
            if m is not None:
                res = CalcMolFP(m, i, opt_noH=opt_noH, f=get_fp, frags=frags, per_atom_fragments=per_atom_fragments,
                                id_field_name = id_field_name)
                if output_format == "txt":
                    mols.update(res)
                if output_format == "svm":
                    for mol_name, descr_dict in res.items():
                        saver.save_mol_descriptors(mol_name, descr_dict)
        if output_format == "txt":
            SaveSimplexes(out_fname, mols, output_format)
    else:
        print("Input file extension should be SDF Current file has %s. Please check it." %
              input_file_extension.upper())
        return None


def GetMorganFingerprint_2(m):
    return AllChem.GetMorganFingerprint(m, radius=2)


def GetMorganFingerprint_2_bin(m):
    return AllChem.GetMorganFingerprintAsBitVect(m, radius=2)


def Get_RDKFP_24(m):
    return AllChem.UnfoldedRDKFingerprintCountBased(m, minPath=2, maxPath=4)


def Get_RDKFP_24_bin(m):
    return AllChem.RDKFingerprint(m, minPath=2, maxPath=4)


def entry_point():
    parser = argparse.ArgumentParser(description='Calculate fingerprint descriptors')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=False, default=None,
                        help='input file ( sdf with standardized structures')
    parser.add_argument('-o', '--out', metavar='output.txt', required=False, default=None,
                        help='output file with calculated descriptors. Can be in text or sparse svm format.')
    parser.add_argument('-b', '--output_format', metavar='output_format', default='txt',
                        help='format of output file with calculated descriptors (txt|svm). '
                             'Txt is ordinary tab-separated text file. Svm is sparse format, two additional files will '
                             'be saved with extensions .colnames and .rownames. Default: txt.')
    parser.add_argument('--fp_type', required=True,
                        help='Type of fingerprints to calculate choose one of: MG2 or bMG2 (Morgan radius 2), '
                             'AP or bAP (atom-pair), RDK or bRDK (2-4 atoms RDK fingerprint), '
                             'TT (topological torsion). Prefix b means binary fingerprint of length 2048.')
    parser.add_argument('-x', '--noH', action='store_true', default=False,
                        help='if set this flag hydrogen atoms will be excluded from the simplexes calculation.')
    parser.add_argument('-f', '--fragments', metavar='fragments.txt', default=None,
                        help='text file containing list of names of single compounds, fragment names and atom '
                             'indexes of fragment to remove (all values are tab-separated).')
    parser.add_argument('--per_atom_fragments', action='store_true', default=False,
                        help='if set this flag input fragments will be omitted and single atoms will be considered '
                             'as fragments.')
    parser.add_argument('-w', '--id_field_name', metavar='field_name', default=None,
                        help='field name of unique ID for compounds. '
                             'If omitted for sdf molecule titles will be used or auto-generated names')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='if set this flag progress will be printed out (may cause decrease in speed).')

    args = vars(parser.parse_args())

    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "output_format": output_format = v
        if o == "fp_type": get_fp = v
        if o == "verbose": opt_verbose = v
        if o == "noH": opt_noH = v
        if o == "fragments": frag_fname = v
        if o == "per_atom_fragments": per_atom_fragments = v
        if o == "id_field_name": id_field_name = v

    main_params(in_fname=in_fname, out_fname=out_fname, output_format=output_format, get_fp=get_fp,
                opt_verbose=opt_verbose, opt_noH=opt_noH, frag_fname=frag_fname,
                per_atom_fragments=per_atom_fragments, id_field_name=id_field_name)


if __name__ == '__main__':
    entry_point()
