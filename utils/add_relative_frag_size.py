#==============================================================================
# author          : Maria Matveiva
# date            : 07-06-2018
# version         : 0.1
# python_version  : 3.5.2
# copyright       : Maria Matveieva 2018
# license         : GPL3
#==============================================================================

from rdkit import Chem
import argparse
import re
import shutil


def main_params(contrib_fname, inp_format, out_fname, sd_fname):
    mol_dict = load_mol_dict(sd_fname)
    if inp_format == "wide":
        fr_dict, fr_lst = load_frags_wide(contrib_fname)
        write_wide(contrib_fname, fr_dict, fr_lst, mol_dict, out_fname)
    else:
        fr_dict = load_frags_long(contrib_fname)
        write_long(contrib_fname, out_fname, fr_dict, mol_dict)


def load_mol_dict(sd_fname):
    sdf = Chem.SDMolSupplier(sd_fname)
    mols_dict = {}
    for i in sdf:
        if i is not None and i.GetProp("_Name") not in mols_dict:  
            mols_dict[i.GetProp("_Name")] = i.GetNumHeavyAtoms()
    return mols_dict


def load_frags_wide(contrib_fname):
    fr_dict = {}
    fr_lst = []
    with open(contrib_fname, "rt") as inp:
        for mol_fr_c in inp.readline().strip().split('\t')[1:]:
            tmp = re.split("###|\|\||#\\d+$", mol_fr_c)
            fr_lst.append(tmp[:2])
            if tmp[1] not in fr_dict:
                fr_dict[tmp[1]] = Chem.MolFromSmarts(tmp[1]).GetNumHeavyAtoms()
    return fr_dict, fr_lst


def load_frags_long(contrib_fname):
    fr_dict = {}
    with open(contrib_fname, "rt") as inp:
        inp.readline()
        for line in inp:
            tmp = line.strip().split("\t")[0:2]
            if tmp[1] not in fr_dict:
                fr_dict[tmp[1]] = Chem.MolFromSmarts(tmp[1]).GetNumHeavyAtoms()
    return fr_dict


def write_wide(contrib_fname, fr_dict, fr_lst, mol_dict, out_fname):

    m = ["mol_size"]
    f = ["frag_size"]
    p = ["relative_frag_size"]

    for mol, frag in fr_lst:
        m.append(mol_dict[mol])
        f.append(fr_dict[frag])
        p.append(round(fr_dict[frag] / mol_dict[mol], 2))
    if out_fname is not None:
        shutil.copyfile(contrib_fname, out_fname)
        with open(out_fname, "a") as out:
            out.write("\t".join(map(str, m)) + "\n")
            out.write("\t".join(map(str, f)) + "\n")
            out.write("\t".join(map(str, p)) + "\n")
    else:
        with open(contrib_fname, "a") as contr_file:
            contr_file.write("\t".join(map(str, m)) + "\n")
            contr_file.write("\t".join(map(str, f)) + "\n")
            contr_file.write("\t".join(map(str, p)) + "\n")


def write_long(contrib_fname, out_fname, fr_dict, mol_dict):
    with open(out_fname, "wt") as out:
        with open(contrib_fname, "rt") as contr_file:
            out.write(contr_file.readline().strip() + "\t" + "\t".join(["mol_size", "frag_size", "relative_frag_size"]) + "\n")
            for line in contr_file:
                tmp = line.strip().split()
                out.write(line.strip() + "\t" + "\t".join(map(str, [mol_dict[tmp[0]],
                                                                    fr_dict[tmp[1]],
                                                                    round(fr_dict[tmp[1]]/mol_dict[tmp[0]], 2)]))
                          + "\n")


def main():

    parser = argparse.ArgumentParser(description='Calculate size (heavy atom count, HAC) and relative size '
                                                 '(HAC of fragment/HAC of molecule) for fragments in molecules '
                                                 '(HAC of molecule will also be calculated and written)')
    parser.add_argument('-i', '--input', metavar='input_contributions.txt', required=True,
                        help='input text file with fragments contributions')
    parser.add_argument('-f', '--input_format', metavar='wide|long', default='wide',
                        help='format of input file wide/long. Default: wide.')
    parser.add_argument('-o', '--output', metavar='output_contributions.txt', default=None,
                        help='output text file will contain the same data as the input file with appended cols '
                             '(long format) or rows (wide format) with calculated data. '
                             'This argument is required only for long format. '
                             'If it is omitted for wide format new data will be appended to the input file '
                             '(the file will be overwritten).')
    parser.add_argument('--sdf', metavar='input.sdf', required=True,
                        help='input sdf with molecules from which frags were derived. Should contain molecule titles '
                             'identical to those used in the input file with contributions.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": contrib_fname = v
        if o == "input_format": inp_format = v
        if o == "output": out_fname = v
        if o == "sdf": sd_fname = v
       
    if inp_format not in ['wide', 'long']:
        print("Wrong input file format - %s. Only wide and long file formats are allowed." % inp_format)
        exit()
    if inp_format == 'long' and not out_fname:
        print("Error. -o argument is missing. For long format -o --out argument is required")
        exit()

    main_params(contrib_fname, inp_format, out_fname, sd_fname)


if __name__ == '__main__':
    main()
