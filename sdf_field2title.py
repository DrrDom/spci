#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3
# copyright       : Pavel Polishchuk 2014
# license         : GPL3
#==============================================================================

import argparse


def extract_field_value(mol, field_name):
    for i, line in enumerate(mol):
        if line.strip().find("> <" + field_name + ">") == 0 or line.strip().find(">  <" + field_name + ">") == 0:
            return mol[i+1].strip()


def set_mol_name(mol, mol_name):
    mol[0] = mol_name + "\n"
    return mol


def main_params(input_sdf_fname, field_name, output_sdf_fname):

    with open(input_sdf_fname, "rt") as fin:
        with open(output_sdf_fname, "wt") as fout:
            mol = []
            id = 1
            for line in fin:
                if line.strip() != "$$$$":
                    mol.append(line)
                else:
                    mol.append(line)
                    if field_name is None:
                        tmp = 'MolID_' + str(id)
                        id += 1
                    else:
                        tmp = extract_field_value(mol, field_name)
                    print(tmp)
                    set_mol_name(mol, tmp)
                    fout.writelines(mol)
                    mol = []


def main():

    parser = argparse.ArgumentParser(description='Plot fragments contributions.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input sdf file.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output sdf file.')
    parser.add_argument('-f', '--field_name', default=None,
                        help='field name in sdf file which will be copied to molecule title, '
                             'if omitted sequence number will be inserted.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fname = v
        if o == "output": output_fname = v
        if o == "field_name": field_name = v

    main_params(input_fname, field_name, output_fname)


if __name__ == '__main__':
    main()