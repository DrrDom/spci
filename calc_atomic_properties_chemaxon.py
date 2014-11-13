#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3
# copyright       : Pavel Polishchuk 2014
# license         : GPL3
#==============================================================================

import os
import argparse
from subprocess import call


def add_pH(params, pH):
    if pH is not None:
        params.extend(["-H", pH])
    return (params)


def main_params(in_fname, out_fname, prop, pH, cxcalc_path):

    start_params = [cxcalc_path]
    start_params.extend(["-o", out_fname])
    start_params.append("-S")
    start_params.append(in_fname)

    # logp, charge, atom_polarizability, refractivity, acc, don
    run_params = start_params[:]
    if "logp" in prop:
        run_params.extend(["logp", "-t", "increments"])
        run_params = add_pH(run_params, pH)
    if "charge" in prop:
        run_params.append("charge")
        run_params = add_pH(run_params, pH)
    if "atompol" in prop:
        run_params.append("atompol")
        run_params = add_pH(run_params, pH)
    if "refractivity" in prop:
        run_params.extend(["refractivity", "-t", "increments"])
    if "acc" in prop:
        run_params.append("acc")
        run_params = add_pH(run_params, pH)
    if "don" in prop:
        run_params.append("don")
        run_params = add_pH(run_params, pH)
    if run_params != start_params:
        call(run_params, shell=True)

    # create new field HB - hydrogen bond parameters
    # A - acceptor
    # D - donor
    # AD - acceptor & donor
    # I - neither acceptor nor donor
    tmp_out_fname = out_fname + ".tmp"
    if "acc" in prop and "don" in prop:
        with open(tmp_out_fname, "wt") as fout:
            with open(out_fname) as fin:
                line = fin.readline()
                while line:
                    if line.strip() != "$$$$":
                        fout.write(line)
                    if line.strip() == ">  <ACC>":
                        tmp = fin.readline()
                        fout.write(tmp)
                        acc = ["I" if el == "0" else "A" for el in tmp.strip().split(";")]
                    if line.strip() == ">  <DON>":
                        tmp = fin.readline()
                        fout.write(tmp)
                        don = ["I" if el == "0" else "D" for el in tmp.strip().split(";")]
                    if line.strip() == "$$$$":
                        hb = []
                        for a, d in zip(acc, don):
                            if a == "A" and d == "D":
                                hb.append("AD")
                            elif a == "A":
                                hb.append("A")
                            elif d == "D":
                                hb.append("D")
                            else:
                                hb.append("I")
                        fout.write(">  <HB>\n")
                        fout.write(";".join(hb) + "\n")
                        fout.write("\n")
                        fout.write(line)
                    line = fin.readline()
        os.remove(out_fname)
        os.rename(tmp_out_fname, out_fname)


def main():

    parser = argparse.ArgumentParser(description='Calculate atomic properties with Chemaxon cxcalc tool.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf file with standardized structures, molecules should have titles')
    parser.add_argument('-o', '--out', metavar='output.sdf', required=True,
                        help='output sdf file with added calculated atomic properties')
    parser.add_argument('-p', '--properties', metavar='', nargs='*',
                        default=['charge', 'logp', 'acc', 'don', 'refractivity'],
                        help="list of atomic properties to calc. Default: 'charge', 'logp', 'acc', 'don', "
                             "'refractivity'. If 'acc' and 'don' specified simultaneously hydrogen bonding "
                             "descriptors (hb) will be generated.")
    parser.add_argument('-H', '--pH', metavar='pH_value', default=None,
                        help='pH value which will be used during calculation of atomic properties')
    parser.add_argument('-c', '--cxcalc_path', metavar='path-to-cxcalc.bat',
                        default="c:\\Program Files (x86)\\ChemAxon\\JChem\\bin\\cxcalc",
                        help='path to cxcalc executable file. It is usually '
                             'c:\\Program Files (x86)\\ChemAxon\\JChem\\bin\\cxcalc on Windows')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "properties": prop = v
        if o == "pH": pH = v
        if o == "cxcalc_path": cxcalc_path = v

    main_params(in_fname, out_fname, prop, pH, cxcalc_path)


if __name__ == '__main__':
    main()
