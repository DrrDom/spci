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
import argparse
from subprocess import call


def quote_str(s):
    return '"%s"' % s


def add_pH(params, pH):
    if pH is not None:
        params.extend(["-H", pH])
    return (params)


def add_HB(molstr):

    for i, line in enumerate(molstr):
        if line.strip() == ">  <ACC>":
            acc = ["I" if el == "0" else "A" for el in molstr[i+1].strip().split(";")]
        if line.strip() == ">  <DON>":
            don = ["I" if el == "0" else "D" for el in molstr[i+1].strip().split(";")]

    hb = []
    for a, d in zip(acc, don):
        if a == "A" and d == "D":
            hb.append("A|D")
        elif a == "A":
            hb.append("A")
        elif d == "D":
            hb.append("D")
        else:
            hb.append("I")

    molstr.insert(len(molstr) - 1, ">  <HB>\n")
    molstr.insert(len(molstr) - 1, ";".join(hb) + "\n")
    molstr.insert(len(molstr) - 1, "\n")

    return molstr


def main_params(in_fname, out_fname, prop, pH, cxcalc_path):

    start_params = [cxcalc_path]   # remove quotes because Win cannot run "cxcalc" only "dir/cxcalc"
    start_params.extend(["-o", quote_str(out_fname)])
    start_params.append("-S")
    log_fname = os.path.splitext(out_fname)[0] + "_cxcalc.log"
    start_params.extend(["--log", quote_str(log_fname)])
    start_params.append(quote_str(in_fname))

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
        call(' '.join(run_params), shell=True)

    # read numbers of molecules which produce errors
    mol_errors = []
    if os.path.isfile(log_fname):
        with open(log_fname) as f:
            for line in f:
                if line.find("<_MOLCOUNT>") >= 0:
                    mol_errors.append(int(f.readline().strip()))

    # create new field HB - hydrogen bond parameters
    # A - acceptor
    # D - donor
    # AD - acceptor & donor
    # I - neither acceptor nor donor
    # and remove error molecules

    tmp_out_fname = out_fname + ".tmp"
    add_hb = "acc" in prop and "don" in prop
    if add_hb or mol_errors:
        with open(tmp_out_fname, "wt") as fout:
            with open(out_fname) as fin:

                molstr = []
                cur_mol_id = 1
                for line in fin:
                    if line.rstrip() != "$$$$":
                        molstr.append(line)
                    else:
                        molstr.append(line)
                        if mol_errors and cur_mol_id in mol_errors:
                            molstr = []
                        elif add_hb:
                            molstr = add_HB(molstr)
                            fout.writelines(molstr)
                            molstr = []
                        cur_mol_id += 1

        os.remove(out_fname)
        os.rename(tmp_out_fname, out_fname)


def main():

    parser = argparse.ArgumentParser(description='Calculate atomic properties with Chemaxon cxcalc tool. '
                                                 'Molecules which produce errors are omitted from output file. '
                                                 'Omitted molecules will be listed in the created log file.')
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
                        default="cxcalc",
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
