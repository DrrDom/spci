#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2014
# license         : GPL3
#==============================================================================

import sys
import argparse


def main_params(in_fname, out_fname, title, field_names):

    with open(in_fname) as ifs:
        with open(out_fname, "w") as ofs:

            header = "\t".join(field_names) + "\n"
            if title:
                header = "Title\t" + header

            ofs.write(header)
            sys.stdout.write(header)

            output_data = [ifs.readline().strip()] if title else []

            for line in ifs:
                if title:
                    start_pos = 1
                else:
                    start_pos = 0
                if line[0] == ">":
                    for field_name in field_names:
                        if line.find("<" + field_name + ">") >= 0:
                            output_data.append(ifs.readline().strip())
                            break
                elif line.find("$$$$") > -1:
                    ofs.write("\t".join(output_data) + "\n")
                    sys.stdout.write("\t".join(output_data) + "\n")
                    output_data = [ifs.readline().strip()] if title else []


def main():
    parser = argparse.ArgumentParser(description='Extract field values from sdf-files.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf file with standardized structures, molecules should have titles')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file with values of specified fields')
    parser.add_argument('-t', '--title', action='store_true', default=False,
                        help='If true molecules titles will be extracted')
    parser.add_argument('-f', '--field_names', metavar='[field_name_1 field_name_2 ...]', required=True, nargs='*',
                        help='space separated list of field names for extraction')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "title": title = v
        if o == "field_names": field_names = v

    main_params(in_fname, out_fname, title, field_names)


if __name__ == '__main__':
    main()
