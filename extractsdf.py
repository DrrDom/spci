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
from itertools import chain


def main_params(in_fname, out_fname, title, field_names, all_fields):

    with open(in_fname) as ifs:

        output = []

        # add title
        if title:
            title_str = ifs.readline().strip()
            if title_str != '':
                output.append({'Title': title_str})
            else:
                output.append({'Title': 'NA'})
        else:
            output.append(dict())

        for line in ifs:

            if line[0] == ">":

                if all_fields:
                    start = line.find("<") + 1
                    field_name = line.strip()[start:-1]
                    output[-1][field_name] = ifs.readline().strip()

                elif field_names is not None:
                    for field_name in field_names:
                        if line.find("<" + field_name + ">") >= 0:
                            output[-1][field_name] = ifs.readline().strip()
                            break

            elif line.find("$$$$\n") > -1:
                # add title
                if title:
                    title_str = ifs.readline().strip()
                    if title_str != '':
                        output.append({'Title': title_str})
                    else:
                        output.append({'Title': 'NA'})
                else:
                    output.append(dict())

        # remove last blank record from the output list
        output.pop()

    with open(out_fname, "w") as ofs:
        # get sorted unique field names
        field_names = sorted(list(set(list(chain.from_iterable([list(s.keys()) for s in output])))))
        if 'Title' in field_names:
            field_names.remove('Title')
            field_names.insert(0, 'Title')
        ofs.write("\t".join(field_names) + "\n")
        sys.stdout.write("\t".join(field_names) + "\n")
        for item in output:
            line = [item.get(f, 'NA') for f in field_names]
            ofs.write("\t".join(line) + "\n")
            sys.stdout.write("\t".join(line) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Extract field values from sdf-files.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf file with standardized structures, molecules should have titles')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file with values of specified fields')
    parser.add_argument('-t', '--title', action='store_true', default=False,
                        help='If true molecules titles will be extracted')
    parser.add_argument('-f', '--field_names', metavar='[field_name_1 field_name_2 ...]',
                        required=False, default=None, nargs='*',
                        help='space separated list of field names for extraction')
    parser.add_argument('-a', '--all_fields', action='store_true', default=False,
                        help='If set (true) all fields will be extracted.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "title": title = v
        if o == "field_names": field_names = v
        if o == "all_fields": all_fields = v

    main_params(in_fname, out_fname, title, field_names, all_fields)


if __name__ == '__main__':
    main()
