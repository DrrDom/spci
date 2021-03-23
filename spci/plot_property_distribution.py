#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2014
# license         : LGPLv3
#==============================================================================

import argparse
from matplotlib import pyplot as plt
from collections import Counter


def readDataField(sdf_fname, field_names):
    d = dict()
    field_names = [">  <" + el.lower() + ">" for el in field_names]
    with open(sdf_fname) as f:
        line = f.readline()
        while line:
            if line.strip().lower() in field_names:
                field_name = line.strip().lower()[4:-1]
                if field_name in d.keys():
                    d[field_name].extend(f.readline().strip().split(";"))
                else:
                    d[field_name] = f.readline().strip().split(";")
            line = f.readline()

    for field_name in d.keys():
        try:
            d[field_name] = [float(el.replace(",", ".")) for el in d[field_name] if el != ""]
        except ValueError:
            d[field_name] = [el for el in d[field_name] if el != ""]

    return (d)


def plot(data):
    def autolabel(splot, rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            splot.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')

    fig = plt.figure()
    h = 1  # number of rows on plot
    w = 1  # number of columns on plot
    for k, i in data.items():
        if i:
            # change w and h values
            if len(fig.axes) >= w * h:
                if w == h:
                    w += 1
                else:
                    h += 1
            # resize previous;y added subplots
            for j in range(len(fig.axes)):
                fig.axes[j].change_geometry(h, w, j + 1)
            # add new subplot
            splot = fig.add_subplot(h, w, len(fig.axes) + 1)
            splot.set_xlabel(k)
            if type(i[0]) is float:
                n, bins, patches = splot.hist(i, 50, facecolor='green', alpha=0.5, histtype="bar")
            else:
                c = Counter(i)
                rect = splot.bar(range(len(c)), list(c.values()), width=0.8)
                splot.set_xticks([el + 0.4 for el in range(len(c))])
                splot.set_xticklabels(list(c.keys()))
                autolabel(splot, rect)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='Plot distribution diagrams for data fields in sdf-file.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf file with standardized structures, molecules should have titles.')
    parser.add_argument('-p', '--properties', metavar='', nargs='*',
                        default=['charge', 'logp', 'acc', 'don', 'atom_polarizability', 'refractivity', 'hb'],
                        help="list of properties (field names) to plot distribution diagrams. Default: 'charge', 'logp', 'acc', 'don', 'atom_polarizability', 'refractivity', 'hb'.")

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "properties": prop = v

    d = readDataField(in_fname, prop)

    plot(d)


if __name__ == '__main__':
    main()
