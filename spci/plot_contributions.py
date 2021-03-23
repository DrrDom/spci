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
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import numpy as np

mol_frag_sep = "###"


def read_contrib_file(fname, model_names, min_M, min_N, contr_names):
    """
    OUTPUT
    {'gbm': {'charge': [0.119, 0.174, 0.083, 0.220, ... ],
            {'overall': [0.201, 0.256, 0.283, 0.420, ... ],
              ...
            },
     'rf': ...}
    """
    d = defaultdict(dict)

    with open(fname) as f:

        names = f.readline().strip().split('\t')[1:]

        # count number of molecules for each fragment
        frag_mol_count = defaultdict(int)
        frag_names = []

        for n in names:
            mol_name, frag_name = n.split(mol_frag_sep)
            frag_mol_count[frag_name] += 1
            frag_names.append(frag_name)

        # count number of each fragment
        frag_count = Counter(frag_names)

        # create new fragment names and create list of filtered indices (according to min_M and min_N)
        keep_ids = []
        for i, v in enumerate(frag_names):
            frag_names[i] = frag_names[i] + " (M=" + str(frag_mol_count[v]) + ", N=" + str(frag_count[v]) + ")"
            if frag_mol_count[v] >= min_M and frag_count[v] >= min_N:
                keep_ids.append(i)

        # filter out frag_names by keep_ids
        frag_names = [frag_names[i] for i in keep_ids]

        for line in f:
            tmp = line.strip().split('\t')
            model_name, prop_name = tmp[0].rsplit("_", 1)
            # skip contributions which are not selected
            if prop_name not in contr_names:
                continue
            if "all" in model_names or model_name in model_names:
                values = list(map(float, tmp[1:]))
                # filter out values by keep_ids
                values = [values[i] for i in keep_ids]
                d[model_name][prop_name] = values

    return frag_names, d


def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length // 2] + sorts[length // 2 - 1]) / 2.0
    return sorts[length // 2]


def add_consensus_average(contr_dict):
    consensus_dict = {}
    for model_name in contr_dict.keys():
        for prop_name in contr_dict[model_name].keys():
            if prop_name not in consensus_dict.keys():
                consensus_dict[prop_name] = [0] * len(contr_dict[model_name][prop_name])
            consensus_dict[prop_name] = [a + b for (a,b) in zip(consensus_dict[prop_name],
                                                                contr_dict[model_name][prop_name])]
    for prop_name in consensus_dict:
        consensus_dict[prop_name] = [v / len(contr_dict) for v in consensus_dict[prop_name]]
    contr_dict["consensus"] = consensus_dict


def prep_data_boxplot(contr_dict, frag_names_list, add_consensus=True):
    """
    INPUT
    {'gbm': {'charge': [0.119, 0.174, 0.083, 0.220, ... ],
             'overall': [0.201, 0.256, 0.283, 0.420, ... ],
              ...
            },
     'rf': ...}

    OUTPUT
    fragment_names_list = ['SO2', 'NO2', ...]  - corresponds to number of rows
    {'gbm': {'charge': [[0.119, 0.174, 0.083, 0.220, ... ],
                        [0.119, 0.174, 0.083, 0.220, ... ],
                        ...
                       ],
            {'overall': [[0.201, 0.256, 0.283, 0.420, ... ],
                         [0.429, 0.345, 0.121, 0.301, ... ],
                         ...
                        ],
     'rf': ...}
    """
    if add_consensus and len(contr_dict) > 1:
        add_consensus_average(contr_dict)
    output_dict = dict()
    for model_name in contr_dict.keys():
        output_dict[model_name] = dict()
        for prop_name in contr_dict[model_name]:
            tmp = defaultdict(list)
            for i, value in enumerate(contr_dict[model_name][prop_name]):
                tmp[frag_names_list[i]].append(value)
            sorted_frag_names = sorted(tmp.keys())
            output_dict[model_name][prop_name] = [tmp[n] for n in sorted_frag_names]
    return sorted_frag_names, output_dict


def sorted_ids(ls, ascending=True):
    if ascending:
        ids = sorted(range(len(ls)), key=lambda k: ls[k])
    else:
        ids = sorted(range(len(ls)), key=lambda k: -ls[k])
    ids = sorted(enumerate(ids), key=lambda x: x[1])
    ids = [i for i, j in ids]
    return ids


def sort_by_median(list_of_lists, ascending=True):
    m = [median(i) for i in list_of_lists]
    ids = sorted_ids(m, ascending)
    return ids


def prep_data_barplot(contr_dict, frag_names_list, add_consensus=True):
    """
    INPUT
    {'gbm': {'charge': [0.119, 0.174, 0.083, 0.220, ... ],
             'overall': [0.201, 0.256, 0.283, 0.420, ... ],
              ...
            },
     'rf': ...}

    OUTPUT
    fragment_names_list = ['SO2', 'NO2', ...]  - corresponds to number of rows
    {'gbm': {'charge': [0.119, 0.219, ... ],   - median contribution values
            {'overall': [0.119, 0.219, ... ],
                         ...
                        ],
     'rf': ...}
    """
    sorted_frag_names, d = prep_data_boxplot(contr_dict, frag_names_list, add_consensus)
    for model_name in d.keys():
        for prop_name in d[model_name].keys():
            d[model_name][prop_name] = [median(v) for v in d[model_name][prop_name]]
    return sorted_frag_names, d


def prep_data_barplot_binary_class(contr_dict, frag_names_list, add_consensus=True):
    """
    INPUT
    {'gbm': {'charge': [0.119, 0.174, 0.083, 0.220, ... ],
             'overall': [0.201, 0.256, 0.283, 0.420, ... ],
              ...
            },
     'rf': ...}

    OUTPUT
    fragment_names_list = ['SO2', 'NO2', ...]  - corresponds to number of rows
    {'gbm': {'charge': [(-2, 8), (-12, 3), ... ],   - negative and positive contribution values (sum of all negatives and sum of all positives)
            {'overall': [(0, 10), (-11, 1), ... ],
                         ...
                        ],
     'rf': ...}
    """
    sorted_frag_names, d = prep_data_boxplot(contr_dict, frag_names_list, add_consensus)
    for model_name in d.keys():
        for prop_name in d[model_name].keys():
            tmp = []
            for v in d[model_name][prop_name]:
                n = sum([vv for vv in v if vv < 0])
                p = sum([vv for vv in v if vv > 0])
                tmp.append((n, p))
            d[model_name][prop_name] = tmp
    return sorted_frag_names, d


def main_params(contr_fname, contr_names, fig_fname, model_names, on_screen, min_M, min_N, model_type):

    frag_names, contr_dict = read_contrib_file(contr_fname, model_names, min_M, min_N, contr_names)

    if model_type == 'reg' or model_type == 'class':

        # plot boxplot for single property (contribution)
        if len(contr_names) == 1:

            sorted_frag_names, data = prep_data_boxplot(contr_dict, frag_names)

            if "consensus" in data.keys():
                x_label_position = sort_by_median(data["consensus"][contr_names[0]])
            else:
                x_label_position = sort_by_median(list(data.values())[0][contr_names[0]])

            fig = plt.figure(1)
            fig.suptitle(contr_names[0] + " contributions", fontsize='x-large', color='blue')

            for i, model_name in enumerate(sorted(data.keys(), reverse=True)):

                ax = fig.add_subplot(1, len(data), len(data) - i)
                ax.grid(True, axis='y')
                plt.title(model_name)
                bp = ax.boxplot(x=list(data[model_name][contr_names[0]]),
                                sym="o",
                                positions=x_label_position,
                                vert=False)
                ax.axvline(0, color='grey', linestyle='--')

                # show x labels only on the first plot
                if i == len(data) - 1:
                    ax.set_yticklabels(sorted_frag_names)
                else:
                    ax.set_yticklabels([])
                for box in bp["boxes"]:
                    box.set(color="blue")
                for w in bp["whiskers"]:
                    w.set(color="blue", linestyle="-")
                for flier in bp["fliers"]:
                    flier.set(markeredgecolor="grey",
                              markeredgewidth=0,
                              markerfacecolor="grey",
                              markersize=4)

            # save figure if specified
            # fig.subplots_adjust(bottom=0, left=0, top=1, right=1)
            # fig.tight_layout()
            if fig_fname is not None:
                w = 3.5 * len(data) + 3.5
                h = 0.5 * len(sorted_frag_names) + 1
                fig.set_size_inches(w, h)
                fig.subplots_adjust(wspace=0.05, bottom=0.01, left=0.01, top=(h-1)/h, right=0.99)
                fig.savefig(fig_fname, dpi=300, bbox_inches='tight')
                fig.set_size_inches(8, 6)
            if on_screen:
                plt.show()
            plt.clf()

        # plot bar plot with medians for several properties (contributions)
        else:
            sorted_frag_names, data = prep_data_barplot(contr_dict, frag_names)

            # remove overall to avoid its plotting
            # for v in data.values():
            #     del v["overall"]

            fig = plt.figure(1)
            bar_width = 0.8 / len(list(data.values())[0])
            colors = ['royalblue', 'tomato', 'yellowgreen', 'darkturquoise', 'darkorchid']

            for i, model_name in enumerate(sorted(data.keys(), reverse=True)):
                ax = fig.add_subplot(1, len(data), len(data) - i)
                plt.title(model_name)
                ind = np.arange(len(sorted_frag_names))
                ax.set_yticks(ind + 0.4)
                ax.grid(True, axis='y')

                for j, prop_name in enumerate(data[model_name].keys()):
                    br = ax.barh(bottom=ind + j * bar_width,
                                 width=data[model_name][prop_name],
                                 height=bar_width,
                                 color=colors[j])
                    if prop_name == "CHARGE":
                        br.set_label("electrostatic")
                    if prop_name == "LOGP":
                        br.set_label("hydrophobic")
                    if prop_name == "HB":
                        br.set_label("hydrogen bonding")
                    if prop_name == "REFRACTIVITY":
                        br.set_label("dispersive")
                    if prop_name == "overall":
                        br.set_label("overall")

                if i == len(data) - 1:
                    # reverse the order
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0, 1.12),
                              ncol=1, mode="expand", borderaxespad=0., fontsize='medium')
                    ax.set_yticklabels(sorted_frag_names, rotation=0)
                else:
                    ax.set_yticklabels([])

            fig.tight_layout()
            if fig_fname is not None:
                w = 3.5 * len(data) + 3.5
                h = 0.5 * len(sorted_frag_names) + 1
                fig.set_size_inches(w, h)
                fig.subplots_adjust(wspace=0.05, bottom=0.01, left=0.01, top=(h-1)/h, right=0.99)
                fig.savefig(fig_fname, dpi=300, bbox_inches='tight')
                fig.set_size_inches(8, 6)
            if on_screen:
                plt.show()
            plt.clf()
                
    elif model_type == 'class':  # this is experimental code for integer contributions -1, 0, 1
                                 # it will be never executed

        sorted_frag_names, data = prep_data_barplot_binary_class(contr_dict, frag_names)

        fig = plt.figure(1)
        bar_width = 0.8 / len(list(data.values())[0])
        colors = ['royalblue', 'tomato', 'yellowgreen', 'darkturquoise', 'darkorchid']

        ind = np.arange(len(sorted_frag_names))

        for i, model_name in enumerate(sorted(data.keys(), reverse=True)):
            ax = fig.add_subplot(1, len(data), len(data) - i)
            plt.title(model_name)
            ax.set_yticks(ind + 0.4)
            ax.grid(True, axis='y')

            for j, prop_name in enumerate(data[model_name].keys()):
                br = ax.barh(bottom=ind + j * bar_width,
                             width=[i[1] for i in data[model_name][prop_name]],
                             height=bar_width,
                             color=colors[j])
                if prop_name == "CHARGE":
                    br.set_label("electrostatic")
                if prop_name == "LOGP":
                    br.set_label("hydrophobic")
                if prop_name == "HB":
                    br.set_label("hydrogen bonding")
                if prop_name == "REFRACTIVITY":
                    br.set_label("dispersive")
                if prop_name == "overall":
                    br.set_label("overall")
                br = ax.barh(bottom=ind + j * bar_width,
                             width=[i[0] for i in data[model_name][prop_name]],
                             height=bar_width,
                             color=colors[j])

            if i == len(data) - 1:
                # reverse the order
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0, 1.12),
                          ncol=1, mode="expand", borderaxespad=0., fontsize='medium')
                ax.set_yticklabels(sorted_frag_names, rotation=0)
            else:
                ax.set_yticklabels([])

        fig.tight_layout()
        if fig_fname is not None:
            fig.set_size_inches(3.5 * len(data) + 3.5, 0.5 * len(sorted_frag_names) + 1)
            fig.subplots_adjust(wspace=0.05)
            fig.savefig(fig_fname, dpi=300, bbox_inches='tight')
            fig.set_size_inches(8, 6)
        if on_screen:
            plt.show()


def main():

    parser = argparse.ArgumentParser(description='Plot fragments contributions.')
    parser.add_argument('-i', '--input', metavar='frag_contributions.txt', required=True,
                        help='input text file with calculated fragment contributions.')
    parser.add_argument('-o', '--output', metavar='output.png', default=None,
                        help='output file name for plot.')
    parser.add_argument('-c', '--contributions', metavar='[overall charge logp hb refractivity]', default=['overall'],
                        nargs='*', help='name of contributions to calculate.')
    parser.add_argument('-m', '--models', metavar='[rf gbm svr pls knn]', default=["all"], nargs='*',
                        help='file name of saved models (without extension)')
    parser.add_argument('-s', '--show_in_separate_window', action='store_true', default=False,
                        help='show plot on the screen in a  separate window')
    parser.add_argument('-N', '--min_fragments_count', metavar='', default=1,
                        help='minimal number of fragments in the dataset compounds')
    parser.add_argument('-M', '--min_molecules_count', metavar='', default=1,
                        help='minimal number of compounds containing a fragment')
    parser.add_argument('-t', '--model_type', metavar='reg|class', default='reg',
                        help='models type: reg for regression and class for classification')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": contr_fname = v
        if o == "contributions": contr_names = v
        if o == "output": fig_fname = v
        if o == "models": model_names = v
        if o == "show_in_separate_window": on_screen = v
        if o == "min_fragments_count":
            try:
                min_N = int(v)
            except ValueError:
                print("Incorrect minimal number of fragments. Not a number. Default value (1) will be used.")
                min_N = 1
        if o == "min_molecules_count":
            try:
                min_M = int(v)
            except ValueError:
                print("Incorrect minimal number of molecules containing a fragment. Not a number. "
                      "Default value (1) will be used.")
                min_M = 1
        if o == "model_type": model_type = v

    main_params(contr_fname, contr_names, fig_fname, model_names, on_screen, min_M, min_N, model_type)


if __name__ == '__main__':
    main()
