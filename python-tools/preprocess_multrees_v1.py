"""
This file is a python prototype of Algorithm 1 from the paper:
Molloy, E.K., Warnow, T. (2020). FastMulRFS: Fast and accurate species tree estimation under generic gene duplication and loss models.
https://doi.org/10.1101/835553

Copyright (c) 2020 Erin K. Molloy
All rights reserved.

License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import argparse
import dendropy
import sys


def build_down_profiles(tree, g2s_map):
    """
    """
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            gene = node.taxon.label
            species = g2s_map[gene]
            node.down = set([species])
            try:
                g2s_map[species] = g2s_map[species] + [gene]
            except KeyError:
                g2s_map[species] = [gene]
        else:
            node.down = set([])
            for child in node.child_nodes():
                node.down = node.down.union(child.down)


def build_up_profiles(tree):
    """
    """
    root = tree.seed_node
    children_of_root = root.child_nodes()
    for node in children_of_root:
        node.up = set([])
        for sibl in children_of_root:
            if node != sibl:
                node.up = node.up.union(sibl.down)
    children_of_root = set(children_of_root)

    for node in tree.preorder_node_iter():
        if (node == root) or (node in children_of_root):
            pass
        elif node.is_leaf():
            pass
        else:
            parent = node.parent_node
            node.up = parent.up
            for sibl in parent.child_nodes():
                if node != sibl:
                    node.up = node.up.union(sibl.down)


def contract_edges_w_invalid_bipartitions(tree):
    """
    """
    L_M = 0
    X = 0
    R = 0
    O = 0

    for node in tree.postorder_node_iter():
        if node == tree.seed_node:
            pass
        elif node.is_leaf():
            node.edge.length = 1.0
            L_M += 1
        else:
            test = node.down.intersection(node.up)
            if len(test) != 0:
                X += 1
                node.edge.length = 0.0
            else:
                if (len(node.down) == 1) or (len(node.up) == 1):
                    R += 1
                else:
                    O += 1
                node.edge.length = 1.0

    tree.collapse_unweighted_edges()

    for edge in tree.edges():
        edge.length = None

    E_M = L_M + X + R + O

    return [L_M, E_M, R, O]


def prune_multiple_copies_of_species(tree, g2s_map):
    """
    """
    found = set([])
    c = 0
    for l in tree.leaf_nodes():
        gene = l.taxon.label
        species = g2s_map[gene]
        if gene != g2s_map[species][0]:
            l.taxon = None
            if not (species in found):
                found.add(species)
                c += 1

    tree.prune_leaves_without_taxa()

    for l in tree.leaf_nodes():
        x = l.taxon.label
        l.taxon.label = g2s_map[x]

    L_MX = len(tree.leaf_nodes())

    return [L_MX, c]


def preprocess_multree(tree, g2s_map):
    """
    Preprocesses MUL-tree as described in the FastMulRFS paper
    """
    tree.is_rooted = False
    tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    build_down_profiles(tree, g2s_map)
    build_up_profiles(tree)

    [L_M, E_M, R, O] = contract_edges_w_invalid_bipartitions(tree)
    [L_MX, c] = prune_multiple_copies_of_species(tree, g2s_map)
    E_MX = O + L_MX

    print("S = %d, c = %d, E_M = %d, E_MX = %d, R = %d, L_M = %d" % (L_MX, c, E_M, E_MX, R, L_M))
    score_shift = L_MX + c + E_M - E_MX - (2 * R) - L_M

    return score_shift


def read_g2s_map(ifil):
    """
    Parameters
    ----------
    ifil : string
           name of input label map file; each row has form:
           species_name:gene_name_1,gene_name_2,...

    Returns
    -------
    g2sm : python dictionary
           maps gene labels to species labels
    """
    g2sm = {}
    with open(ifil, 'r') as f:
        for line in f.readlines():
            [species, genes] = line.split(':')
            genes = genes.split(',')
            genes[-1] = genes[-1].replace('\n', '')
            for gene in genes:
                g2sm[gene] = species
    return g2sm


def read_preprocess_and_write_multrees(ifil, mfil, ofil):
    """
    Creates file with preprocessed MUL-trees for FastRFS

    Parameters
    ----------
    ifil : string
           name of input gene family tree file (one newick string per line)
    mfil : string
           name of input label map file; each row has form:
           species_name:gene_name_1,gene_name_2,...
    ofil : string
           name of output file (one newick string per line)
    """
    g2s_map = read_g2s_map(mfil)

    with open(ifil, 'r') as fi, open(ofil, 'w') as fo:
        for line in fi.readlines():
            temp = "".join(line.split())
            tree = dendropy.Tree.get(data=temp,
                                     schema="newick",
                                     rooting='force-unrooted',
                                     preserve_underscores=True)

            score_shift = preprocess_multree(tree, g2s_map)

            if len([l for l in tree.leaf_nodes()]) > 3:
                fo.write(tree.as_string(schema="newick")[5:])
            else:
                f.write("\n")


def main(args):
    if args.output is None:
        base = args.input.rsplit('.', 1)
        prefix = base[0]
        suffix = base[1]
        output = base[0] + "-for-fastrfs." + base[1]
    else:
        output = args.output

    read_preprocess_and_write_multrees(args.input, args.map, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file containing gene family trees "
                             "(one newick string per line)",
                        required=True)
    parser.add_argument("-a", "--map", type=str,
                        help="Input file containing label map (each row has form: "
                             "'species_name:gene_name_1,gene_name_2,...')",
                        required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file name",
                        required=False)

    main(parser.parse_args())
