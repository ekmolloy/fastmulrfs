"""
This file is a python prototype of Algorithm 1 from the paper:
Molloy, E.K., Warnow, T. (2020). FastMulRFS: Fast and accurate species tree
    estimation under generic gene duplication and loss models.
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
    Annotates edge above each node with an 'down profile', i.e., the set of
    species below the edge

    Parameters
    ----------
    tree : dendropy tree object
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
    Annotates edge above each node with an 'up profile', i.e., the set of
    species above the edge

    NOTE: Must be called after build_down_profiles()

    Parameters
    ----------
    tree : dendropy tree object
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
    Contracts edges that do not induce valid bipartitions

    NOTE: Must be called after build_down_profiles() and build_up_profiles()

    Parameters
    ----------
    tree : dendropy tree object
    """
    nLM = 0
    nX = 0
    nR = 0
    nO = 0

    for node in tree.postorder_node_iter():
        if node == tree.seed_node:
            pass
        elif node.is_leaf():
            node.edge.length = 1.0
            nLM += 1
        else:
            test = node.down.intersection(node.up)
            if len(test) != 0:
                nX += 1
                node.edge.length = 0.0
            else:
                if (len(node.down) == 1) or (len(node.up) == 1):
                    nR += 1
                else:
                    nO += 1
                node.edge.length = 1.0

    tree.collapse_unweighted_edges()

    for edge in tree.edges():
        edge.length = None

    nEM = nLM + nX + nR + nO

    return [nLM, nEM, nR, nO]


def prune_multiple_copies_of_species(tree, g2s_map):
    """
    Removes all but one leaf with the same species label

    Parameters
    ----------
    tree : dendropy tree object
    g2s_map : dictionary
              maps gene copy labels to species labels
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

    nLMX = 0
    for l in tree.leaf_nodes():
        temp = l.taxon.label
        l.taxon.label = g2s_map[temp]
        nLMX += 1

    return [nLMX, c]


def preprocess_multree(tree, g2s_map):
    """
    Preprocesses MUL-tree as described in the FastMulRFS paper

    Parameters
    ----------
    tree : dendropy tree object
    g2s_map : dictionary
              maps gene copy labels to species labels
    """
    tree.is_rooted = False
    tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    build_down_profiles(tree, g2s_map)
    build_up_profiles(tree)

    [nLM, nEM, nR, nO] = contract_edges_w_invalid_bipartitions(tree)
    [nLMX, c] = prune_multiple_copies_of_species(tree, g2s_map)
    nEMX = nO + nLMX

    score_shift = nLMX + c + nEM - nEMX - (2 * nR) - nLM

    return score_shift


def read_g2s_map(ifile):
    """
    Reads file containing map from gene copy to species labels into dictionary

    Parameters
    ----------
    ifile : string
            name of file containing map between gene copy and species labels
            each row has form:
            species_name:gene_name_1,gene_name_2,...

    Returns
    -------
    g2s_map : python dictionary
              maps gene copy labels to species labels
    """
    g2s_map = {}
    with open(ifile, 'r') as f:
        for line in f.readlines():
            [species, genes] = line.split(':')
            genes = genes.split(',')
            genes[-1] = genes[-1].replace('\n', '')
            for gene in genes:
                g2s_map[gene] = species
    return g2s_map


def read_preprocess_and_write_multrees(ifile, mfile, ofile):
    """
    Creates file with preprocessed MUL-trees for FastRFS

    Parameters
    ----------
    ifile : string
            name of file containing gene family trees
            (one newick string per line)
    mfile : string
            name of file containing label map file; each row has form:
            species_name:gene_name_1,gene_name_2,...
    ofile : string
            name of output file (one newick string per line)
    """
    g2s_map = read_g2s_map(mfile)

    with open(ifile, 'r') as fi, open(ofile, 'w') as fo:
        g = 1
        for line in fi.readlines():
            temp = "".join(line.split())
            tree = dendropy.Tree.get(data=temp,
                                     schema="newick",
                                     rooting="force-unrooted",
                                     preserve_underscores=True)

            score_shift = preprocess_multree(tree, g2s_map)

            if len(tree.leaf_nodes()) > 3:
                fo.write(tree.as_string(schema="newick")[5:])
            else:
                sys.stdout.write("Removing gene tree on line %d "
                                 "as it has three or fewer taxa!\n" % g)

            g += 1


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
                        help="Input file containing label map; "
                             "each row has form: "
                             "'species_name:gene_name_1,gene_name_2,...')",
                        required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file name",
                        required=False)

    main(parser.parse_args())
