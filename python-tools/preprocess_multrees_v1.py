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


def is_binary(tree):
    """
    """
    root = tree.seed_node
    if len(root.child_nodes()) != 3:
        sys.exit("Tree is not binary!")

    for node in tree.preorder_node_iter():
        if (node != root) and (not node.is_leaf()):
            if len(node.child_nodes()) != 2:
                sys.exit("Tree is not binary!")


def build_down_profiles(tree, g2s_map):
    """
    """
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            gene = node.taxon.label
            species = g2s_map[gene]
            node.down_profile = set([species])
            try:
                g2s_map[species] = g2s_map[species] + [gene]
            except KeyError:
                g2s_map[species] = [gene]
        else:
            children = node.child_nodes()
            if len(children) == 2:
                [l, r] = children
                node.down_profile = l.down_profile.union(r.down_profile)
            else:
                [l, m, r] = children
                lr = l.down_profile.union(r.down_profile)
                node.down_profile = lr.union(m.down_profile)


def build_up_profiles(tree):
    """
    """
    root = tree.seed_node
    [rootl, rootm, rootr] = root.child_nodes()
    rootl.up_profile = rootm.down_profile.union(rootr.down_profile)
    rootm.up_profile = rootl.down_profile.union(rootr.down_profile)
    rootr.up_profile = rootl.down_profile.union(rootm.down_profile)

    for node in tree.preorder_node_iter():
        if (node == root) or (node == rootl) or (node == rootm) or (node == rootr):
            pass
        elif node.is_leaf():
            pass
        else:
            parent = node.parent_node
            [pl, pr] = parent.child_nodes()
            if node == pl:
                sibling = pr.down_profile
            else:
                sibling = pl.down_profile
            node.up_profile = parent.up_profile.union(sibling)


def contract_edges_w_invalid_bipartitions(tree):
    """
    """
    for node in tree.postorder_node_iter():
        if node == tree.seed_node:
            pass
        elif node.is_leaf():
            node.edge.length = 1.0
        else:
            test = node.down_profile.intersection(node.up_profile)
            if len(test) != 0:
                node.edge.length = 0.0
            else:
                node.edge.length = 1.0

    tree.collapse_unweighted_edges()

    for edge in tree.edges():
        edge.length = None


def prune_multiple_copies_of_species(tree, g2s_map):
    """
    """
    for l in tree.leaf_nodes():
        x = l.taxon.label
        s = g2s_map[x]
        if x != g2s_map[s][0]:
            l.taxon = None

    tree.prune_leaves_without_taxa()

    # Relabel leaves
    for l in tree.leaf_nodes():
        x = l.taxon.label
        l.taxon.label = g2s_map[x]


def preprocess_multree(tree, g2s_map):
    """
    Preprocesses MUL-tree as described in the FastMulRFS paper
    """
    tree.is_rooted = False
    tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    is_binary(tree)

    build_down_profiles(tree, g2s_map)

    build_up_profiles(tree)

    contract_edges_w_invalid_bipartitions(tree)

    prune_multiple_copies_of_species(tree, g2s_map)


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
    ostr = ""

    g2s_map = read_g2s_map(mfil)

    with open(ifil, 'r') as f:
        for line in f.readlines():
            temp = "".join(line.split())
            tree = dendropy.Tree.get(data=temp,
                                     schema="newick",
                                     rooting='force-unrooted',
                                     preserve_underscores=True)

            preprocess_multree(tree, g2s_map)

            if len([l for l in tree.leaf_nodes()]) > 3:
                ostr += tree.as_string(schema="newick")[5:]

    with open(ofil, 'w') as f:
        f.write(ostr)


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
