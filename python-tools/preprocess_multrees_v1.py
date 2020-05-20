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


def prune_multiple_copies_of_species(tree, g2s_map, s2g_map):
    """
    Removes all but one leaf with the same species label

    Parameters
    ----------
    tree : dendropy tree object
    g2s_map : dictionary
              maps gene copy labels to species labels
    s2g_map : dictionary
              maps species label to gene copy labels
    """
    found_duplicate = set([])
    nLMX = 0
    c = 0

    for leaf in tree.leaf_nodes():
        gene = leaf.taxon.label
        species = g2s_map[gene]
        all_genes = s2g_map[species]

        if gene == all_genes[0]:
            leaf.taxon.label = species
            nLMX += 1
        else:
            leaf.taxon = None
            if not (species in found_duplicate):
                found_duplicate.add(species)
                c += 1

    tree.prune_leaves_without_taxa()

    return [nLMX, c]


def preprocess_multree(tree, g2s_map, s2g_map):
    """
    Preprocesses MUL-tree as described in the FastMulRFS paper

    Parameters
    ----------
    tree : dendropy tree object
    g2s_map : dictionary
              maps gene copy labels to species labels
    s2g_map : dictionary
              maps species label to gene copy labels
    """
    tree.is_rooted = False
    tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    build_down_profiles(tree, g2s_map)
    build_up_profiles(tree)

    [nLM, nEM, nR, nO] = contract_edges_w_invalid_bipartitions(tree)
    [nLMX, c] = prune_multiple_copies_of_species(tree, g2s_map, s2g_map)

    nEMX = nO + nLMX

    return [nEM, nLM, nR, c, nEMX, nLMX]


def compute_score_shift(nEM, nLM, nR, c, nEMX, nLMX):
    """
    Compute constant shift for RF score as described in FastMulRFS paper

    Parameters
    ----------
    nEM : int
          Number of edges in MUL-tree
    nLM : int
          Number of leaves in MUL-tree
    nR : int
         Number of edges in MUL-tree that induce invalid bipartitions
         (i.e., edges split the label set into two non-disjoint sets)
    c : int
        Number of species with multiple copies in the MUL-tree
    nEMX : int
           Number of edges in preprocessed MUL-tree
    nLMX : int
           Number of leaves in preprocessed MUL-tree,
           which is the same as the number of species

    Returns
    -------
    Constant shift for RF score as described in FastMulRFS paper
    """
    return nLMX + c + nEM - nEMX - (2 * nR) - nLM


def read_label_map(ifile):
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
    g2s_map : dictionary
              maps gene copy labels to species labels
    s2g_map : dictionary
              maps species label to gene copy labels
    """
    g2s_map = {}
    s2g_map = {}

    with open(ifile, 'r') as f:
        for line in f.readlines():
            [species, genes] = line.split(':')
            genes = genes.split(',')
            genes[-1] = genes[-1].replace('\n', '')

            s2g_map[species] = genes

            for gene in genes:
                if gene == species:
                    sys.exit("Error: Gene copy label cannot be the species "
                             "label!")
                g2s_map[gene] = species

    return [g2s_map, s2g_map]


def read_preprocess_and_write_multrees(ifile, mfile, ofile, verbose):
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
    [g2s_map, s2g_map] = read_label_map(mfile)

    with open(ifile, 'r') as fi, open(ofile, 'w') as fo:
        g = 1
        for line in fi.readlines():
            if verbose:
                sys.stdout.write("Analyzing gene tree on line %d...\n" % g)
                sys.stdout.flush()

            temp = "".join(line.split())
            tree = dendropy.Tree.get(data=temp,
                                     schema="newick",
                                     rooting="force-unrooted",
                                     preserve_underscores=True)

            [nEM, nLM, nR, c, nEMX, nLMX] = preprocess_multree(tree,
                                                               g2s_map,
                                                               s2g_map)

            score_shift = compute_score_shift(nEM, nLM, nR, c, nEMX, nLMX)

            if nLMX > 3:
                fo.write(tree.as_string(schema="newick")[5:])
            else:
                if verbose:
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

    read_preprocess_and_write_multrees(args.input,
                                       args.map,
                                       output,
                                       args.verbose)


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
    parser.add_argument("--verbose", action="store_true")

    main(parser.parse_args())
