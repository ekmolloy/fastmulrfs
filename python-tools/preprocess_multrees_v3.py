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
import sys
import treeswift


def count_leaves(tree):
    """
    Count number of leaves in tree

    Parameters
    ----------
    tree : treeswift tree object

    Returns number of leaves in tree
    """
    return len([l for l in tree.traverse_leaves()])


def unroot(tree):
    """
    Unroot tree

    Parameters
    ----------
    tree : treeswift tree object
    """
    tree.suppress_unifurcations()

    for node in tree.traverse_preorder():
        if node.is_root():
            root = node
            break

    if root.num_children() == 2:
        [left, right] = root.child_nodes()

        if left.is_leaf():
            right.contract()
        else:
            left.contract()


def build_down_profiles(tree):
    """
    Annotates edge above each node with an 'down profile', i.e., the set of
    species below the edge

    Parameters
    ----------
    tree : treeswift tree object
    """
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.down = set([node.get_label()])
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
    tree : treeswift tree object
    """
    # Find root node
    for node in tree.traverse_preorder():
        if node.is_root():
            root = node
            root.skip = True
        elif node.is_leaf():
            node.skip = True
        else:
            node.skip = False

    # Compute up profiles for children of root
    children_of_root = root.child_nodes()
    for node in children_of_root:
        node.up = set([])
        node.skip = True
        for sibl in children_of_root:
            if node != sibl:
                node.up = node.up.union(sibl.down)
    children_of_root = set(children_of_root)

    # Compute up profiles for remaining nodes
    for node in tree.traverse_preorder():
        label = node.get_label()
        if not node.skip:
            parent = node.get_parent()
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
    tree : treeswift tree object
    """
    nLM = 0
    nX = 0
    nR = 0
    nO = 0

    for node in tree.traverse_postorder():
        if node.is_root():
            node.edge_length = None
        elif node.is_leaf():
            node.edge_length = None
            nLM += 1
        else:
            test = node.down.intersection(node.up)
            if len(test) != 0:
                nX += 1
                node.edge_length = 0.0
            else:
                if (len(node.down) == 1) or (len(node.up) == 1):
                    nR += 1
                else:
                    nO += 1
                node.edge_length = None

    for node in tree.traverse_postorder():
        if node.edge_length is not None:
            node.contract()

    nEM = nLM + nX + nR + nO

    return [nLM, nEM, nR, nO]


def prune_multiple_copies_of_species(tree):
    """
    Removes all but one leaf with the same species label

    Parameters
    ----------
    tree : treeswift tree object
    """
    found = set([])
    found_duplicate = set([])
    nLMX = 0
    c = 0

    for leaf in tree.traverse_leaves():
        species = leaf.get_label()

        if species not in found:
            found.add(species)
            nLMX += 1
        else:
            leaf.contract()
            if species not in found_duplicate:
                found_duplicate.add(species)
                c += 1

    tree.suppress_unifurcations()

    return [nLMX, c]


def preprocess_multree(tree):
    """
    Preprocesses MUL-tree as described in the FastMulRFS paper

    Parameters
    ----------
    tree : treeswift tree object
    """
    unroot(tree)

    build_down_profiles(tree)

    build_up_profiles(tree)

    [nLM, nEM, nR, nO] = contract_edges_w_invalid_bipartitions(tree)

    [nLMX, c] = prune_multiple_copies_of_species(tree)

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


def read_preprocess_and_write_multrees(ifile, ofile, verbose):
    """
    Creates file with preprocessed MUL-trees for FastRFS

    Parameters
    ----------
    ifile : string
            name of file containing gene family trees
            (one newick string per line)
    ofile : string
            name of output file (one newick string per line)
    """
    with open(ifile, 'r') as fi, open(ofile, 'w') as fo:
        g = 1
        for line in fi.readlines():
            if verbose:
                sys.stdout.write("Preprocessing gene tree on line %d...\n" % g)
                sys.flush()

            temp = "".join(line.split())

            if not temp:
                donot = 1
            else:
                tree = treeswift.read_tree_newick(temp)

                toofew = False
                if count_leaves(tree) < 4:
                    dotnot = 2
                else:
                    [nEM, nLM, nR, c, nEMX, nLMX] = preprocess_multree(tree)
                    score_shift = compute_score_shift(nEM, nLM, nR, c, nEMX,
                                                      nLMX)

                    if nLMX < 4:
                        donot = 3
                    else:
                        fo.write(tree.newick() + '\n')

                if donot and verbose:
                    sys.stdout.write("...did not write tree as ")
                    if donot == 1:
                        sys.stdout.write("as line is empty!")
                    elif donot == 2:
                        sys.stdout.write("as tree has <4 leaves before "
                                         "preprocessing!")
                    elif donot == 3:
                        sys.stdout.write("as tree has <4 leaves after "
                                         "preprocessing!")
                    sys.stdout.write('\n')
                    sys.stdout.flush()

            g += 1


def main(args):
    read_preprocess_and_write_multrees(args.input, args.output, args.verbose)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file containing gene family trees "
                             "(one newick string per line)",
                        required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file name",
                        required=True)
    parser.add_argument("--verbose", action="store_true")

    main(parser.parse_args())
