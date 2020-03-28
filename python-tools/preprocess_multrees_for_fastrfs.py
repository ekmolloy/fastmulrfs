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


def read_g2s_map(ifil):
    """
    Parameters
    ----------
    ifil : string
           name of input gene to species label map file (ASTRAL-multi)

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


def is_binary(tree):
    nodes = [n for n in tree.preorder_node_iter()]
    for node in nodes[1:]:
        if not node.is_leaf():
            children = node.child_nodes()
            if len(children) != 2:
                sys.exit("Tree is not binary!")


def transform_multrees(ifil, mfil, ofil):
    """
    Creates file for FastRFS as described in FastMulRFS paper

    Parameters
    ----------
    ifil : string
           name of input gene tree file (one newick string per line)
    mfil : string
           name of input gene to species label map file (ASTRAL-multi)
    ofil : string
           name of output file (one newick string per line)
    """
    ostr = ""

    g2sm = read_g2s_map(mfil)
    s2gm = {}

    with open(ifil, 'r') as f:
        for line in f.readlines():
            temp = "".join(line.split())
            tree = dendropy.Tree.get(data=temp,
                                     schema="newick",
                                     rooting='force-unrooted',
                                     preserve_underscores=True)

            # Randomly root tree!
            is_binary(tree)
            tree.resolve_polytomies(limit=2, update_bipartitions=False)

            # Create down profiles
            for node in tree.postorder_node_iter():
                if node.is_leaf():
                    gene = node.taxon.label
                    species = g2sm[gene]
                    node.down_profile = set([species])
                    try:
                        s2gm[species] = s2gm[species] + [gene]
                    except KeyError:
                        s2gm[species] = [gene]
                else:
                    [l, r] = node.child_nodes()
                    node.down_profile = l.down_profile.union(r.down_profile)

            # Create up profiles
            nodes = [n for n in tree.preorder_node_iter()]
            root = tree.seed_node
            [rootl, rootr] = root.child_nodes()
            rootl.up_profile = rootr.down_profile
            rootr.up_profile = rootl.down_profile

            #test = rootl.down_profile.intersection(rootr.down_profile)
            #if len(test) != 0:
            #    sys.exit("Duplication at root!")

            for node in nodes:
                if (node == root) or (node == rootl) or (node == rootr):
                    pass
                elif node.is_leaf():
                    pass
                else:
                    p = node.parent_node
                    [pl, pr] = p.child_nodes()
                    if node == pl:
                        x = pr.down_profile
                    else:
                        x = pl.down_profile
                    node.up_profile = p.up_profile.union(x)

            # Contract edges
            for node in tree.postorder_node_iter():
                if node == root:
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

            # Prune leaves
            for l in tree.leaf_nodes():
                x = l.taxon.label
                s = g2sm[x]
                if x != s2gm[s][0]:
                    l.taxon = None

            tree.prune_leaves_without_taxa()

            # Relabel leaves
            for l in tree.leaf_nodes():
                x = l.taxon.label
                l.taxon.label = g2sm[x]

            # Unroot tree
            tree.is_rooted = False
            tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)

            if len([l for l in tree.leaf_nodes()]) > 3:
                ostr += tree.as_string(schema="newick")[5:]

    with open(ofil, 'w') as f:
        f.write(ostr)


def main(args):
    base = args.input.rsplit('.', 1)
    prefix = base[0]
    suffix = base[1]
    output = base[0] + "-for-fastrfs." + base[1]
    transform_multrees(args.input, args.map, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file", required=True)
    parser.add_argument("-a", "--map", type=str,
                        help="Input file", required=True)

    main(parser.parse_args())
