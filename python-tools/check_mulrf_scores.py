import argparse
from copy import deepcopy
import dendropy
from preprocess_multrees_v1 import preprocess_multree, read_g2s_map
import os
import sys


def relabel_multree(tree, g2s_map):
    """
    """
    for l in tree.leaf_nodes():
        x = l.taxon.label
        l.taxon.label = g2s_map[x]


def score_with_MulRF(mulrf, stree, gtree, tmpname):
    """
    """
    otree = tmpname + ".tree"
    ofile = tmpname + ".out"
    lfile = tmpname + ".log"

    with open(otree, 'w') as f:
        f.write(stree.as_string(schema="newick").replace("'", ""))
        f.write(gtree.as_string(schema="newick")[5:].replace("'", ""))

    os.system(mulrf + " -i " + otree + " -o " + ofile + " &> " + lfile)

    with open(ofile, 'r') as f:
        line = f.readline()
        words = line.split()
    score = float(words[-1].replace(']', ''))

    os.system("rm " + otree)
    os.system("rm " + ofile)
    os.system("rm " + lfile)

    return score


def check_mulrf_scores(sfile, gfile, mfile, mulrf):
    """
    """
    # Read species tree
    Stree = dendropy.Tree.get(path=sfile,
                              schema="newick",
                              preserve_underscores=True)

    for node in Stree.preorder_node_iter():
        node.label = None
        node.edge.length = None

    # Read gene to species name map
    g2s_map = read_g2s_map(mfile)

    with open(gfile, 'r') as f:
        g = 1
        for line in f.readlines():
            temp = "".join(line.split())

            # Build MUL-tree
            Mtree = dendropy.Tree.get(data=temp,
                                      schema="newick",
                                      rooting="force-unrooted",
                                      preserve_underscores=True)
            Mtree.is_rooted = False
            Mtree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
            relabel_multree(Mtree, g2s_map)

            # Build pre-processed MUL-tree
            MXtree = dendropy.Tree.get(data=temp,
                                       schema="newick",
                                       rooting="force-unrooted",
                                       preserve_underscores=True)
            MXtree.is_rooted = False
            MXtree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
            score_shift = preprocess_multree(MXtree, g2s_map)

            # Compute MulRF scores
            tmpname = gfile.rsplit('.', 1)[0]
            MulRF_M = score_with_MulRF(mulrf, Stree, Mtree, tmpname + "-scored")
            MulRF_MX = score_with_MulRF(mulrf, Stree, MXtree, tmpname + "-preprocessed-and-scored")

            # Check scores match!
            if MulRF_MX + score_shift != MulRF_M:
                sys.stdout.write("Gene tree on line %d failed!" % g)

            g += 1


def main(args):
    check_mulrf_scores(args.stree, args.gtree, args.map, args.mulrf)


if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--stree", type=str,
                        help="Input file containing singly-labeled tree",
                        required=True)
    parser.add_argument("-g", "--gtree", type=str,
                        help="Input file containing gene family trees "
                             "(one newick string per line)",
                        required=True)
    parser.add_argument("-a", "--map", type=str,
                        help="Input file containing label map (each row has form: "
                        "'species_name:gene_name_1,gene_name_2,...')",
                        required=True)
    parser.add_argument("-x", "--mulrf", type=str,
                        help="MulRF binary (e.g., MulRFScorerMac) with full path",
                        required=True)

    main(parser.parse_args())
