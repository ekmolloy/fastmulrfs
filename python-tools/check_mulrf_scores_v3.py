import argparse
from preprocess_multrees_v3 import compute_score_shift
from preprocess_multrees_v3 import preprocess_multree
from preprocess_multrees_v3 import unroot
import os
import sys
import treeswift


def score_with_MulRF(mulrf, stree, gtree, temp):
    """
    Uses MulRF to compute the RF distance between a species tree and a gene
    family tree

    Parameters
    ----------
    mulrf : string
            name including full path of MulRFScorer binary
    stree : treeswift tree object
            species tree
    gtree : treeswift tree object
            gene family tree
    temp : string
           name for creating temporary files

    Returns
    -------
    score : integer
            RF distance between species tree and gene tree
    """
    ifile = temp + ".tree"
    ofile = temp + ".out"
    lfile = temp + ".log"

    with open(ifile, 'w') as f:
        f.write(stree.newick())
        f.write('\n')
        f.write(gtree.newick())
        f.write('\n')

    os.system(mulrf + " -i " + ifile + " -o " + ofile + " &> " + lfile)

    with open(ofile, 'r') as f:
        line = f.readline()
    words = line.split()
    score = float(words[-1].replace(']', ''))

    os.system("rm " + ifile)
    os.system("rm " + ofile)
    os.system("rm " + lfile)

    return score


def remove_internal_node_labels(tree):
    """
    Remove internal node label from tree before running MulRF

    Parameters
    ----------
    tree : treeswift tree object
    """
    for node in tree.traverse_preorder():
        node.edge_length = None
        if not node.is_leaf():
            node.label = None


def check_mulrf_scores(sfile, gfile, mulrf):
    """
    Checks RF scores are the same regardless of preprocessing gene family trees

    Parameters
    ----------
    sfile : string
            name of file containing species tree
    gfile : string
            name of file containing gene family trees
    mulrf: string
           name including full path of MulRFScorer binary
    """
    # Read species tree
    stree = treeswift.read_tree(sfile, "newick")
    remove_internal_node_labels(stree)
    stree.suppress_unifurcations()

    total_rf = 0

    with open(gfile, 'r') as f:
        g = 1
        for line in f.readlines():
            temp = "".join(line.split())

            # Build MUL-tree
            mtree = treeswift.read_tree(temp, "newick")
            remove_internal_node_labels(mtree)
            unroot(mtree)

            # Build pre-processed MUL-tree
            mxtree = treeswift.read_tree(temp, "newick")
            remove_internal_node_labels(mxtree)

            [nEM, nLM, nR, c, nEMX, nLMX] = preprocess_multree(mxtree)

            score_shift = compute_score_shift(nEM, nLM, nR, c, nEMX, nLMX)

            # Compute MulRF scores
            temp = gfile.rsplit('.', 1)[0]
            mscore = score_with_MulRF(mulrf, stree, mtree,
                                      temp + "-scored")
            mxscore = score_with_MulRF(mulrf, stree, mxtree,
                                       temp + "-preprocessed-and-scored")

            # Check scores match!
            if mxscore + score_shift != mscore:
                sys.exit("Gene tree on line %d failed!\n" % g)

            total_rf += mscore

            g += 1

    sys.stdout.write('%d\n' % total_rf)
    sys.stdout.flush()
    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


def main(args):
    check_mulrf_scores(args.stree, args.gtree, args.mulrf)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--stree", type=str,
                        help="Input file containing singly-labeled tree",
                        required=True)
    parser.add_argument("-g", "--gtree", type=str,
                        help="Input file containing gene family trees "
                             "(one newick string per line)",
                        required=True)
    parser.add_argument("-x", "--mulrf", type=str,
                        help="MulRFScorer binary including full path",
                        required=True)

    main(parser.parse_args())
