import argparse
from preprocess_multrees_v2 import compute_score_shift
from preprocess_multrees_v2 import preprocess_multree
from preprocess_multrees_v2 import read_label_map
import os
import sys
import treeswift


def relabel_tree_by_species(tree, g2s_map):
    """
    Relabels leaves in tree by species (previously labeled by gene copies)

    Parameters
    ----------
    tree : treeswift tree object
    g2s_map : dictionary
              maps gene copy labels to species labels
    """
    for leaf in tree.traverse_leaves():
        temp = leaf.get_label()
        leaf.set_label(g2s_map[temp])


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


def strip_extra(tree):
    """
    Remove edge lengths and internal node label from tree

    Parameters
    ----------
    tree : treeswift tree object
    """
    for node in tree.traverse_preorder():
        node.edge_length = None
        if not node.is_leaf():
            node.label = None


def check_mulrf_scores(sfile, gfile, mfile, mulrf):
    """
    Checks RF scores are the same regardless of preprocessing gene family trees

    Parameters
    ----------
    sfile : string
            name of file containing species tree
    gfile : string
            name of file containing gene family trees
    mfile : string
            name of file containing map between gene copy and species labels
    mulrf: string
           name including full path of MulRFScorer binary
    """
    # Read species tree
    stree = treeswift.read_tree(sfile, "newick")
    strip_extra(stree)

    # Read gene to species name map
    [g2s_map, s2g_map] = read_label_map(mfile)

    total_rf = 0

    with open(gfile, 'r') as f:
        g = 1
        for line in f.readlines():
            temp = "".join(line.split())

            # Build MUL-tree
            mtree = treeswift.read_tree(temp, "newick")
            mtree.deroot()
            strip_extra(mtree)

            relabel_tree_by_species(mtree, g2s_map)

            # Build pre-processed MUL-tree
            mxtree =  treeswift.read_tree(temp, "newick")
            mxtree.deroot()
            strip_extra(mxtree)

            [nEM, nLM, nR, c, nEMX, nLMX] = preprocess_multree(mxtree,
                                                               g2s_map,
                                                               s2g_map)

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
    check_mulrf_scores(args.stree, args.gtree, args.map, args.mulrf)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--stree", type=str,
                        help="Input file containing singly-labeled tree",
                        required=True)
    parser.add_argument("-g", "--gtree", type=str,
                        help="Input file containing gene family trees "
                             "(one newick string per line)",
                        required=True)
    parser.add_argument("-a", "--map", type=str,
                        help="Input file containing label map; "
                             "each row has form: "
                             "'species_name:gene_name_1,gene_name_2,...'",
                        required=True)
    parser.add_argument("-x", "--mulrf", type=str,
                        help="MulRFScorer binary including full path",
                        required=True)

    main(parser.parse_args())
