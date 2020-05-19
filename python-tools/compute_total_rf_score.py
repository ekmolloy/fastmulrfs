import argparse
import dendropy
from compare_two_trees import compare_trees
import os
import sys


def main(args):
    with open(args.stree, 'r') as f:
        temp = f.read()

    total_fp = 0
    total_fn = 0
    total_rf = 0

    with open(args.gtreelist, 'r') as f:
        for l, line in enumerate(f.readlines()):
            taxa = dendropy.TaxonNamespace()

            stree = dendropy.Tree.get(string=temp,
                                      schema='newick',
                                      rooting='force-unrooted',
                                      taxon_namespace=taxa)

            gtree = dendropy.Tree.get(string=line,
                                      schema='newick',
                                      rooting='force-unrooted',
                                      taxon_namespace=taxa)

            [nl, ei1, ei2, fn, fp, rf] = compare_trees(stree, gtree)

            total_fp += fp
            total_fn += fn
            total_rf += rf

    sys.stdout.write('%d,%d,%d\n' % (total_fn, total_fp, total_rf))
    sys.stdout.flush()
    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--stree", type=str,
                        help="Input file containing species tree",
                        required=True)
    parser.add_argument("-g", "--gtreelist", type=str,
                        help="Input file containing gene trees "
                             "(one newick string per line)",
                        required=True)

    main(parser.parse_args())
