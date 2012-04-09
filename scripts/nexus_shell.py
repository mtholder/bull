#!/usr/bin/env python

import sys
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-n", "--ntax", dest="ntax", default=4, 
        type="int",
        help="The number of taxa to generate")
    parser.add_option("-t", "--trees", dest="treesb", default=False, 
        action="store_true",
        help="If true a trees block will be opened (with a translate block).  The block will have no trees commands and will need to be terminated.")
        
    (options, args) = parser.parse_args()
    if args:
        sys.exit("No arguments accepted.  Use the -h flag to see the list of command line options (flags) that are accepted")
    ntax = options.ntax
    if ntax < 0:
        sys.exit("ntax must be greater than 0")
    trees = options.treesb
    taxa_names = ["t%d" % i for i in range(1, ntax + 1)]
    print("""#NEXUS
BEGIN TAXA;
  Dimensions NTax = %d ;
  TaxLabels %s ;
END;
""" % (ntax, " ".join(taxa_names)))
    if trees:
        sys.stdout.write("""BEGIN TREES;
  Translate""")
        for i, tn in enumerate(taxa_names):
            if i > 0:
                sys.stdout.write(" ,")
            sys.stdout.write("\n    %d %s" % (i+1, tn))
        print(" ;")
