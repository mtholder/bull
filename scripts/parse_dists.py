#!/usr/bin/env python
from optparse import OptionParser
import sys
parser = OptionParser(usage="", add_help_option=True)
parser.add_option('-c','--columns',  
              dest='column',
              default=False,
              action="store_true",
              help="column output"
              )
parser.add_option('-n','--num',  
              dest='num',
              default=0,
              type="int",
              help="""The number of branches in the model tree.  
                      Assumes the tree dist layout was is:
                          model tree
                          model tree
                          star tree
                          "num" x lines with just 1 branch
                          model tree
                          result tree
                          star tree
                          a tree  for each branch in the result tree
                      followed by the trees being tested.""")   
(options, args) = parser.parse_args()
dfn=args[0]
df = open(dfn, "rU")
lines = df.readlines()
source = iter(lines)
spl = source.next().strip().split()
if spl[0] != "tree":
    sys.exit('Expecting first line to start with "tree", found "%s"' % spl[0])
spl = [int(i) for i in source.next().strip().split()]
if spl[1] != 0:
    sys.exit('Second line to contain just a treedist of 0')
    
spl = [int(i) for i in source.next().strip().split()]
if spl[1] != options.num:
    sys.exit('Expecting the tree 2 to tree 1 distance to be %d, but found %d' % (options.num, spl[1]))

first_branch_offset = 3
star_column_index = first_branch_offset - 1

# lines corresponding to the trees that are the true branches
for j in range(options.num):
    spl = [int(i) for i in source.next().strip().split()]
    if spl[1] != options.num - 1:
        sys.exit('Expecting the tree %d to tree 1 distance to be %d, but found %d' % (first_branch_offset + j, options.num, spl[1]))

spl = [int(i) for i in source.next().strip().split()]
if spl[1] != 0:
    sys.exit('Expecting the tree %d to tree 1 distance to be 0, but found %d' % (options.num+ first_branch_offset, spl[1]))
try:
    spl = [int(i) for i in source.next().strip().split()]
except StopIteration:
    sys.exit("No trees returned.")

rf = spl[1]
nedges_returned = spl[star_column_index]
#print "nedges_returned = ", nedges_returned
missed_edge_rf = nedges_returned + 1
found_edge_rf = nedges_returned - 1
espl = spl[star_column_index+1:-2]
fn = n_missed_edge_rf = espl.count(missed_edge_rf)
fp = rf - fn
n_found_edge_rf = espl.count(found_edge_rf)
if n_found_edge_rf + n_missed_edge_rf != options.num:
    sys.exit("%d + %d != %d "%(n_found_edge_rf, n_missed_edge_rf, options.num))
if options.column:
    sys.stdout.write("%d %d %d" % (rf, fn, fp))
else:
    print "RF", rf
    print "FN", fn
    print "FP", fp
nfec = 0
for n, rfe in enumerate(espl):
    if rfe == found_edge_rf:
        if options.column:
            sys.stdout.write(" 1")
        else:
            print "T%d 1" % (n+1)
        nfec += 1
    else:
        if options.column:
            sys.stdout.write(" 0")
        else:
            print "T%d 0" % (n+1)
assert n_found_edge_rf == nfec
spl = [int(i) for i in source.next().strip().split()]
if spl[1] != options.num:
    sys.exit('Expecting the tree %d to tree 1 distance to be %d, but found %d' % (options.num + 6, options.num, spl[1]))
true_edge_rf = options.num - 1

nfec = 0
for n in range(nedges_returned):
    rfe = int(source.next().strip().split()[1])
    if rfe == true_edge_rf:
        if options.column:
            sys.stdout.write(" 1")
        else:
            print "R%d 1" % (n+1)
        nfec += 1
    else:
        if options.column:
            sys.stdout.write(" 0")
        else:
            print "R%d 0" % (n+1)
assert n_found_edge_rf == nfec
