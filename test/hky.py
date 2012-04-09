#!/usr/bin/env python
import math, sys
mito_code = {"F": ["UUU", "UUC"],
"T": ["ACU", "ACC", "ACA", "ACG"],
"D": ["GAU", "GAC"],
"L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
"A": ["GCU", "GCC", "GCA", "GCG"],
"E": ["GAA", "GAG",],
"Y": ["UAU", "UAC"],
"C": ["UGU", "UGC"],
"I": ["AUU", "AUC"],
"*": ["UAA", "UAG", "AGA", "AGG"],
"W": ["UGA", "UGG"],
"M": ["AUA", "AUG"],
"H": ["CAU", "CAC"],
"R": ["CGU", "CGC", "CGA", "CGG"],
"V": ["GUU", "GUC", "GUA", "GUG"],
"Q": ["CAA", "CAG"],
"S": ["AGU", "AGC", "UCU", "UCC", "UCA", "UCG"],
"N":["AAU", "AAC"],
"P":[ "CCU", "CCC", "CCA", "CCG"],
"K" :["AAA", "AAG"],
"G":["GGU", "GGC", "GGA", "GGG"],
}
c = [(i,j) for i,j in mito_code.iteritems()]
c.sort()
sorted_code = c[1:] + [c[0]]



from optparse import OptionParser
parser = OptionParser(usage="newick string with numbers as argument", add_help_option=True)
parser.add_option('-k','--kappa',
              dest='kappa',
              default=2.0,
              type="float",
              help="Kappa parameter")
parser.add_option('-a','--pi-a',
              dest='a',
              default=0.25,
              type="float",
              help="Freq of A")
parser.add_option('-c','--pi-c',
              dest='c',
              default=0.25,
              type="float",
              help="Freq of C")
parser.add_option('-g','--pi-g',
              dest='g',
              default=0.25,
              type="float",
              help="Freq of G")
parser.add_option('-b','--branch-length',
              dest='b',
              default=0.25,
              type="float",
              help="Branch length")
parser.add_option('-n','--num-codons',
              dest='n',
              default=100,
              type="int",
              help="Number of codons in codlikestartval command")
parser.add_option('-f','--four-fold-deg',
              dest='f',
              default=False,
              action="store_true",
              help="Output")
parser.add_option('-p','--pmat',
              dest='pmat',
              default=False,
              action="store_true",
              help="Print Prob Mat (if not specified, a Bull CodlikeStartVal command is printed")
parser.add_option('-t','--triplets',
              dest='triplets',
              default=False,
              action="store_true",
              help="Calculate the  Prob Mat for codons")
(options, args) = parser.parse_args()

class Nuc:
    A = 0
    C = 1
    G = 2
    T = 3

pi = tuple([options.a, options.c, options.g, 1.0 - options.a - options.c - options.g])
pi_str = { "A" : pi[Nuc.A],
         "C" : pi[Nuc.C],
         "G" : pi[Nuc.G],
         "U" : pi[Nuc.T],
        }
t = options.b

if pi[Nuc.T] < 0.0:
    sys.exit("The frequency of A, C, and G cannot be greater than 1.0")
if pi[Nuc.A] <= 0.0:
    sys.exit("The frequency of A must be greater than 0.0")
if pi[Nuc.C] <= 0.0:
    sys.exit("The frequency of C must be greater than 0.0")
if pi[Nuc.G] <= 0.0:
    sys.exit("The frequency of G must be greater than 0.0")
kappa = options.kappa
if kappa <= 0.0:
    sys.exit("kappa cannot be negative")


# class freq
classFreq = (pi[Nuc.A] + pi[Nuc.G], 
             pi[Nuc.C] + pi[Nuc.T], 
             pi[Nuc.A] + pi[Nuc.G], 
             pi[Nuc.C] + pi[Nuc.T])


p = [[0.25, 0.25, 0.25, 0.25], 
     [0.25, 0.25, 0.25, 0.25], 
     [0.25, 0.25, 0.25, 0.25], 
     [0.25, 0.25, 0.25, 0.25], ]

bigPiInv = [ 1.0/i for i in classFreq]
beta = 0.5 / ((classFreq[Nuc.G])*(classFreq[Nuc.T]) + kappa*((pi[Nuc.A]*pi[Nuc.G]) + (pi[Nuc.C]*pi[Nuc.T])))
x = math.exp(-beta*t)


td = -beta*(1 + classFreq[Nuc.A]*(kappa - 1.0))
y = math.exp(t*td)
ta = pi[Nuc.A]*(bigPiInv[Nuc.A] - 1.0)
tb = pi[Nuc.G]*bigPiInv[Nuc.A]
tc = pi[Nuc.A]*bigPiInv[Nuc.A]
p[Nuc.A][Nuc.A] = pi[Nuc.A] + (x*ta) + (y*tb)
p[Nuc.C][Nuc.A] = pi[Nuc.A]*(1.0 - x)
p[Nuc.G][Nuc.A] = pi[Nuc.A] + (x*ta) - (y*tc)
p[Nuc.T][Nuc.A] = p[Nuc.C][Nuc.A]


td = -beta*(1 + classFreq[Nuc.C]*(kappa - 1.0))
y = math.exp(t*td)
ta = pi[Nuc.C]*(bigPiInv[Nuc.C] - 1.0)
tb = pi[Nuc.T]*bigPiInv[Nuc.C]
tc = pi[Nuc.C]*bigPiInv[Nuc.C]
p[Nuc.A][Nuc.C] = pi[Nuc.C]*(1.0 - x)
p[Nuc.C][Nuc.C] = pi[Nuc.C] + (x*ta) + (y*tb)
p[Nuc.G][Nuc.C] = p[Nuc.A][Nuc.C]
p[Nuc.T][Nuc.C] = pi[Nuc.C] + (x*ta) - (y*tc)


td = -beta*(1 + classFreq[Nuc.G]*(kappa - 1.0))
y = math.exp(t*td)
ta = pi[Nuc.G]*(bigPiInv[Nuc.G] - 1.0)
tb = pi[Nuc.A]*bigPiInv[Nuc.G]
tc = pi[Nuc.G]*bigPiInv[Nuc.G]
p[Nuc.A][Nuc.G] = pi[Nuc.G] + (x*ta) - (y*tc)
p[Nuc.C][Nuc.G] = pi[Nuc.G]*(1.0 - x)
p[Nuc.G][Nuc.G] = pi[Nuc.G] + (x*ta) + (y*tb)
p[Nuc.T][Nuc.G] = p[Nuc.C][Nuc.G]


td = -beta*(1 + classFreq[Nuc.T]*(kappa - 1.0))
y = math.exp(t*td)
ta = pi[Nuc.T]*(bigPiInv[Nuc.T] - 1.0)
tb = pi[Nuc.C]*bigPiInv[Nuc.T]
tc = pi[Nuc.T]*bigPiInv[Nuc.T]
p[Nuc.A][Nuc.T] = pi[Nuc.T]*(1.0 - x)
p[Nuc.C][Nuc.T] = pi[Nuc.T] + (x*ta) - (y*tc)
p[Nuc.G][Nuc.T] = p[Nuc.A][Nuc.T]
p[Nuc.T][Nuc.T] = pi[Nuc.T] + (x*ta) + (y*tb)



if options.pmat:
    if options.triplets:
        nuc_range = (Nuc.A, Nuc.C, Nuc.G, Nuc.T)
        codon_pmat = []
        for from_f in nuc_range:
            for from_s in nuc_range:
                for from_t in nuc_range:
                    codon_p_row = []
                    for to_f in nuc_range:
                        for to_s in nuc_range:
                            for to_t in nuc_range:
                                x = p[from_f][to_f]*p[from_s][to_s]*p[from_t][to_t]
                                codon_p_row.append(x)
                    codon_pmat.append(codon_p_row)
        fmt = " ".join(["%9.8f"]*64)
        for row in codon_pmat:
            print  fmt % tuple(row)
        
    else:
        for row in p:
            print "%7.6f %7.7f %7.6f %7.6f" % tuple(row)
        
else:
    aa_freq = []
    for i in sorted_code:
        slcode, triplets = i
        f = 0.0
        for triplet in triplets:
            first_base = triplet[0]
            second_base = triplet[1]
            third_base = triplet[2]
            f += pi_str[first_base] * pi_str[second_base] * pi_str[third_base]
        if slcode != "*":
            aa_freq.append(f)
        else:
            stop_freq = f
    print "[ non-stop = %f ] [ all = %f ]" % ( sum(aa_freq), sum(aa_freq) + stop_freq)
    aa_freq_str = "(%s)" %" ".join(["%f" % i for i in aa_freq])
    ncods = options.n
    aa_freq_str_list = [aa_freq_str] * ncods
    aa_freq_str = "%d %s" % (ncods, "\n                ".join(aa_freq_str_list))
    print """codlikestartval currbranchlengths groupskip=0
            lastParamImprovement = 0.0120391 lastBranchImprovement = 0.00940639
            PARAMIMPROVEThisround = 0.00378464
            nst=6 basefreq = ( %f %f %f )
            rMatrix = ( 1.0 %f 1.0 1.0 %f 1.0)
            treescale = 6.40052
            aafreq = (
                %s
                );""" % (pi[Nuc.A], pi[Nuc.C], pi[Nuc.G], kappa, kappa, aa_freq_str)

