#!/usr/bin/env python
# 2019-7-17
# Find A/B compartment according RNA-seq data!
import sys
#from os import system
from optparse import OptionParser
import csv

parser = OptionParser()
parser.add_option("-g", "--gtf", action="store", type="string", dest="gtf_file", help="gtf file")
parser.add_option("-e", "--eigen", action="store", type="string", dest="eigen_file",help="eigen file")
parser.add_option("-o","--output", action="store", type="string", dest="output_file", help="output file")
(options, args) = parser.parse_args()

if not (options.gtf_file and options.output_file and options.eigen_file):
  parser.print_help()
  sys.exit(0)

#------------------------------------------------------------------------------------------------
# gene number in each chromosome
def align(Eigen, Gene):
    number = 0
    for start, end, value in Eigen:
        for promoter in Gene:
            if start <= promoter <= end:
                number += 1
    return number

#------------------------------------------------------------------------------------------------
# correct
def corr(Bool, Bed):
    Cor = []
    if Bool == 0:
        for start, end, value in Bed:
            Cor.append((start, end, value))
    if Bool == 1:
        for start, end, value in Bed:
            Cor.append((start, end, value*-1))
    return Cor

#------------------------------------------------------------------------------------------------
# correct eigen in each chromosome
def ab(EIGEN, GENE, CHR):
    Eigen_NaN = []
    Eigen_Neg = []
    Eigen_Pos = []
    Gene = []
    for line in EIGEN:
        if line[0] == CHR and line[3] == 0:
            Eigen_NaN.append((line[1], line[2], "NaN"))
        if line[0] == CHR and line[3] > 0:
            Eigen_Pos.append((line[1], line[2], line[3]))
        if line[0] == CHR and line[3] < 0:
            Eigen_Neg.append((line[1], line[2], line[3]))
    for gline in  GENE:
        if gline[0] == CHR:
            Gene.append(gline[1])
        
    pos_num = align(Eigen_Pos, Gene)
    neg_num = align(Eigen_Neg, Gene)
    if pos_num >= neg_num:
        Bool = 0
    if pos_num < neg_num:
        Bool = 1
        
    correct_eigen_pos = corr(Bool, Eigen_Pos)
    correct_eigen_neg = corr(Bool, Eigen_Neg)
        
    Eigenvector = correct_eigen_pos + correct_eigen_neg + Eigen_NaN
    return Eigenvector
        
    
#------------------------------------------------------------------------------------------------
# gtf file to get chromosome, start, and end
with open(options.gtf_file, 'r') as gtf:
    gene = []
    chromosome = []
    lines = csv.reader(gtf, delimiter = '\t')
    #next(lines)
    for line in lines:
        if line[0][0] != "#" and line[0] not in chromosome:
            chromosome.append(line[0])
        if line[0][0] != "#" and  line[2] == "gene" and line[6] == "+":
            gene.append((line[0], int(line[3])))
        if line[0][0] != "#" and  line[2] == "gene" and line[6] == "-":
            gene.append((line[0], int(line[4])))

#------------------------------------------------------------------------------------------------
# engenvector file to get chromosome, start, end, and value
with open(options.eigen_file, 'r') as fe:
    eigen = []
    lines = csv.reader(fe, delimiter = '\t')
    for line in lines:
        if line[3] == "NaN":
            line[3] = 0
        eigen.append((line[0], int(line[1]), int(line[2]), float(line[3])))    

#------------------------------------------------------------------------------------------------
# write result
Result = open(options.output_file, 'w')

for chrom in chromosome:
    EIGEN = ab(eigen, gene, chrom)
    EIGEN_sort =  sorted(EIGEN, key=lambda x: x[1])
    for line in EIGEN_sort:
        Result.write("%s\t%s\t%s\t%s\n" %(chrom, line[0], line[1], line[2]))

Result.close()

