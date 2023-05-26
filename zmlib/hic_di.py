#!/usr/bin/env python

# transfer Hi-C matrix to the form that can be used by Ren Bin pip!
import os, sys
from os import system
from optparse import OptionParser
import csv
parser = OptionParser()
parser.add_option("-a", "--abs", action="store", type="string", dest="abs_file", help="chromosome bin file that is generated in HiC-Pro")
parser.add_option("-o", "--output", action="store", type="string", dest="output_file", help="File for store result")
parser.add_option("-m", "--matrix", action="store", type="string", dest="matrix_file", help="hic density matrix file")
parser.add_option("-c", "--chrom", action="store", type="string", dest="chrom_position", help="the chromsome need for processing")
(options, args) = parser.parse_args()
usage = """
  -h, --help            show this help message and exit
  -a, --abs_file        chromosome bin file that is generated in HiC-Pro
  -m, --matrix          hic density matrix file
  -c, --chromsome       the chromsome need for processing
  -o, --output          the final output file
"""
if not(options.abs_file and options.output_file and options.matrix_file and options.chrom_position):
  print(usage)
  print('Example:','\033[1;31mpython hic_di.py -a WT_abs.bed -m chr5.matrix -c chr 5 -o di_density.matrix\033[0m')
  sys.exit(0)

chromosome = []
with open(options.abs_file, 'r') as fbed:
    bedlines = csv.reader(fbed, delimiter= '\t')
    for line in bedlines:
        if line[0] == options.chrom_position:
            chromosome.append((line[0], line[1], line[2]))   

di_density = []    
with open(options.matrix_file, 'r') as fd:
    density_lines = fd.readlines()
    if len(density_lines) != len(chromosome):
        print("some thing wrang!")
        sys.exit(0)
    for i in range(len(chromosome)):
        chrom = chromosome[i][0] 
        start = chromosome[i][1]
        end = chromosome[i][2]
        mat = density_lines[i].strip()
        di_density.append((chrom, start, end, mat))
        
with open(options.output_file,'w') as fo:
    for line in di_density:
        fo.write("%s\t%s\t%s\t%s\n" %(line[0], line[1], line[2], line[3]))
