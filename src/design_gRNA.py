#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:53:04 2019

@author: pjvonk
"""
#this program requires the following things:
#connectFGDB in PYTHONPATH
#biopython installed
#bowtie and bowtie_build in PATH (NOTE: this is bowtie, not bowtie2)
#dump_data.py in PATH

import connectFGDB
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
import glob
import shutil

#Getting all the information
parser = argparse.ArgumentParser()
parser.add_argument('-m', help='Mysql config file', required = True)
parser.add_argument('-d', help='Database name of species', required = True)
parser.add_argument('-p', help='target_proteinId', required = True)
parser.add_argument('-o', help='Output directory (cannot exist yet)', required = True)
parser.add_argument('-r', help='Comma separated coordinates of target for gRNA relative to start and end of the gene (default -10,100). If using a negative number, place between '' and add a space before the - sign', required = False, default = '-10,100')
parser.add_argument('--upoffset', help='Offset of the upflank. Positive values shift search range to the right, negative values to the left. Recommended to set the same as in design_KO', required = False, default = 0, type = int)
parser.add_argument('--downoffset', help='Offset of the downflank. Positive values shift search range to the right, negative values to the left. Recommended to set the same as in design_KO', required = False, default = 0, type = int)
#parser.add_argument('-b', help='Location of bowtie (default /usr/bin)', required = False, default = '/usr/bin')

#Define variables
args = parser.parse_args()
mysql_config = args.m
database = args.d
gene = args.p
location = args.r.split(',')
location = [int(x) for x in location]
output_dir = args.o
up_offset = args.upoffset
down_offset = args.downoffset
#bowtie_location = args.b

# Confirm that bowtie is in the PATH. If it is, then extract the directory of the executable. This is needed by CCTOP. NOTE: this is bowtie, not bowtie2!
bowtie_path = shutil.which("bowtie")
if not bowtie_path:
    raise Exception("the executable 'bowtie' is not in your PATH variable")
bowtie_location = os.path.dirname(bowtie_path)

script_dir = os.path.dirname(os.path.realpath(__file__))

#create output folder and temportary folder
if os.path.exists(output_dir):
    raise Exception('Output directory %s already exists' % output_dir)
os.mkdir(output_dir)
os.chdir(output_dir)
os.mkdir('run_files')

#make build if build is not yet available
if not glob.glob('%s/index/%s*' % (script_dir, database)):
    #os.system('python %s/build_genome_index.py -m %s -d %s' % (script_dir, mysql_config, database))
    os.system('dump_data.py -d %s -m %s -t assembly -o %s/index/%s.fasta' % (database, mysql_config, script_dir, database))
    os.system('bowtie-build %s/index/%s.fasta %s/index/%s' % (script_dir, database, script_dir, database))

#create the required fasta file
fasta_file = open('run_files/target_regions.fasta', 'w')
db = connectFGDB.FGDB(database, mysql_config, preload_sequences=False, preload_assembly=True)

#Create fasta file of the target regions
gene_start = db.genes[gene]['cdsCoordinates'][0][0]
gene_end = db.genes[gene]['cdsCoordinates'][-1][1]
scaffold = db.genes[gene]['scaffold']
strand = db.genes[gene]['strand']
gene_sequence = db.assembly[scaffold][gene_start:gene_end]

#get sequence_region for prediction
left_candidate_region = db.assembly[scaffold][gene_start+location[0]+up_offset:gene_start+location[1]+up_offset]
right_candidate_region = db.assembly[scaffold][gene_end-location[1]+down_offset:gene_end-location[0]+down_offset]
gene_underscore = gene.split('|')[1]
fasta_file.write('>%s\n%s\n\n>%s\n%s\n\n' % (gene_underscore + '_left', left_candidate_region, gene_underscore + '_right', right_candidate_region))

# write the genomic locus genbank file
adjecent_genes = db.find_adjecent_genes(gene, 2, 5000)
left_border = db.genes[adjecent_genes[0]]['cdsCoordinates'][0][0] - 500
right_border = db.genes[adjecent_genes[-1]]['cdsCoordinates'][-1][1] + 500

if strand == '-':
    gene_strand = -1
else:
    gene_strand = +1
adjecent_genes.remove(gene)

#make the genbank file and annotate it
sequence_genomic_locus = db.assembly[scaffold][left_border:right_border]
sequence_object = Seq(sequence_genomic_locus, IUPAC.ambiguous_dna)
sequence_record = SeqRecord(sequence_object, id = gene_underscore, name = gene_underscore, description = 'Scaffold region of %s' % gene_underscore)
feature = SeqFeature(FeatureLocation(start = gene_start-left_border, end = gene_end - left_border),
                                         type = 'promotor', strand = gene_strand, qualifiers = {'label': gene_underscore})
sequence_record.features.append(feature)
for adjecent_gene in adjecent_genes:
    adjecent_gene_start = db.genes[adjecent_gene]['cdsCoordinates'][0][0]
    adjecent_gene_end = db.genes[adjecent_gene]['cdsCoordinates'][-1][1]
    adjecent_strand = db.genes[adjecent_gene]['strand']
    if adjecent_strand == '-':
        adjecent_gene_strand = -1
    else:
        adjecent_gene_strand = +1    
    adjecent_gene_underscore = adjecent_gene.split('|')[0] + '_' + adjecent_gene.split('|')[1]
    feature = SeqFeature(FeatureLocation(start = adjecent_gene_start-left_border, end = adjecent_gene_end - left_border),
                                         type = 'gene', strand = adjecent_gene_strand, qualifiers = {'label': adjecent_gene_underscore})
    sequence_record.features.append(feature)
genbank_output = open('%s.gb' % gene_underscore, 'w')
SeqIO.write(sequence_record, genbank_output, 'genbank')
genbank_output.close()  
fasta_file.close()

#run cctop
cctop_command = 'python %s/cctop_standalone/CCTop.py --input run_files/target_regions.fasta --index %s/index/%s --bowtie %s --pam NGG --sgRNA5 NN --totalMM 4 --output run_files' % (script_dir, script_dir, database, bowtie_location)
os.system(cctop_command)

#parse results and output the required files
result_file = open('result.txt', 'w')

result_file.write('ID\tproteinId\tSide\tScore\toff_targets\tgRNA\tfw_primer\trv_primer\n')
i = 1
gene_underscore = gene.split('|')[1]
left_result = open('run_files/%s_left.xls' % gene_underscore, 'r')
right_result = open('run_files/%s_right.xls' % gene_underscore, 'r')
for line in left_result.readlines()[9:]:
    if line.split('\t')[0] == 'T' + str(i):
        off_target_count = -1
        gRNA = line.split('\t')[1][0:20]
        seq_gRNA = Seq(gRNA)
        rc_gRNA = str(seq_gRNA.reverse_complement())
        score = line.split('\t')[4]
        fw_primer = 'TAATACGACTCACTATAG' + gRNA
        rv_primer = 'TTCTAGCTCTAAAAC' + rc_gRNA
        i+=1
    if line.split('\t')[0].split('_')[0] == 'scaffold':
        off_target_count += 1
    if line == '\n':
        result_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('T'+str(i-1), gene, 'left', score.split(' ')[0], str(off_target_count), gRNA, fw_primer, rv_primer))
result_file.write('\n')
result_file.write('ID\tproteinId\tSide\tScore\toff_targets\tgRNA\tfw_primer\trv_primer\n')
i = 1
for line in right_result.readlines()[9:]:
    if line.split('\t')[0] == 'T' + str(i):
        off_target_count = -1
        gRNA = line.split('\t')[1][0:20]
        seq_gRNA = Seq(gRNA)
        rc_gRNA = str(seq_gRNA.reverse_complement())
        score = line.split('\t')[4]
        fw_primer = 'TAATACGACTCACTATAG' + gRNA
        rv_primer = 'TTCTAGCTCTAAAAC' + rc_gRNA
        i+=1
    if line.split('\t')[0].split('_')[0] == 'scaffold':
        off_target_count += 1
    if line == '\n':
        result_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('T'+str(i-1), gene, 'right', score.split(' ')[0], str(off_target_count), gRNA, fw_primer, rv_primer))
result_file.write('\n')
left_result.close()
right_result.close()
result_file.close()
