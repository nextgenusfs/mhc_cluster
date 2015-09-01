#!/usr/bin/env python

#This script runs trims FASTQ to set length, Quality Filters, filters for DRB locus, and clusters using USEARCH8
#written by Jon Palmer palmer.jona at gmail dot com

import os
import argparse
import subprocess
import inspect
import re
import multiprocessing
from Bio import SeqIO
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='mhc-OTU_cluster.py', usage="%(prog)s [options] -f file.demux.fq\n%(prog)s -h for help menu",
    description='''Script processes data from mhc-process_reads.  This script runs quality filtering on the reads, filters the reads for contamination by using an MHC DRB HMM model, and then clusters into OTUs using USEARCH.
    Requires USEARCH: http://drive5.com/usearch and HMMER3: http://hmmer.janelia.org''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-f','--fastq', dest="FASTQ", required=True, help='FASTQ file (Required)')
parser.add_argument('-o','--out', default='out', help='Base output name')
parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
parser.add_argument('-l','--length', default='auto', help='Trim Length')
parser.add_argument('-p','--pct_otu', default=99, help="OTU Clustering Percent")
parser.add_argument('-m','--minsize', default='2', help='Min size to keep for clustering')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
parser.add_argument('--translate', action='store_true', help='Translate OTUs to protein space')
args=parser.parse_args()

#find cpus, use 1 less than total
cpus = multiprocessing.cpu_count() - 1
cpus = str(cpus)

#open log file for usearch8 stderr redirect
log_name = args.out + '.EE' + args.maxee + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)
log_file = open(log_name, 'ab')

usearch = args.usearch
try:
    subprocess.call([usearch, '--version'], stdout = log_file, stderr = log_file)
except OSError:
    print "%s not found in your PATH, exiting." % usearch 
    os._exit(1)
    
#now run usearch8 fastq filtering step, output to fasta
filter_out = args.out + '.EE' + args.maxee + '.filter.fa'
print "\nCMD: Quality Filtering\n%s -fastq_filter %s -fastq_maxee %s -fastq_qmax 45 -fastaout %s\n" % (usearch, args.FASTQ, args.maxee, filter_out)
subprocess.call([usearch, '-fastq_filter', args.FASTQ, '-fastq_maxee', args.maxee, '-fastq_qmax', '45', '-fastaout', filter_out], stdout = log_file, stderr = log_file)

#now run HMMer3 to filter contaminant DRB sequences out.
hmm_out = args.out + '.EE' + args.maxee + '.DRB.hmm.txt'
hmm = script_path + '/lib/MHC_DNA.hmm'
FNULL = open(os.devnull, 'w')
print "CMD: Running HMMER3 using MHC DRB HMM model (using %s cpus)\nhmmscan --cpu %s --domtblout %s %s %s\n" % (cpus, cpus, hmm_out, hmm, filter_out)
subprocess.call(['hmmscan', '--cpu', cpus, '--domtblout', hmm_out, hmm, filter_out], stdout = FNULL, stderr = FNULL)

#now parse HMMer results
print "CMD: Filtering HMM results"
#now parse the hmmer results
q_list = []
l_list = []
counter = 0
hmmer_results = open(hmm_out, 'r')
for qresult in SearchIO.parse(hmmer_results, "hmmscan3-domtab"):
    q_list.append(qresult.id)
    l_list.append(qresult.seq_len)
    counter += 1
hmmer_results.close()

sum = sum(l_list)
avg_len = sum / counter

print "%10u sequences passed filter" % counter
print "%10u bp is average length" % avg_len

if args.length == 'auto':
    trim_len = avg_len
else:
    trim_len = int(args.length)
    
print "\nCMD: Retrieving results and trimming/padding to %i bp for clustering\n" % trim_len
#now retrieve seqs in the "pass" list by looping through the query list, make same length
pass_out = args.out + '.EE' + args.maxee + '.hmm.pass.fa'
pass_handle = open(pass_out, 'wb')
filtered = SeqIO.parse(filter_out, "fasta")
for rec in filtered:
    if rec.id in q_list:
        L = len(rec.seq)
        if L < trim_len:
            Seq = rec.seq + (trim_len - L)*'N'
        else:
            T = trim_len - 1
            Seq = rec.seq[:T]
        pass_handle.write(">%s\n%s\n" % (rec.id, Seq))
pass_handle.close()

#now run usearch8 full length dereplication
derep_out = args.out + '.EE' + args.maxee + '.derep.fa'
print "CMD: De-replication\n%s -derep_fulllength %s -sizeout -fastaout %s\n" % (usearch, pass_out, derep_out)
subprocess.call([usearch, '-derep_fulllength', pass_out, '-sizeout', '-fastaout', derep_out], stdout = log_file, stderr = log_file)

#now run usearch 8 sort by size
sort_out = args.out + '.EE' + args.maxee + '.sort.fa'
print "CMD: Sorting by Size\n%s -sortbysize %s -minsize %s -fastaout %s\n" % (usearch, derep_out, args.minsize, sort_out)
subprocess.call([usearch, '-sortbysize', derep_out, '-minsize', args.minsize, '-fastaout', sort_out], stdout = log_file, stderr = log_file)

#now run clustering algorithm
radius = str(100 - int(args.pct_otu))
otu_out = args.out + '.EE' + args.maxee + '.otus.fa'
print "CMD: Clustering OTUs\n%s -cluster_otus %s -sizein -sizeout -relabel MHC_ -otu_radius_pct %s -otus %s" % (usearch, sort_out, radius, otu_out)
subprocess.call([usearch, '-cluster_otus', sort_out, '-sizein', '-sizeout', '-relabel', 'MHC_', '-otu_radius_pct', radius, '-otus', otu_out], stdout = log_file, stderr = log_file)

#Fix OTUs, remove trailing N's
fix_otus = args.out + '.EE' + args.maxee + '.fixed.otus.fa'
fix_handle = open(fix_otus, 'wb')
fix = open(otu_out, 'rb')
otu_count = 0
for rec in SeqIO.parse(fix, "fasta"):
    otu_count += 1
    Seq = re.sub('[^GATC]', "", str(rec.seq).upper())
    fix_handle.write(">%s\n%s\n" % (rec.id, Seq))
fix_handle.close()
fix.close()
print "%10u total OTUs\n" % otu_count
    
#now map reads back to OTUs
uc_out = args.out + '.EE' + args.maxee + '.mapping.uc'
mapping_pct = str(float(args.pct_otu) / 100)
#mapping_pct = '0.97'
print "CMD: Mapping Reads to OTUs\n%s -usearch_global %s -strand plus -id %s -db %s -uc %s\n" % (usearch, pass_out, mapping_pct, fix_otus, uc_out)
subprocess.call([usearch, '-usearch_global', pass_out, '-strand', 'plus', '-id', mapping_pct, '-db', fix_otus, '-uc', uc_out], stdout = log_file, stderr = log_file)

#Build OTU table
otu_table = args.out + '.EE' + args.maxee + '.otu_table.txt'
uc2tab = script_path + "/lib/uc2otutab.py"
print "CMD: Creating OTU Table\npython %s %s > %s" % (uc2tab, uc_out, otu_table)
os.system('%s %s %s %s %s' % ('python', uc2tab, uc_out, '>', otu_table))

#translate to protein space (optional)
if args.translate:
    trans_out = args.out + '.EE' + args.maxee + '.proteins.fa'
    print "\nCMD: Translating to Protein Space\n%s -fastx_findorfs %s -aaout %s -orfstyle 7 -mincodons 40\n" % (usearch, fix_otus, trans_out)
    subprocess.call([usearch, '-fastx_findorfs', fix_otus, '-aaout', trans_out, '-orfstyle', '7', '-mincodons', '40'], stdout = log_file, stderr = log_file)
    
    #HMM against translated amino acids
    trans_hmm = args.out + '.EE' + args.maxee + '.proteins.hmm.txt'
    hmm_prot = script_path + '/lib/MHC.hmm'
    print "CMD: Running HMMER3 MHC HMM model (using %s cpus)\nhmmscan --cpu %s --tblout %s %s %s" % (cpus, cpus, trans_hmm, hmm_prot, trans_out)
    subprocess.call(['hmmscan', '--cpu', cpus, '--tblout', trans_hmm, hmm_prot, trans_out], stdout = FNULL, stderr = FNULL)
    
    #now filter results for best hit
    hmmer_prots = open(trans_hmm, 'r')
    hit_list = []
    for qresult in SearchIO.parse(hmmer_prots, "hmmer3-tab"):
        hits = qresult.hits
        if len(hits) > 0:
            hit_list.append(hits[0].query_id)
    hmmer_prots.close()
    
    #now filter the fasta file to get only hits that have MHC domain
    pass_prot = args.out + '.EE' + args.maxee + '.proteins.pass.fa'
    pass_file = open(pass_prot, 'wb')
    prot_filtered = SeqIO.parse(trans_out, "fasta")
    for rec in prot_filtered:
        if rec.id in hit_list:
            pass_file.write(">%s\n%s\n" % (rec.id, rec.seq))
    pass_file.close()

#Print location of files to STDOUT
print "\n------------------------------------------------"
print "OTU Clustering Script has Finished Successfully"
print "------------------------------------------------"
print ("Input FASTQ:           %s" % (args.FASTQ))
print ("Filtered FASTQ:        %s" % (filter_out))
print ("HMM Pass FASTA:        %s" % (pass_out))
print ("Dereplicated FASTA:    %s" % (derep_out))
print ("Sorted FASTA:          %s" % (sort_out))
print ("Clustered OTUs:        %s" % (fix_otus))
if args.translate:
    print ("Translated OTUs:       %s" % (pass_prot))
print ("UCLUST Mapping file:   %s" % (uc_out))
print ("OTU Table:             %s" % (otu_table))
print ("USEARCH LogFile:       %s" % (log_name))
print "------------------------------------------------"

 