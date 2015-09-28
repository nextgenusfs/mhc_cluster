#!/usr/bin/env python

#This script runs trims FASTQ to set length, Quality Filters, filters for DRB locus, and clusters using USEARCH8
#written by Jon Palmer palmer.jona at gmail dot com

import os
import argparse
import subprocess
import inspect
import re
import multiprocessing
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO
    from Bio import SeqIO
    
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
parser.add_argument('-n','--num_diff', default='False', help="OTU Clustering Number of differences")
parser.add_argument('-m','--minsize', default='2', help='Min size to keep for clustering')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
parser.add_argument('--translate', action='store_true', help='Translate OTUs')
args=parser.parse_args()

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]


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

if args.num_diff != 'False':
    print "\nWarning:  --num_diff %s was specified, this will override the --pct_otu option" % args.num_diff
            
    
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

#run UNOISE on dereplicated data
#unoise_out = args.out + '.EE' + args.maxee + '.denoised.fa'
#print "CMD: Denoising Data with UNOISE\n%s -cluster_fast %s -centroids %s -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size\n" % (usearch, derep_out, unoise_out)
#subprocess.call([usearch, '-cluster_fast', derep_out, '-centroids', unoise_out, '-id', '0.9', '-maxdiffs', '5', '-abskew', '10', '-sizein', '-sizeout', '-sort', 'size'], stdout = log_file, stderr = log_file)

#now run usearch 8 sort by size
sort_out = args.out + '.EE' + args.maxee + '.sort.fa'
print "CMD: Sorting by Size\n%s -sortbysize %s -minsize %s -fastaout %s\n" % (usearch, derep_out, args.minsize, sort_out)
subprocess.call([usearch, '-sortbysize', derep_out, '-minsize', args.minsize, '-fastaout', sort_out], stdout = log_file, stderr = log_file)

#now run clustering algorithm
if args.num_diff == 'False':
    radius = str(100 - float(args.pct_otu))
    mapping_pct = str(float(args.pct_otu) / 100)
else:
    diff = int(args.length) - int(args.num_diff)
    percent = float(diff) / float(args.length)
    radius = str(100 - (100 * percent))
    mapping_pct = "{:0.3f}".format(percent)
otu_out = args.out + '.EE' + args.maxee + '.otus.fa'
print "CMD: Clustering OTUs\n%s -cluster_otus %s -sizein -sizeout -relabel MHC_ -otu_radius_pct %s -otus %s" % (usearch, sort_out, radius, otu_out)
subprocess.call([usearch, '-cluster_otus', sort_out, '-sizein', '-relabel', 'MHC_', '-otu_radius_pct', radius, '-otus', otu_out], stdout = log_file, stderr = log_file)

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
print "CMD: Mapping Reads to OTUs\n%s -usearch_global %s -strand plus -id %s -db %s -uc %s\n" % (usearch, pass_out, mapping_pct, fix_otus, uc_out)
subprocess.call([usearch, '-usearch_global', pass_out, '-strand', 'plus', '-id', mapping_pct, '-db', fix_otus, '-uc', uc_out], stdout = log_file, stderr = log_file)

#Build OTU table
otu_table = args.out + '.EE' + args.maxee + '.otu_table.txt'
uc2tab = script_path + "/lib/uc2otutab.py"
print "CMD: Creating OTU Table\npython %s %s > %s" % (uc2tab, uc_out, otu_table)
os.system('%s %s %s %s %s' % ('python', uc2tab, uc_out, '>', otu_table))

#translate to protein space (optional)
if args.translate:
    print "\nCMD: Translating to Amino Acid Sequence\n"
    trans_out = args.out + '.EE' + args.maxee + '.proteins.fa'
    trans_file = open(trans_out, "w")
    trans_in = open(fix_otus, 'r')
    for rec in SeqIO.parse(trans_in, "fasta"):
        trans_file.write(">%s;frame_1\n" % (rec.id))
        trans_file.write("%s\n" % (rec.seq.translate(to_stop=True)))
        trans_file.write(">%s;frame_2\n" % (rec.id))
        trans_file.write("%s\n" % (rec.seq[1:].translate(to_stop=True)))
        trans_file.write(">%s;frame_3\n" % (rec.id))
        trans_file.write("%s\n" % (rec.seq[2:].translate(to_stop=True)))
    trans_in.close()
    trans_file.close()
 
    #HMM against translated amino acids
    trans_hmm = args.out + '.EE' + args.maxee + '.proteins.hmm.txt'
    hmm_prot = script_path + '/lib/MHC.hmm'
    print "\nCMD: Running HMMER3 MHC HMM model (using %s cpus)\nhmmscan --cpu %s --tblout %s %s %s" % (cpus, cpus, trans_hmm, hmm_prot, trans_out)
    subprocess.call(['hmmscan', '--cpu', cpus, '--domtblout', trans_hmm, hmm_prot, trans_out], stdout = FNULL, stderr = FNULL)
    
    #now filter results for best hit, and get alignment coordinates
    hmmer_prots = open(trans_hmm, 'r')
    hit_list = {}
    for qresult in SearchIO.parse(hmmer_prots, "hmmscan3-domtab"):
        hits = qresult.hits
        if len(hits) > 0:
            hit_list['%s' % hits[0].query_id] = []
            hit_list['%s' % hits[0].query_id].append(hits[0].hsps[0].query_start)
            hit_list['%s' % hits[0].query_id].append(hits[0].hsps[0].query_end)
    hmmer_prots.close()

    #now filter the fasta file to get only hits that have MHC domain
    pass_prot = args.out + '.EE' + args.maxee + '.proteins.pass.fa'
    pass_file = open(pass_prot, 'wb')
    SeqRecords = SeqIO.to_dict(SeqIO.parse(trans_out, "fasta"))
    for key in sorted(hit_list.keys(), key=natural_sort_key):
        start = hit_list[key][0]
        end = hit_list[key][1]
        subseq = SeqRecords[key][start:end].seq
        pass_file.write(">%s;%s-%s\n%s\n" % (key, start, end, subseq))
    pass_file.close()
    
    #finally concatenate duplicated sequences
    concat_out = args.out + '.EE' + args.maxee + '.proteins.unique.fa'
    keep = {}
    total_count = 0
    concat_file = open(concat_out, 'wb')
    Seq = SeqIO.parse(pass_prot, "fasta")
    for rec in Seq:
        total_count += 1
        sequence=str(rec.seq).upper()
        if sequence not in keep:
            keep[sequence]=rec.id
        else:
            keep[sequence]+="|"+rec.id
    flipKeep = {y:x for x,y in keep.iteritems()}
    unique_count = 0
    for key in sorted(flipKeep, key=natural_sort_key):
        unique_count += 1
        concat_file.write(">"+key+"\n"+flipKeep[key]+"\n")
    concat_file.close()
    print "%10u total proteins" % total_count
    print "%10u unique proteins" % unique_count
    
#Print location of files to STDOUT
print "\n------------------------------------------------"
print "OTU Clustering Script has Finished Successfully"
print "------------------------------------------------"
print ("Input FASTQ:           %s" % (args.FASTQ))
print ("Filtered FASTQ:        %s" % (filter_out))
print ("HMM Pass FASTA:        %s" % (pass_out))
print ("Dereplicated FASTA:    %s" % (derep_out))
#print ("De-noised FASTA:       %s" % (unoise_out))
print ("Sorted FASTA:          %s" % (sort_out))
print ("Clustered OTUs:        %s" % (fix_otus))
if args.translate:
    print ("Translated OTUs:       %s" % (pass_prot))
    print ("Translated Unique:     %s" % (concat_out))
print ("USEARCH Mapping file:  %s" % (uc_out))
print ("OTU Table:             %s" % (otu_table))
print ("USEARCH LogFile:       %s" % (log_name))
print "------------------------------------------------"

 