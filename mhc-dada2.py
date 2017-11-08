#!/usr/bin/env python

import sys, os, argparse, logging, shutil, subprocess, inspect, warnings, itertools, math
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,currentdir)
import lib.lib as lib
import lib.revcomp_lib as revcomp_lib
import numpy as np
import pandas as pd
from natsort import natsorted
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO, SearchIO

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 0.5]"%(x,))
    return x

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-dada2.py',
    description='''Script takes output from AMPtk pre-processing and runs DADA2 MHC clustering''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', required=True, help='Input Demuxed containing FASTQ')
parser.add_argument('-o','--out', default='dada2', help='Output Basename')
parser.add_argument('-l','--length', type=int, help='Length to truncate reads')
parser.add_argument('-e','--maxee', default='1.0', help='MaxEE quality filtering')
parser.add_argument('-m','--mapping_file', help='Metadata file: QIIME mapping format can have extra meta data columns')
parser.add_argument('--pool', default='on', choices = ['on', 'off'], help='Pool all sequences together for DADA2')
parser.add_argument('--platform', default='ion', choices=['ion', 'illumina', '454'], help='Sequencing platform')
parser.add_argument('--uchime_ref', help='Run UCHIME REF, fasta db file')
parser.add_argument('--min_protlen', default=75, type=int, help='Minimum length of amino acid sequence length')
parser.add_argument('--min_reads_sample', default=500, type=int, help='Minimum number of reads per sample')
parser.add_argument('--filter_allele', default='2', help='Remove alelles if not found in more than X samples')
parser.add_argument('--index_bleed', default=0.001, type=restricted_float, help='Remove counts below threshold. 0.005 = 0.5%%')
parser.add_argument('--rev_comp', action='store_true', help='Reverse complement reads')
parser.add_argument('--debug', action='store_true', help='Keep all intermediate files')
args=parser.parse_args()

dada2script = os.path.join(currentdir, 'lib', 'dada2_pipeline_nofilt.R')

def folder2list(input, ending):
    names = []
    if not os.path.isdir(input):
        return False
    else:
        for x in os.listdir(input):
            if x.endswith(ending):
                x = os.path.join(input, x)
                names.append(x)
    return names

def splitDemux2(input, outputdir):
    for title, seq, qual in FastqGeneralIterator(open(input)):
        sample = title.split('barcodelabel=')[1]
        sample = sample.replace(';', '')
        if not args.length:
            with open(os.path.join(outputdir, sample+'.fastq'), 'ab') as output:
                output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        else:
            if len(seq) >= int(args.length):
                with open(os.path.join(outputdir, sample+'.fastq'), 'ab') as output:
                    output.write("@%s\n%s\n+\n%s\n" % (title, seq[:int(args.length):], qual[:int(args.length)]))
def getAvgLength(input):
    AvgLength = []
    for title, seq, qual in FastqGeneralIterator(open(input)):
        AvgLength.append(len(seq))
    Average = sum(AvgLength) / float(len(AvgLength))
    Min = min(AvgLength)
    Max = max(AvgLength)
    a = np.array(AvgLength)
    nintyfive = np.percentile(a, 5)
    return (Average, Min, Max, int(nintyfive))

#remove logfile if exists
log_name = args.out + '.mhc-dada2.log'
if os.path.isfile(log_name):
    lib.removefile(log_name)

lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
#initialize script, log system info and usearch version
lib.SystemInfo()

#check dependencies
programs = ['Rscript', 'vsearch']
lib.CheckDependencies(programs)

#create tmpdir to store intermediate files
tmpdir = args.out+'_'+str(os.getpid())
os.makedirs(tmpdir)

#reverse complement data if argument passed
SeqIn = os.path.join(tmpdir, args.out+'.input.fq')
if args.rev_comp:
    lib.log.info("Reverse complementing FASTQ records")
    with open(SeqIn, 'w') as revcompout:
        for title, seq, qual in FastqGeneralIterator(open(args.fastq)):
            SeqRC = revcomp_lib.RevComp(seq)
            QualRC = qual[::-1]
            revcompout.write("@%s\n%s\n+\n%s\n" % (title, SeqRC, QualRC))
else:
    shutil.copyfile(args.fastq, SeqIn)

#Count FASTQ records and remove 3' N's as dada2 can't handle them
lib.log.info("Loading FASTQ Records")
orig_total = lib.countfastq(SeqIn)
size = lib.checkfastqsize(SeqIn)
readablesize = lib.convertSize(size)
lib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

#quality filter
lib.log.info("Quality Filtering, expected errors < %s" % args.maxee)
derep = os.path.join(tmpdir, args.out+'.qual-filtered.fq')
filtercmd = ['vsearch', '--fastq_filter', SeqIn, '--fastq_maxee', str(args.maxee), '--fastqout', derep, '--fastq_qmax', '55', '--fastq_maxns', '0']
lib.runSubprocess(filtercmd, lib.log)
total = lib.countfastq(derep)
lib.log.info('{0:,}'.format(total) + ' reads passed')

'''
#Get Average length without any N's
averageLen = getAvgLength(derep)
lib.log.info("DADA2 compatible read lengths, avg: %i bp, min: %i bp, max: %i bp, top 95%%: %i bp" % (averageLen[0], averageLen[1], averageLen[2], averageLen[3]))
if averageLen[0] < int(args.length):
    TruncLen = int(averageLen[3])
    lib.log.error('Warning: Average length of reads %i bp, is less than specified truncation length %s bp' % (averageLen[0], args.length))
    lib.log.error('Resetting truncation length to %i bp (keep > 95%% of data) ' % TruncLen)
else:
    TruncLen = int(args.length)
'''
#now split into individual files
lib.log.info("Splitting FASTQ file by Sample")
filtfolder = tmpdir+'_filtered'
if os.path.isdir(filtfolder):
    shutil.rmtree(filtfolder)
os.makedirs(filtfolder)
splitDemux2(derep, filtfolder)

#Loop through each file and filter for number of sequences
removed = []
for file in os.listdir(filtfolder):
    f = os.path.join(filtfolder, file)
    numrecs = lib.countfastq(f)
    if numrecs < args.min_reads_sample:
        os.remove(f)
        removed.append(file+': '+str(numrecs))
lib.log.info("Removed %i samples with read number < %i (--min_reads_sample)" % (len(removed), args.min_reads_sample))
lib.log.debug("%s\n" % ', '.join(removed))

#now run DADA2 on filtered folder
lib.log.info("Running DADA2 pipeline")
dada2log = os.path.join(tmpdir, args.out+'.dada2.Rscript.log')
dada2out = os.path.join(tmpdir, args.out+'.dada2.csv')
CORES = str(lib.getCPUS())
if args.pool == 'off':
    POOL = 'FALSE'
elif args.pool == 'on':
    POOL = 'TRUE'
with open(dada2log, 'w') as logfile:
    subprocess.call(['Rscript', '--vanilla', dada2script, filtfolder, dada2out, args.platform, POOL, CORES], stdout = logfile, stderr = logfile)

#check for results
if not os.path.isfile(dada2out):
    lib.log.error("DADA2 run failed, please check %s logfile" % dada2log)
    sys.exit(1)
    
#now process the output, pull out fasta, rename, etc
fastaout = os.path.join(tmpdir, args.out+'.otus.tmp')
counter = 1
DADA2tmp = os.path.join(tmpdir, args.out+'.dada2.tmp')
with open(DADA2tmp, 'w') as dada2temp:
    with open(fastaout, 'w') as writefasta:
        with open(dada2out, 'rU') as input:
            #next(input)
            for line in input:
                line = line.replace('\n', '')                    
                line = line.replace('"', '')
                cols = line.split(',')
                if cols[0] == '':
                    cols[0] = "#OTU ID"
                    dada2temp.write('%s\n' % '\t'.join(cols))
                    continue
                Seq = cols[0]
                ID = 'iSeq'+str(counter)
                writefasta.write(">%s\n%s\n" % (ID, Seq))
                dada2temp.write('%s\t%s\n' % (ID, '\t'.join(cols[1:])))
                counter += 1

#get number of bimeras from logfile
with open(dada2log, 'rU') as bimeracheck:
    for line in bimeracheck:
        if line.startswith('Identified '):
            bimeraline = line.split(' ')
            bimeras = int(bimeraline[1])
            totalSeqs = int(bimeraline[5])
        if line.startswith('[1] "dada2'):
            dada2version = line.split(' ')[-1].replace('"\n', '').rstrip()
        if line.startswith('[1] "R '):
            Rversion = line.split(' ')[-1].replace('"\n', '').rstrip()
validSeqs = totalSeqs - bimeras
lib.log.info("R v%s, DADA2 v%s" % (Rversion, dada2version))
lib.log.info('{0:,}'.format(totalSeqs) + ' total inferred sequences (iSeqs)')
lib.log.info('{0:,}'.format(bimeras) + ' denovo chimeras removed')
lib.log.info('{0:,}'.format(validSeqs) + ' valid iSeqs')

#optional UCHIME Ref
uchime_out = os.path.join(tmpdir, args.out+'.nonchimeras.fa')
chimeraFreeTable = os.path.join(tmpdir, args.out+'.otu_table.txt')
iSeqs = os.path.join(tmpdir, args.out+'.iSeqs.fa')
if not args.uchime_ref:
    os.rename(fastaout, iSeqs)
else:
    #check if file is present, remove from previous run if it is.
    if os.path.isfile(uchime_out):
        lib.removefile(uchime_out)
    if os.path.isfile(args.uchime_ref):
        uchime_db = os.path.abspath(args.uchime_ref)
    else:
        lib.log.error("%s is not a valid file, skipping reference chimera filtering" % args.uchime_ref)
        uchime_out = fastaout
    #now run chimera filtering if all checks out
    if not os.path.isfile(uchime_out):
        lib.log.info("Chimera Filtering (VSEARCH) using %s DB" % args.uchime_ref)
        cmd = ['vsearch', '--mindiv', '1.0', '--uchime_ref', fastaout, '--db', uchime_db, '--nonchimeras', uchime_out]
        lib.runSubprocess(cmd, lib.log)
        total = lib.countfasta(uchime_out)
        uchime_chimeras = validSeqs - total
        lib.log.info('{0:,}'.format(total) + ' iSeqs passed, ' + '{0:,}'.format(uchime_chimeras) + ' ref chimeras removed')

    #now reformat OTUs and OTU table, dropping chimeric OTUs from table, sorting the output as well
    nonchimeras = lib.fasta2list(uchime_out)
    inferredSeqs = SeqIO.index(uchime_out, 'fasta')
    with open(iSeqs, 'w') as iSeqout:
        for x in natsorted(nonchimeras):
            SeqIO.write(inferredSeqs[x], iSeqout, 'fasta')
    if not args.debug:
        #clean up chimeras fasta
        lib.removefile(uchime_out)
        if os.path.isfile(fastaout):
            lib.removefile(fastaout)

#filter the sequences based on HMM model
lib.log.info("Translating to Amino Acid Sequence in 3 frames")
trans_out = os.path.join(tmpdir, args.out+'.proteins.fa')
with open(trans_out, 'w') as trans_file:
    with open(iSeqs, 'rU') as trans_in:
        for rec in SeqIO.parse(trans_in, "fasta"):
            Seq1 = rec.seq.translate(to_stop=True)
            if len(Seq1) >= args.min_protlen:
                trans_file.write(">%s;frame_1\n%s\n" % (rec.id, Seq1))
            Seq2 = rec.seq[1:].translate(to_stop=True)
            if len(Seq2) >= args.min_protlen:
                trans_file.write(">%s;frame_2\n%s\n" % (rec.id, Seq2))
            Seq3 = rec.seq[2:].translate(to_stop=True)
            if len(Seq3) >= args.min_protlen:
                trans_file.write(">%s;frame_3\n%s\n" % (rec.id, Seq3))
lib.log.info('{0:,}'.format(lib.countfasta(trans_out))+ ' total translated sequences > %i aa' % args.min_protlen)

#HMM against translated amino acids
trans_hmm = os.path.join(tmpdir, args.out + '.proteins.hmm.txt')
hmm_prot = os.path.join(currentdir, 'lib', 'MHC.hmm')
lib.log.info("Running HMMER3 MHC HMM model")
subprocess.call(['hmmscan', '--cpu', CORES, '--domtblout', trans_hmm, hmm_prot, trans_out], stdout = FNULL, stderr = FNULL)

#now filter results for best hit, and get alignment coordinates
hit_list = {}
with open(trans_hmm, 'rU') as hmmer_prots:
    for qresult in SearchIO.parse(hmmer_prots, "hmmscan3-domtab"):
        hits = qresult.hits
        if len(hits) > 0:
            hitLength = int(hits[0].hsps[0].query_end) - int(hits[0].hsps[0].query_start)
            if hitLength >= args.min_protlen:
                hit_list['%s' % hits[0].query_id] = []
                hit_list['%s' % hits[0].query_id].append(hits[0].hsps[0].query_start)
                hit_list['%s' % hits[0].query_id].append(hits[0].hsps[0].query_end)

#now filter the fasta file to get only hits that have MHC domain and are longer than minimum protein length
pass_prot = os.path.join(tmpdir, args.out + '.proteins.pass.fa')
SeqRecords = SeqIO.to_dict(SeqIO.parse(trans_out, "fasta"))
with open(pass_prot, 'w') as pass_file:
    for key, value in natsorted(hit_list.items()):
        subseq = SeqRecords[key][value[0]:value[1]].seq
        pass_file.write(">%s;%s-%s\n%s\n" % (key, value[0], value[1], subseq))

#finally concatenate duplicated sequences
concat_out = os.path.join(tmpdir, args.out +'.proteins.unique.fa')
keep = {}
total_count = 0
concat_file = open(concat_out, 'wb')
with open(concat_out, 'w') as concat_file:
    for rec in SeqIO.parse(pass_prot, "fasta"):
        total_count += 1
        sequence=str(rec.seq).upper()
        if sequence not in keep:
            keep[sequence]=rec.id
        else:
            keep[sequence]+="|"+rec.id
    flipKeep = {y:x for x,y in keep.iteritems()}
    unique_count = 0
    for key, value in natsorted(flipKeep.items()):
        unique_count += 1
        concat_file.write(">%s\n%s\n" % (key, value))

lib.log.info('{0:,}'.format(total_count)+ ' pass HMM filter, '+'{0:,}'.format(unique_count)+ ' are unique AA sequence.')

#now build OTU table from only those seqs that match the HMM model to some degree
#get keeper list for DNA sequences
keepers = []
for k,v in hit_list.iteritems():
    hit = k.split(';')[0]
    keepers.append(hit)

#get only keepers at DNA level to map
keep_out = os.path.join(tmpdir, args.out + '.hmm.filtered.fa')
with open(keep_out, 'w') as output:
    with open(iSeqs, 'rU') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            if rec.id in keepers:
                SeqIO.write(rec, output, 'fasta')
'''
#I originally wanted to map like I do for other amplicons, but this gets tricky with MHC
#because we want to enforce the 1 bp differences, so then what do you map at, 100%???
#for that reason I think safer to use DADA2 output directly
#setup output files
dadademux = os.path.join(tmpdir, args.out+'.dada2.map.uc')
demuxtmp = os.path.join(tmpdir, args.out+'.original.fa')
uctmp = os.path.join(tmpdir, args.out+'.map.uc')

#map reads to DADA2 filtered iSeqs
lib.log.info("Mapping reads to DADA2 filtered iSeqs")
cmd = ['vsearch', '--fastq_filter', os.path.abspath(SeqIn),'--fastq_qmax', '55', '--fastaout', demuxtmp]
lib.runSubprocess(cmd, lib.log)
cmd = ['vsearch', '--usearch_global', demuxtmp, '--db', keep_out, '--id', '0.97', '--uc', dadademux, '--strand', 'plus', '--otutabout', chimeraFreeTable ]
lib.runSubprocess(cmd, lib.log)
total = lib.line_count(dadademux)
lib.log.info('{0:,}'.format(total) + ' reads mapped to iSeqs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))
'''

#load OTU table into pandas to do some filtering
#output reads per Barcode for original, filtered, and hmm passed files
#now loop through data and find barcoded samples, counting each.....
BarcodeCountA = {}
with open(SeqIn, 'rU') as input:
    header = itertools.islice(input, 0, None, 4)
    for line in header:
        ID = line.split("=")[-1].split(";")[0]
        if ID not in BarcodeCountA:
            BarcodeCountA[ID] = 1
        else:
            BarcodeCountA[ID] += 1
BarcodeCountB = {}
with open(derep, 'rU') as input:
    header = itertools.islice(input, 0, None, 4)
    for line in header:
        ID = line.split("=")[-1].split(";")[0]
        if ID not in BarcodeCountB:
            BarcodeCountB[ID] = 1
        else:
            BarcodeCountB[ID] += 1

bc_count = args.out + '.barcode.counts.csv'
with open(bc_count, 'w') as output:
    output.write("Barcode,Demux_total,Filter_total\n")               
    for k,v in natsorted(BarcodeCountA.items(), key=lambda (k,v): v, reverse=True):
        bc_name = str(k)
        allCount = str(BarcodeCountA[k])
        filtCount = BarcodeCountB.get(k)
        output.write("%s,%s,%s\n" % (bc_name, allCount, filtCount))

#load in original DADA2 OTU table, we want to use that as basepoint
#now make a normalized OTU table using the count data, likely from the filtered dataset, use pandas it should be easier
df = pd.read_csv(DADA2tmp, sep='\t')
df.set_index('#OTU ID', inplace=True)

#sort the table, and keep only those alleles in keepers
df2 = df.reindex(index=natsorted(keepers))
sorted_headers = natsorted(df.columns)
df = df2.reindex(columns=natsorted(df2.columns))

#setup the filter_alleles option, if > 1 then use numbers, if less than 1 then treat as percentage of total samples
if float(args.filter_allele) >= 1:
    allelefilter = int(args.filter_allele)
elif float(args.filter_allele) > 0:
    #that means between 0 and 1, so treat as percentage
    allelefilter = int(round(float(args.filter_allele)*len(keepers)))
else:
    allelefilter = 1
#get sums of columns, filter to minimum reads per sample from command line
reads_sample = args.out+'.reads.per.sample.csv'
ss = df.sum(axis=0)
fs = ss[ss >= args.min_reads_sample]
fs.to_csv(reads_sample)
keep = fs.index
filtered = pd.DataFrame(df, columns=keep)
filt2 = filtered.loc[(filtered != 0).any(1)]
los = filt2.sum(axis=1)
fotus = los[los > 2] #valid allele must be found by more than 2 reads, i.e. no singletons
keep = fotus.index
filt3 = pd.DataFrame(filt2, index=keep)

#to combat barcode switching, remove counts less than 0.1% for each OTU
cleaned_OTUtab = os.path.join(tmpdir, args.out+'.otu_table.csv')
cleaned = []
for row in filt3.itertuples():
    result = [row[0]]
    total = sum(row[1:])
    sub = total * args.index_bleed
    subround = math.ceil(sub) #round up to nearest integer
    for i in row[1:]:
        if i < subround:
            i = 0
        result.append(i)
    cleaned.append(result)

header = ['#OTU ID']
for i in filt3.columns:
    header.append(i)
filt4 = pd.DataFrame(cleaned, columns=header)
filt4.set_index('#OTU ID', inplace=True)
filt4.to_csv(cleaned_OTUtab)

normalTab = os.path.join(tmpdir, args.out+'.normalized.otu_table.csv')
normalTabFilt = os.path.join(tmpdir, args.out+'.filtered.normalized.otu_table.csv')
normal = filt4.truediv(fs)
normal.to_csv(normalTab)
normal[normal < 0] = 0  
normal.to_csv(normalTabFilt)

#convert to binary, apply allele per individual filter here if passed via argument, after filtering, check for empty samples
normal[normal > 0] = 1
los = normal.sum(axis=1)
fotus = los[los >= allelefilter]
keep = fotus.index
filt5 = pd.DataFrame(normal, index=keep)
filt5.to_csv(args.out+'.otu_table.txt', sep='\t', float_format='%.f')

#get the actual read counts from binary table
merge = {}
for index, row in filt5.iteritems():
	merge[index] = []
	for i in range(0, len(row)):
		if row[i] == 0:
			merge[index].append(row[i])
		else:
			merge[index].append(df[index][row.index[i]])

FiltTable = pd.DataFrame(merge, index=list(filt5.index))
FiltTable.index.name = '#OTU ID'
FiltTable.to_csv(args.out+'.otu_table.counts.txt', sep='\t', float_format='%.f')

#get total number of reads for each OTU
TotalCounts = FiltTable.sum(axis=1)
RelativeCounts = TotalCounts.div(FiltTable.values.sum())
RelativeCounts = RelativeCounts.multiply(100)
RelativeCounts.to_csv(args.out+'.iSeqs.relative_counts.csv', float_format='%.4f%%')

#get final sequences, should be what is left in filt5
finalseqs = args.out+'.alleles.fa'
with open(finalseqs, 'w') as output:
    with open(iSeqs, 'rU') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            if rec.id in list(filt5.index.values):
                SeqIO.write(rec, output, 'fasta')
finalprots = args.out+'.alleles.proteins.fa'
with open(finalprots, 'w') as output:
    with open(pass_prot, 'rU') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            if rec.id.split(';')[0] in list(filt5.index.values):
                rec.id = rec.id.split(';')[0]
                rec.name = ''
                rec.description = ''
                SeqIO.write(rec, output, 'fasta')

alleles = []
for row in filt5.T.itertuples():
    count = 0
    for i in row[1:]:
        if i > 0:
            count +=1
    result = [row[0], count]
    alleles.append(result)
counts = pd.DataFrame(alleles)
counts.columns = ['Sample', 'Num_alleles']
counts.to_csv(args.out+'.counts.per.sample.csv', index=False)
avg_alleles = counts.mean(axis=0)
stddev = counts.std(axis=0)
minimum = counts.min(axis=0)
maximum = counts.max(axis=0)

#make a faux taxonomy file for biom format
FauxTax = args.out+'.taxonomy.txt'
with open(FauxTax, 'w') as outfile:
    outfile.write('#OTUID\ttaxonomy\n')
    with open(args.out+'.iSeqs.relative_counts.csv', 'rU') as inputfile:
        for line in inputfile:
            ID = line.split(',')[0]
            Tax = ['k__Animalia', 'p__Chordata', 'c__', 'o__', 'f__', 'g__', 's__']
            outfile.write('%s\t%s\n' % (ID, ';'.join(Tax)))

#convert to biological matrix format (biom)
outBiom = args.out + '.biom'
tmpBiom = args.out + '.biom.tmp'
if args.mapping_file:
    if lib.which('biom'):
        lib.removefile(outBiom)
        lib.removefile(tmpBiom)
        cmd = ['biom', 'convert', '-i', args.out+'.otu_table.counts.txt', '-o', tmpBiom, '--table-type', "OTU table", '--to-json']
        lib.runSubprocess(cmd, lib.log)
        cmd = ['biom', 'add-metadata', '-i', tmpBiom, '-o', outBiom, '--observation-metadata-fp', FauxTax, '--sc-separated', 'taxonomy', '--output-as-json', '-m', args.mapping_file]
        lib.runSubprocess(cmd, lib.log)
        lib.log.info("BIOM OTU table created: %s" % outBiom)
    else:
        lib.log.info("biom program not installed, install via `pip install biom-format`")

#infer phylogeny, here we will use protein sequences to get relationships
lib.log.info("Generating phylogenetic tree")
tree_out = args.out + '.tree.phy'
aln_out = args.out + '.mafft.aln'
cmd = ['mafft', finalprots]
#lib.runSubprocess2(cmd, lib.log, aln_out)
cmd = ['usearch9', '-cluster_agg', finalprots, '-treeout', tree_out]
lib.runSubprocess(cmd, lib.log)
        
if not args.debug:
    shutil.rmtree(tmpdir)
    shutil.rmtree(filtfolder)
    lib.removefile(dada2out)
    lib.removefile(derep)
    #lib.removefile(demuxtmp)
    #lib.removefile(uctmp)
    #lib.removefile(dadademux)
    lib.removefile(FauxTax)
    lib.removefile(tmpBiom)

#Print location of files to STDOUT
print "-------------------------------------------------------"
print "MHC DADA2 Script has Finished Successfully"
print "-------------------------------------------------------"
print "%i samples\n%i total alleles\n%f +/- %f average alleles per sample\n%i - %i range of alleles" % (len(counts), len(filt5.index), round(float(avg_alleles),2), round(stddev, 2), minimum[1], maximum[1])
print "-------------------------------------------------------"
if args.debug:
    print "Tmp Folder of files: %s" % filtfolder
    print "Intermediate files: %s" % tmpdir
print "MHC Alleles: %s" % finalseqs
print "MHC Protein Alleles : %s" % finalprots
print "Allelic OTU Table: %s" % args.out+'.otu_table.txt'
print "MHC phylogeny : %s" % tree_out
if args.mapping_file:
    print "Allelic BIOM table: %s" % outBiom
print "-------------------------------------------------------"
        