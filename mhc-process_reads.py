#!/usr/bin/env python

import sys, os, inspect, argparse, shutil
from Bio import SeqIO
import lib.fasta as fasta
import lib.fastq as fastq
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.progress as progress
import lib.die as die

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

#get script path and barcode file name
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pgm_barcodes = os.path.join(script_path, 'lib', 'pgm_barcodes.fa')

parser=argparse.ArgumentParser(prog='mhc-process_reads.py', usage="%(prog)s [options] file.fastq > out.fastq\n%(prog)s -h for help menu",
    description='''This script processes Ion Torrent PGM data by demultiplexing barcodes and trimming primer sequences.  Only full length sequences are kept, i.e. that have both a forward and reverse primer.  By default, 2 mismatches in primer sequence is allowed.''',
    epilog="""Written by Jon Palmer (2015) palmer.jona@gmail.com, modified from script by Robert Edgar""",
    formatter_class=MyFormatter)

parser.add_argument('fastq', help='FASTQ file')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='AGGGTGCTCCTCACAGCCCTGTG', help='Forward Primer (drbF)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='TGTCCCCGCRGCRMATTTCCTG', help='Reverse Primer (drbR)')
parser.add_argument('-b','--list_barcodes', dest="barcodes", default='all', help='Enter Barcodes used separated by commas')
parser.add_argument('-n','--name_prefix', dest="prefix", default='R_', help='Prefix for renaming reads')
parser.add_argument('-m','--min_len', default='50', help='Minimum read length to keep')
parser.add_argument('--rev_comp', action='store_true', help='Reverse complement reads')
parser.add_argument('--mismatches', default=2, type=int, help='Num of primer mismatches (integer)')
parser.add_argument('--mult_samples', dest="multi", default='False', help='Combine multiple samples (i.e. FACE1)')
args=parser.parse_args()

MAX_PRIMER_MISMATCHES = args.mismatches

FileName = args.fastq
FwdPrimer = args.F_primer
RevPrimer = args.R_primer
LabelPrefix = args.prefix
SampleLabel = args.multi
MinLen = int(args.min_len)


#get base name of input file
base = args.fastq.split(".")
base = base[0]

#get barcode list
barcode_file = base + ".barcodes_used.fa"
if os.path.isfile(barcode_file):
    os.remove(barcode_file)
if args.barcodes == "all":
    shutil.copyfile(pgm_barcodes, barcode_file)
else:
    bc_list = args.barcodes.split(",")
    inputSeqFile = open(pgm_barcodes, "rU")
    SeqRecords = SeqIO.to_dict(SeqIO.parse(inputSeqFile, "fasta"))
    for rec in bc_list:
        name = "BC_" + rec
        seq = SeqRecords[name].seq
        outputSeqFile = open(barcode_file, "a")
        outputSeqFile.write(">%s\n%s\n" % (name, seq))
    outputSeqFile.close()
    inputSeqFile.close()
BarcodeFileName = barcode_file
RevPrimer = revcomp_lib.RevComp(RevPrimer)

print >> sys.stderr, "Foward primer: ", FwdPrimer
print >> sys.stderr, "Rev comp'd rev primer ", RevPrimer

SeqCount = 0
OutCount = 0
BarcodeMismatchCount = 0
FwdPrimerMismatchCount = 0
RevPrimerStrippedCount = 0
TooShortCount = 0
LongRead = 0
PadCount = 0

PL = len(FwdPrimer)

Barcodes = fasta.ReadSeqsDict(BarcodeFileName)

def MatchesPrimer(Seq, Primer):
    return primer.MatchPrefix(Seq, Primer)

def FindBarcode(Seq):
    global Barcodes
    for BarcodeLabel in Barcodes.keys():
        Barcode = Barcodes[BarcodeLabel]
        if Seq.startswith(Barcode):
            return Barcode, BarcodeLabel
    return "", ""
    
def OnRec(Label, Seq, Qual):
    global PL, LabelPrefix, Barcode, SeqCount, OutCount, TooShortCount, LongRead, PadCount
    global BarcodeMismatchCount, FwdPrimerMismatchCount, RevPrimerStrippedCount
    global FwdPrimer, RevPrimer

    if SeqCount == 0:
        progress.InitFile(fastq.File)

    progress.File("%u reads, %u outupt, %u bad barcode, %u bad fwd primer, %u rev primer stripped, %u too short." % \
      (SeqCount, OutCount, BarcodeMismatchCount, FwdPrimerMismatchCount, RevPrimerStrippedCount, TooShortCount))

    SeqCount += 1
    Barcode, BarcodeLabel = FindBarcode(Seq)
    if Barcode == "":
        BarcodeMismatchCount += 1
        return

    BarcodeLength = len(Barcode)
    Seq = Seq[BarcodeLength:]
    Qual = Qual[BarcodeLength:]

    Diffs = MatchesPrimer(Seq, FwdPrimer)
    if Diffs > MAX_PRIMER_MISMATCHES:
        FwdPrimerMismatchCount += 1
        return

    if SampleLabel == "False":
        Label = LabelPrefix + str(OutCount) + ";barcodelabel=" + BarcodeLabel + ";"
    else:
        Label = LabelPrefix + str(OutCount) + ";barcodelabel=" + SampleLabel + "_" + BarcodeLabel + ";"

# Strip fwd primer
    Seq = Seq[PL:]
    Qual = Qual[PL:]

    BestPosRev, BestDiffsRev = primer.BestMatch2(Seq, RevPrimer, MAX_PRIMER_MISMATCHES)
    if BestPosRev > 0:
        OutCount += 1
        # Strip rev primer
        RevPrimerStrippedCount += 1
        StrippedSeq = Seq[:BestPosRev]
        StrippedQual = Qual[:BestPosRev]

        # correctness checks
        if 1:
            Tail = Seq[BestPosRev:]
            Diffs2 = primer.MatchPrefix(Tail, RevPrimer)
            if Diffs2 != BestDiffsRev:
                print >> sys.stderr
                print >> sys.stderr, " Seq=" + Seq
                print >> sys.stderr, "Tail=" + Tail
                print >> sys.stderr, "RevP=" + RevPrimer
                die.Die("BestPosRev %u Diffs2 %u BestDiffsRev %u" % (BestPosRev, Diffs2, BestDiffsRev))
            assert StrippedSeq + Tail == Seq
        if args.rev_comp:
            Seq = revcomp_lib.RevComp(StrippedSeq)
            Qual = revcomp_lib.RevComp(StrippedQual)
        else:
            Seq = StrippedSeq
            Qual = StrippedQual
        L = len(Seq)
        if L < MinLen:
            TooShortCount += 1
            OutCount -= 1
            return
        if L > LongRead:
            LongRead = L
        fastq.WriteRec(sys.stdout, Label, Seq, Qual)

fastq.ReadRecs(FileName, OnRec)
progress.FileDone("%u reads, %u outupt, %u bad barcode, %u bad fwd primer, %u rev primer stripped, %u too short" % \
      (SeqCount, OutCount, BarcodeMismatchCount, FwdPrimerMismatchCount, RevPrimerStrippedCount, TooShortCount))
print >> sys.stderr, "%10u seqs" % SeqCount
print >> sys.stderr, "%10u barcode mismatches" % BarcodeMismatchCount
print >> sys.stderr, "%10u fwd primer mismatches (%.1f%% discarded)" % (FwdPrimerMismatchCount, FwdPrimerMismatchCount*100.0/SeqCount)
print >> sys.stderr, "%10u rev primer stripped (%.2f%% kept)" % (RevPrimerStrippedCount, RevPrimerStrippedCount*100.0/SeqCount)
print >> sys.stderr, "%10u too short (%.2f%%)" % (TooShortCount, TooShortCount*100.0/SeqCount)
print >> sys.stderr, "%10u output (%.1f%%)" % (OutCount, OutCount*100.0/SeqCount)
print >> sys.stderr, "%10u Longest read" % (LongRead)            