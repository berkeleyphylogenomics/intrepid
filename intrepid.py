#!/usr/bin/python

import os
import sys
import cPickle
from Bio import SeqIO
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML
from optparse import OptionParser
from datetime import date
import time
import psycopg2
import psycopg2.extras

### db connection constants ###

dbname = 'pfacts003_test'
user = 'webuser'
password = 'w3zx19Ko'

### end db connection constants  ###

"""
this function just formats a time interval into minutes and seconds,
so we can display the elapsed time for users
"""
def getTimeStr(starttime, endtime):
    elapsed_secs = endtime - starttime

    if elapsed_secs > 60:
        mins = int(elapsed_secs / 60)
        elapsed_secs = elapsed_secs - (mins * 60)

        if mins == 1:
            return '[1 min, %.2f secs]' % elapsed_secs

        return '[%d mins, %.2f secs]' % (mins, elapsed_secs)

    return '[%.2f secs]' % elapsed_secs
# end getTimeStr

beginningtime = time.time()

#
# stage 0: read user input, just a seed sequence filename
#

if len(sys.argv) < 2:
    print 'usage: intrepid.py seedfile'
    print '  where seedfile is a fasta-format sequence file with one sequence'
    sys.exit(0)

seedfname = sys.argv[1]

handle = open(seedfname, 'r')
myseedrecord = SeqIO.read(handle, 'fasta')
handle.close()

myseedid = myseedrecord.description

print 'starting INTREPID functional residue prediction for "%s".' % myseedid
print ''

#
# stage 1: blast seed sequence against uniprot using psiblast
#

starttime = time.time()

print 'retrieving homologous sequences from UniProt using PSI-BLAST...',
sys.stdout.flush()

## execute blastpgp ##

blastexe = '/usr/bin/blastpgp'
blastdb  = '/clusterfs/ohana/external/UniProt/current/protein'
iters = 4
eval = 0.0001
maxseqs = 1000

results,errors = NCBIStandalone.blastpgp(blastexe, blastdb, seedfname, expectation=eval, alignments=maxseqs, npasses=iters)

## parse psiblast hits to get ids ##

blasthits = set([])

for blast_record in NCBIXML.parse(results):
    for alignment in blast_record.alignments:
        blasthits.add(alignment.hit_id)

if len(blasthits) < 3:
    print 'Sorry, only %d homologs were retrieved from UniProt, too few sequences to determine patterns of evolutionary conservation.' % len(blasthits)
    sys.exit(0)

handle = open('intrepid-psiblast-ids.txt', 'w')

for x in blasthits:
    print >>handle, x

handle.close()

## add query sequence to sequence list ##

cmd = 'cat %s > intrepid-psiblast-seqs.fa' % seedfname
os.system(cmd)

## retrieve full-length sequences from UniProt ##

cmd = 'fastacmd -d %s -i intrepid-psiblast-ids.txt >> intrepid-psiblast-seqs.fa' % blastdb
os.system(cmd)

endtime = time.time()

print 'done. %s, retrieved %d homologs.' % (getTimeStr(starttime, endtime), len(blasthits))

#
# stage 2: align query and all homologs using muscle
#

starttime = time.time()

print 'aligning query sequence and all homologs...',
sys.stdout.flush()

## align using muscle ##

cmd = 'muscle -maxiters 2 -stable -in intrepid-psiblast-seqs.fa -out intrepid-muscle.align &> /dev/null'
os.system(cmd)

## read alignment ##

query_header = ''
query_sequence = ''
orig_query_sequence = ''

headers = []
sequences = []
orig_sequences = []

handle = open('intrepid-muscle.align', 'r')

for seqrec in SeqIO.parse(handle, 'fasta'):
    if seqrec.description == myseedid:
        query_header = seqrec.description
        orig_query_sequence = seqrec.seq.tostring()

    else:
        headers.append(seqrec.description)
        orig_sequences.append(seqrec.seq.tostring())
        sequences.append('')

handle.close()

## parse alignment 1: remove columns with gaps in query ##

for i in range(len(orig_query_sequence)):
    if orig_query_sequence[i] != '-':
        query_sequence += orig_query_sequence[i]

        for j in range(len(orig_sequences)):
            sequences[j] += orig_sequences[j][i]

## parse alignment 2: remove sequences consisting of 100% gaps ##

kept_headers = []
kept_sequences = []

for i in range(len(headers)):
    if len(sequences[i].replace('-', '')) > 0:
        kept_headers.append(headers[i])
        kept_sequences.append(sequences[i])

## print alignment to output file ##

handle = open('intrepid-alignment.fa', 'w')

print >>handle, '>%s' % query_header
print >>handle, query_sequence

for i in range(len(kept_headers)):
    print >>handle, '>%s' % kept_headers[i]
    print >>handle, kept_sequences[i]

handle.close()

endtime = time.time()

print 'done. %s, %d sequences remaining after removing sequences with only gaps in the alignment.' % (getTimeStr(starttime, endtime), len(kept_headers)+1)

#
# stage 3: construct a phylofacts protein family,
#          including the tree(s) and identification of homologous PDBs
# 

starttime = time.time()

print ''
print 'constructing PhyloFacts book...'
print ''
sys.stdout.flush()

## run buildFamily ##

cmd = 'buildFamily.py --noml --nonjboot intrepid-alignment.fa'
os.system(cmd)

## print notes file ##

handle = open('intrepid-alignment.buildnotes', 'w')

print >>handle, 'gathering_method PSI-BLAST'
print >>handle, 'build_clustering_method PSI-BLAST'
print >>handle, 'build_database_source UniProt'
print >>handle, 'build_alignment_notes Removed columns with gaps in seed'
print >>handle, 'family_specific_evalue_criterion 1.0e-4'
print >>handle, 'family_specific_sw_method local-local'
print >>handle, 'status normal'
print >>handle, 'author bpg'
print >>handle, 'notes intrepid build process'

handle.close()

## run inputFamily ##
'''
cmd = 'inputFamily.py %s intrepid-alignment.fa' % seedfname
os.system(cmd)

## run find_orthologs ##

# get family "book" id from inputFamily
handle = open('intrepid-alignment.familyid', 'r')
bookid = handle.readline().strip()
handle.close()

cmd = 'find_orthologs --book %s --tree nj > find_orthologs.nj.log' % bookid
os.system(cmd)

cmd = 'find_orthologs --book %s --tree sciphy > find_orthologs.sciphy.log' % bookid
os.system(cmd)
'''
endtime = time.time()

sys.stdout.flush()
print ''
print 'done constructing PhyloFacts book. %s' % getTimeStr(starttime, endtime)
print ''

#
# stage 4: run INTREPID to get scores
#

## write special INTREPID-okay alignment   ##
## replacing ambiguity codes with gaps     ##

starttime = time.time()

print 'running INTREPID algorithm for functional residue prediction...',
sys.stdout.flush()

myseqs = {}

handle = open('intrepid-alignment.aln.sciphy', 'r')

for seqrec in SeqIO.parse(handle, 'fasta'):
    myid = seqrec.description
    myseq = seqrec.seq.tostring()

    myseq = myseq.replace('X', '-')
    myseq = myseq.replace('J', '-')
    myseq = myseq.replace('O', '-')
    myseq = myseq.replace('U', '-')
    myseq = myseq.replace('B', '-')
    myseq = myseq.replace('Z', '-')

    myseqs[myid] = myseq

handle.close()

handle = open('intrepid-alignment.aln.intrepid', 'w')

print >>handle, '>SEQ1'
print >>handle, myseqs['SEQ1']

for k in myseqs.keys():
    if k != 'SEQ1':
        print >>handle, '>%s' % k
        print >>handle, myseqs[k]

handle.close()

## run intrepid ##

handle = open('intrepid_config.txt', 'w')

print >>handle, 'msa_file\tintrepid-alignment.aln.intrepid'
print >>handle, 'tree_file\tintrepid-alignment.nj.rooted.tre'
print >>handle, 'sequence_id\tSEQ1'

handle.close()

#run intrepid

cmd = 'intrepid.pl intrepid_config.txt &> intrepid.app.log'
os.system(cmd)
os.system('mv output.aux intrepid.output.aux')

## put scores in db ##

# connect to db
conn = psycopg2.connect("dbname='%s' user='%s' host='db' password='%s'" % (dbname, user, password))
cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

# parse intrepid output file, putting scores in db

handle = open('intrepid.output.aux', 'r')

line = handle.readline()
line = handle.readline()

residue_id = 1

while line:
    myscore = line.split('|')[6]

    sql = 'INSERT INTO residue_score (residue_number, score, method) VALUES (%d, %s, \'intrepid\')' % (residue_id, myscore)
    cur.execute(sql)

    residue_id += 1
    line = handle.readline()

handle.close()

# commit db transaction
conn.commit()

## print that we're done ##

endtime = time.time()

print 'done. %s' % getTimeStr(starttime, endtime)

endingtime = time.time()

print ''
print 'INTREPID finished. %s' % getTimeStr(beginningtime, endingtime)

