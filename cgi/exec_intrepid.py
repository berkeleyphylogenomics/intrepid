#!/usr/bin/python

import os
import sys
import re
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.PDB.PDBParser import PDBParser
from operator import itemgetter
from build_model import build_model



# install directory of intrepid webserver part
web_intrepidhome = '/intrepid_alpha'
os_intrepidhome = '/var/www/html' + web_intrepidhome



def sendEmail():

    try:
        # read email address file
        handle = open('email_addr.txt', 'r')
        email_address = handle.readline()
        email_address = email_address.strip()
        handle.close()

        email_subject = 'INTREPID analysis finished'

        email_msg  = 'Results are available at:\n\n'
        email_msg += 'http://phylogenomics.berkeley.edu/cgi-bin/intrepid/run_intrepid.py?workdir=%s\n\n' % os.getcwd()
        email_msg += 'for the following query sequence:\n\n'
        
        handle = open('seed.fa', 'rU')

        for line in handle:
            line = line.strip()
            email_msg += line + '\n'

        handle.close()

        cmd  = "source /etc/profile\n"
        cmd += "qsub_mail.py -s \"%s\" %s << 'EOF_MAIL'\n" % (email_subject, email_address)
        cmd += email_msg + "\n"
        cmd += "EOF_MAIL\n"
        os.system(cmd)
    
    except:
        print 'Email failed.'
#end sendEmail



def dieError(error_str):
    handle = open('intrepid_finished.txt', 'w')
    print >>handle, error_str
    handle.close()

    if os.path.exists('email_addr.txt'):
        sendEmail()

    sys.exit(0)
#end dieError



def parsePsiBlast(psiblastfilename, max_evalue):

    try:
        results_dict = {}

        handle = open(psiblastfilename, 'r')

        for blast_record in NCBIXML.parse(handle):
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= max_evalue:
                        subjid = alignment.title

                        if subjid in results_dict:
                            if hsp.expect < results_dict[subjid]:
                                results_dict[subjid] = hsp.expect

                        else:
                            results_dict[subjid] = hsp.expect

        handle.close()

        return results_dict

    except:
        dieError('ERROR: PSI-BLAST failed.')
# end parsePsiBlast



def reduceRedundancy(prop_identical_cutoff):

    cmd = 'uniqueseq final.parsed.nr -alignfile final.parsed.a2m -percentid %f' % prop_identical_cutoff
    os.system(cmd)
    os.system('mv final.parsed.nr.a2m final.parsed.a2m')

    handle = open('final.parsed.a2m', 'r')
                        
    new_ids = []
    new_seqs = []
                                                
    for record in SeqIO.parse(handle, 'fasta'):
        new_ids.append(record.id)
        new_seqs.append(record.seq.tostring())
                                                                                            
    handle.close()

    return (new_ids, new_seqs)
# end reducRedundancy



def buildMSA():
    try:
        # constuct HMM from query sequence using w0.5
        cmd = 'w0.5 good_seed.fa good_seed.mod'
        os.system(cmd)

        # align query sequence and all homologs to HMM
        cmd = 'align2model final -modelfile good_seed.mod -db good_seed.fa -db good_homologs.fa -sw 2'
        os.system(cmd)

        # now we need to parse the .a2m file

        # read sequences
        orig_seqs = []

        handle = open('final.a2m', 'r')

        for seqrecord in SeqIO.parse(handle, 'fasta'):
            id = seqrecord.description
            seq = seqrecord.seq.tostring()
            orig_seqs.append((id, seq))

        handle.close()

        # parse sequences

        new_seqs = []

        for (i,s) in orig_seqs:
            new_seqs.append('')

        (qid, qseq) = orig_seqs[0]

        for i in range(len(qseq)):
            # query is lower-case, transfer it to upper, and other seqs are gaps
            if qseq[i].islower():
                new_seqs[0] += qseq[i].upper()

                for j in range(1,len(new_seqs)):
                    new_seqs[j] += '-'

            # query is upper-case, transfer it, and transfer other seqs
            if qseq[i].isupper():
                new_seqs[0] += qseq[i]

                for j in range(1,len(new_seqs)):
                    (sid, sseq) = orig_seqs[j]
                    new_seqs[j] += sseq[i]

        # print new MSA in fasta format
        seq_count = 0

        handle = open('final.parsed.a2m', 'w')

        for i in range(len(new_seqs)):
            (id, oldseq) = orig_seqs[i]
            theseq = new_seqs[i]

            # replace all ambiguity codes with gaps for now
            # this should really be handled in intrepid
            theseq = theseq.replace('X', '-')
            theseq = theseq.replace('J', '-')
            theseq = theseq.replace('O', '-')
            theseq = theseq.replace('U', '-')
            theseq = theseq.replace('B', '-')
            theseq = theseq.replace('Z', '-')
            theseq = theseq.replace('?', '-')
            theseq = theseq.replace('.', '-')
            theseq = theseq.replace('*', '-')

            gaps = theseq.count('-')

            if i == 0 or float(gaps) / float(len(theseq)) < 0.5:
                print >>handle, '>%s' % id
                print >>handle, theseq
                seq_count += 1

        handle.close()

        # drop redundant sequences
        cmd = 'uniqueseq final.parsed.nr100 -alignfile final.parsed.a2m -percentid 1.0'
        os.system(cmd)
        os.system('mv final.parsed.nr100.a2m final.parsed.a2m')

        # read in non-redundant sequences
        handle = open('final.parsed.a2m', 'r')

        new_ids = []
        new_seqs = []

        for record in SeqIO.parse(handle, 'fasta'):
            new_ids.append(record.id)
            new_seqs.append(record.seq.tostring())

        handle.close()

        # error if too few sequences
        if len(new_ids) < 4:
            dieError('Only %d sequences remained after removing duplicate and poorly-aligned sequences; too few sequences to run INTREPID analysis.' % len(new_ids))

        prop_id_cutoff = 0.9

        while len(new_ids) > 1000 and prop_id_cutoff > 0.6:
            (my_new_ids, my_new_seqs) = reduceRedundancy(prop_id_cutoff)

            if len(my_new_ids) >= 4:
                new_ids = my_new_ids
                new_seqs = my_new_seqs
        
            else:
                break;

            prop_id_cutoff -= 0.1

        # print new MSA in fasta format
        handle = open('final.parsed.a2m', 'w')

        for i in range(len(new_ids)):
            print >>handle, '>%s' % new_ids[i]
            print >>handle, new_seqs[i]

        handle.close()

        # print new MSA in phylip format
        ntax = len(new_ids)
        nchar = len(new_seqs[0])

        handle = open('final.parsed.aln', 'w')

        print >>handle, '%d %d' % (ntax, nchar)

        for i in range(len(new_seqs)):
            id = new_ids[i]
            theseq = new_seqs[i]

            theline = id

            for k in range(len(theline),11):
                theline += ' '

            print >>handle, '%s %s' % (theline, theseq)

        handle.close()

        return ntax
    except:
        dieError('ERROR: Construction of multiple sequence alignment failed.')
# end buildMSA



def inferTree():
    
    try:
        # calculate distance matrix
        os.system('cp final.parsed.aln infile')
        os.system('echo "y\n" | /usr/lib/phylip/bin/protdist >& protdist.out')
        os.system('mv outfile final.parsed.dist')
        os.system('cp final.parsed.dist infile')

        # build tree
        os.system('echo "y\n" | /usr/lib/phylip/bin/neighbor >& neighbor.out')
        os.system('mv outtree final.parsed.tre')
        os.system('perl -pi -e \'s/\n//\' final.parsed.tre')

        # reroot tree at midpoint
        os.system('mpreroot final.parsed.tre final.parsed.rooted.tre')

    except:
        dieError('ERROR: Phylogenetic tree inference failed.')
# end inferTree



def runIntrepid():

    try:
        # print intrepid config file
        handle = open('intrepid_config.txt', 'w')

        print >>handle, 'msa_file\tfinal.parsed.a2m'
        print >>handle, 'tree_file\tfinal.parsed.rooted.tre'
        print >>handle, 'sequence_id\tQURY_SEQ'

        handle.close()

        # run intrepid
        cmd = 'intrepid.pl intrepid_config.txt > intrepid.screendump'
        os.system(cmd)
        os.system('mv output.aux intrepid.output.aux')
    
        if not os.path.exists('intrepid.output.aux'):
            dieError('ERROR: INTREPID analysis of evolutionary conservation failed.')

    except:
        dieError('ERROR: INTREPID analysis of evolutionary conservation failed.')
# end runIntrepid



def scorePDB(pdb_eval):

    try:
        # BLAST pre-screen: use blastp with large e-value cutoff
        #cmd = 'blastall -p blastp -d pdbaa -i seed.fa -m 7 -v 10000 -b 10000 -e 100 -o pdb.blastres'
        #os.system(cmd)

        # parse blast results to get ids

        # get blast results using fastacmd, to get sequences

        # create new HMM from trimmed alignment using w0.5
        cmd = 'w0.5 final.parsed.a2m final.parsed.mod'
        os.system(cmd)

        # score PDB with HMM
        cmd = 'hmmscore pdb -modelfile final.parsed.mod -db /usr/local/blastdb/pdbaa -sw 2 &> hmmscore.log'
        os.system(cmd)

        # parse results for e-values <= pdb_eval
        pdbpatt = re.compile('gi\|(\d+?)\|pdb\|(.+?)\|(.)\s+\d+\s+.+?\s+.+?\s+(.+?)\s')

        handle = open('pdb.dist', 'r')
        outf1 = open('pdbgis.txt', 'w')
        outf2 = open('pdbhits.txt', 'w')

        hitcount = 0

        for line in handle:
            line = line.strip()
            match = pdbpatt.match(line)

            if match:
                gi = match.group(1)
                id = match.group(2)
                chain = match.group(3)
                evalue = float(match.group(4))

                if evalue <= pdb_eval:
                    hitcount += 1
                    print >>outf1, gi
                    print >>outf2, '%s %s %s %.4e' % (gi, id, chain, evalue)

        outf1.close()
        outf2.close()
        handle.close()

        # get pdb sequences using fastacmd
        if hitcount > 0:
            cmd = 'fastacmd -i pdbgis.txt -d pdbaa > pdbhits.fa'
            os.system(cmd)

            # transform IDs for pdb sequences
            handle = open('pdbhits.fa', 'r')
            idf = open('id_map.txt', 'a')
            outf = open('good_pdbhits.fa', 'w')

            curr_pdb = 1

            for line in handle:
                line = line.strip()
            
                if line[0] == '>':
                    myid = 'PDB_%d' % curr_pdb
                    curr_pdb += 1

                    print >>idf, '%s%s' % (myid, line)
                    print >>outf, '>%s' % myid

                else:
                    print >>outf, line

            handle.close()
            idf.close()
            outf.close()

        return hitcount

    except:
        dieError('ERROR: Scoring PDB failed.')
# end scorePDB



def checkAlignment():

    try:
        contig_gap = re.compile('-+')
        handle = open('structure.aligned.fa', 'r')

        for record in SeqIO.parse(handle, 'fasta'):
            theseq = record.seq.tostring()

            # check that there aren't too many gaps
            numgaps = theseq.count('-')
            propgaps = float(numgaps) / float(len(theseq))

            if propgaps > 0.3:
                return False

            gappy_regions = contig_gap.findall(theseq)
        
            if len(gappy_regions) > 0:
                longest_gap = max(len(r) for r in gappy_regions)

                if longest_gap > 30:
                    return False

        handle.close()
        return True
    
    except:
        return False
# end checkAlignment



def build3DModel():

    pdb_lines = []

    try:
        handle = open('pdbhits.txt', 'r')

        for line in handle:
            line = line.strip()
            pdb_lines.append(line)

        handle.close()

    except:
        print 'build3DModel failed on reading pdbhits.txt'

    for line in pdb_lines:
        try:
            pdbid = line.split(' ')[1]
            pdbchain = line.split(' ')[2]

            lcpdbid = pdbid.lower()

            pdbfile = '/home/bpg/pdb_files/%s.pdb' % lcpdbid
            pdbheader = lcpdbid + pdbchain

            # get PDB sequence from ATOM part

            aamap = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
                     'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I',
                     'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
                     'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

            parser = PDBParser()
            structure = parser.get_structure(lcpdbid, pdbfile)

            chain = structure[0][pdbchain]

            pdbseq = ''

            for residue in chain:
                if residue.resname in aamap:
                    pdbseq += aamap[residue.resname]

            # get query sequence
            handle = open('good_seed.fa', 'r')

            line = handle.readline()
            line = line.strip()
            seed_id = line[1:len(line)]

            line = handle.readline()

            seedseq = ''

            while line and line[0] != '>':
                line = line.strip()
                seedseq += line
                line = handle.readline()
    
            handle.close()

            # print PDB sequence, then query seuquence
            handle = open('structure.unaligned.fa', 'w')

            print >>handle, '>%s' % pdbheader
            print >>handle, pdbseq

            theseq = seedseq
            theseq = theseq.replace('X', '-')
            theseq = theseq.replace('J', '-')
            theseq = theseq.replace('O', '-')
            theseq = theseq.replace('U', '-')
            theseq = theseq.replace('B', '-')
            theseq = theseq.replace('Z', '-')
            theseq = theseq.replace('?', '-')
            theseq = theseq.replace('.', '-')
            theseq = theseq.replace('*', '-')

            print >>handle, '>%s' % seed_id
            print >>handle, theseq

            handle.close()

            # align sequences
            cmd = 'mafft --maxiterate 1000 --genafpair structure.unaligned.fa > structure.aligned.fa'
            os.system(cmd)

            # check alignment
            # should throw error if alignment is poor
            if checkAlignment():

                # build model
                build_model('structure.aligned.fa')

                # move model to structure.pdb
                if os.path.exists('QURY_SEQ.B99990005.pdb'):
                    os.system('cp QURY_SEQ.B99990005.pdb structure.pdb')

                elif os.path.exists('QURY_SEQ.B99990004.pdb'):
                    os.system('cp QURY_SEQ.B99990004.pdb structure.pdb')

                elif os.path.exists('QURY_SEQ.B99990003.pdb'):
                    os.system('cp QURY_SEQ.B99990003.pdb structure.pdb')

                elif os.path.exists('QURY_SEQ.B99990002.pdb'):
                    os.system('cp QURY_SEQ.B99990002.pdb structure.pdb')

                elif os.path.exists('QURY_SEQ.B99990001.pdb'):
                    os.system('cp QURY_SEQ.B99990001.pdb structure.pdb')

                else:
                    raise Exception()

                mypaths = os.getcwd().split('/')
                mypath = mypaths[len(mypaths)-1]

                seqpath = '%s/structures/%s' % (os_intrepidhome, mypath)

                if not os.path.exists(seqpath):
                    os.mkdir(seqpath)

                os.system('cp structure.pdb %s/structure.pdb' % seqpath)
                return True

        except:
            print 'build3DModel failed on building model'
        
    return False
#end build3DModel



logfname = 'intrepid.log'
donefname = 'intrepid_finished.txt'

# perform initial 2 rounds of psi-blast
handle = open(logfname, 'a')
print >>handle, 'Retrieving homologous sequences from UniProt using PSI-BLAST...'
handle.close()

psiblastfilename = 'psiblast_1.xml'
cmd = 'blastpgp -i good_seed.fa -d UniProt/current/protein -m 7 -j 4 -v 10000 -b 10000 -C psiblast_1.chk -Q psiblast_1.mat -o %s' % psiblastfilename
os.system(cmd)

# parse psi-blast results
max_eval = 1.0e-3
homologs = parsePsiBlast(psiblastfilename, max_eval)

handle = open(logfname, 'a')
print >>handle, 'Retrieved %d homologs from UniProt with e-value < %.1e.' % (len(homologs), max_eval)
handle.close()

# bail out if there are too few homologs
if len(homologs) < 4:
    dieError('Too few homologs were retrieved; INTREPID cannot infer evolutionary conservation for this sequence at this time.')

# retrieved enough homologs, continue with INTREPID

# sort homologs by e-value
homolog_evalue_list = homologs.items()
homolog_evalue_list.sort(key=itemgetter(1), reverse=False)

# get list of UniProt IDs and get sequences
handle = open('homologs.ids', 'w')

for (id, eval) in homolog_evalue_list:
    uniprot_id = id.split('|')[1]
    print >>handle, uniprot_id

handle.close()

# use fastacmd to get sequences
cmd = 'fastacmd -i homologs.ids -d UniProt/current/protein > homologs.fa'
os.system(cmd)

# generate new IDs for homologs
current_id = 1

handle = open('homologs.fa', 'r')
outf = open('good_homologs.fa', 'w')
idf = open('id_map.txt', 'a')

for line in handle:
    line = line.strip()
    
    if line[0] == '>':
        newid = 'HOM_%d' % current_id
        current_id += 1
        print >>outf, '>%s' % newid
        print >>idf, '%s%s' % (newid, line)

    else:
        print >>outf, line

handle.close()
outf.close()
idf.close()

# generate MSA for INTREPID
handle = open(logfname, 'a')
print >>handle, ''
print >>handle, 'Generating multiple sequence alignment of query sequence and all homologs...'
handle.close()

remaining_seqs = buildMSA()

# generate phylogeny for INTREPID
handle = open(logfname, 'a')
print >>handle, 'Multiple sequence alignment computed; %d sequences remaining after removing redundancy.' % remaining_seqs
print >>handle, ''
print >>handle, 'Inferring protein family phylogeny...'
handle.close()

inferTree()

# score PDB for homologous structures
handle = open(logfname, 'a')
print >>handle, 'Protein family phylogeny computed.'
print >>handle, ''
print >>handle, 'Retrieving homologous structures from PDB using HMM scoring...'
handle.close()

pdb_eval = 1.0e-3
struct_count = scorePDB(pdb_eval)

# build comparative model if possible
handle = open(logfname, 'a')
print >>handle, 'Retrieved %d homologous structures from PDB with e-value < %.1e.' % (struct_count, pdb_eval)
print >>handle, ''

if struct_count < 1:
    print >>handle, 'No homologous structures found, skipping construction of 3D model of query sequence.'
    handle.close()

else:
    print >>handle, 'Constructing 3D structure of query sequence...'
    handle.close()

    model_built = build3DModel()

    if model_built:
        handle = open(logfname, 'a')
        print >>handle, '3D structure of query sequence computed.'
        handle.close()

    else:
        handle = open(logfname, 'a')
        print >>handle, '3D model of query sequence could not be constructed, continuing without it.'
        handle.close()

# run INTREPID to infer evolutionary conservation
handle = open(logfname, 'a')
print >>handle, ''
print >>handle, 'Inferring evolutionary conservation using INTREPID...'
handle.close()

runIntrepid()

handle = open(logfname, 'a')
print >>handle, 'Evolutionary conservation scores computed.'
handle.close()

# print intrepid raw results
try:
    # read sequence
    handle = open('good_seed.fa', 'r')

    seq = handle.readline()
    seq = handle.readline()
    seq = seq.strip()

    handle.close()

    # read scores
    handle = open('intrepid.output.aux', 'r')

    line = handle.readline()
    line = handle.readline()

    scores = []

    while line:
        line = line.strip()

        thescore = line.split('|')[6]
        scores.append(thescore)

        line = handle.readline()

    handle.close()

    # print sequence and scores
    handle = open('intrepid.results.txt', 'w')

    print >>handle, 'position\tresidue\tscore'

    for i in range(len(scores)):
        print >>handle, '%d\t%s\t%s' % ((i+1), seq[i], scores[i])

    handle.close()

    # copy results file
    mypaths = os.getcwd().split('/')
    mypath = mypaths[len(mypaths)-1]

    seqpath = '%s/results/%s' % (os_intrepidhome, mypath)

    if not os.path.exists(seqpath):
        os.mkdir(seqpath)

    os.system('cp intrepid.results.txt %s/intrepid.results.txt' % seqpath)

except:
    dieError('ERROR: INTREPID analysis of evolutionary conservation failed.')

# print email if user supplied an address
if os.path.exists('email_addr.txt'):
    sendEmail()

# finished, print finished file
handle = open('intrepid_finished.txt', 'w')
print >>handle, 'INTREPID finished.'
handle.close()
sys.exit(0)

