#!/bin/env python

import os,sys
import pickle
import re

if len(sys.argv) < 5:
    print "Usage: %s <sequence ID as in fasta header(e.g UniProt accession)> <MSA file path> <tree file path> <idmap file path>" % sys.argv[0]
    sys.exit(0)

inputID=sys.argv[1]
MSAfile=sys.argv[2]
treefile=sys.argv[3]
idmapfile=sys.argv[4]
intrepid_input_name = 'reformatted_MSA.afa'
dir=inputID+'_intrepid'


def main():

      
    # make the directory
    if not os.path.exists(dir):
        os.mkdir(dir)

    # enter it
    os.chdir(dir)

    # read the original sequences
    h = open(MSAfile)
    intext = h.read()
    h.close()


    
    # read the mapping of SEQ#s to headers
    h = open(idmapfile)
    idmap = pickle.load(h)
    h.close()

    # get the seqID for the inputID and  replace the headers with SEQ#s
    for seqno,header in idmap.items():
        intext = intext.replace('>%s\n' % header, '>%s\n' % seqno)
        if re.search(inputID,header):
            seqID=seqno

    # cut into a dict of header: sequence pairs
    inseqs = dict([i.split('\n',1) for i in intext.strip().lstrip('>').split('\n>')])
    
    # remove newline chars (coded better, one might remove whitespace)
    for key,val in inseqs.items():
        inseqs[key] = val.replace('\n', '')

    # load the example sequence.  IS THIS ALWAYS SEQ1??? (most probably)
    seed_seq = inseqs[seqID]

    # remove columns where the seed is gapped
    keys, seqs = zip(*inseqs.items())
    newcols = []
    for num,col in enumerate(zip(*seqs)):
        if not (seed_seq[num] == '-' or seed_seq[num] == '.'):
            newcols.append(col)
    newseqs = map(lambda x: ''.join(x), zip(*newcols))
    
    # write the new alignment file
    h = open(intrepid_input_name, 'w')
    for head,seq in sorted(zip(keys, newseqs)):
        h.write('>%s\n%s\n' % (head, seq))
    h.close()

    # write the config file
    h = open('config.txt', 'w')
    h.write('msa_file %s\ntree_file %s\nsequence_id %s\n' %(intrepid_input_name,treefile,seqID))
    h.close()

    # run intrepid
    os.system('intrepid.pl config.txt 1> intrepid.out 2> intrepid.err')

    os.chdir('..')


if __name__ == '__main__':
    main()
