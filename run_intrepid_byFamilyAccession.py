#! /usr/bin/python
'''
This script takes a bpg accession and a uniprot accession, and run interpid to
score the residues in the sequence corresponding to the uniprot accession 
'''
import sys
from bpg.common.utils.dir_of_family import get_dir_of_family_accession as get_dir
import glob
import os

def get_file(family_dir,extention):
  '''get the msa and tree file'''
  file=glob.glob(os.path.join(family_dir,extention))
  if len(file) != 1:
    print "Incorrect number of %s file in %s" % (extention, family_dir)
    exit()
  return file[0]

def main():
  if len(sys.argv)<2:
    print 'Usage: %s "<uniprot accession> <bpg accession>"' % sys.argv[0]
    exit()

  bpg_acc = sys.argv[1].split()[1]
  seq_acc = sys.argv[1].split()[0]
  family_dir = get_dir(bpg_acc)
  msa_file = get_file(family_dir, '*mafft.afa')
  tree_file = get_file(family_dir, '*ml.rooted.tre')
  idmap_file = get_file(family_dir, '*ungapped.idmap')
  cmd = "python "
  cmd+="/home/yshen/ohana_repository/bpg/intrepid/run_intrepid_against_manuscripts.py "
  cmd+="%s %s %s %s" % (seq_acc,msa_file,tree_file,idmap_file)  
  print cmd
  os.system(cmd)

if __name__=='__main__':
  main()
