#!/usr/bin/python

import os
import glob
import cgi
import cgitb
cgitb.enable()



# install directory of intrepid webserver part
web_intrepidhome = '/intrepid_alpha'
os_intrepidhome = '/var/www/html' + web_intrepidhome



def isValid(sequence):

    valid_chars = set([])

    valid_chars.add('A')
    valid_chars.add('R')
    valid_chars.add('N')
    valid_chars.add('D')
    valid_chars.add('C')
    valid_chars.add('E')
    valid_chars.add('Q')
    valid_chars.add('G')
    valid_chars.add('H')
    valid_chars.add('I')
    valid_chars.add('L')
    valid_chars.add('K')
    valid_chars.add('M')
    valid_chars.add('F')
    valid_chars.add('P')
    valid_chars.add('S')
    valid_chars.add('T')
    valid_chars.add('W')
    valid_chars.add('Y')
    valid_chars.add('V')

    valid_chars.add('B')
    valid_chars.add('Z')
    valid_chars.add('X')
    valid_chars.add('J')
    valid_chars.add('O')
    valid_chars.add('U')

# we're not going to allow these:
#    valid_chars.add('.')
#    valid_chars.add('*')
#    valid_chars.add('-')
#    valid_chars.add('?')

    sequence = sequence.strip()
    sequence = sequence.rstrip()

    sequence = sequence.replace('\r\n', '\n')
   
    if sequence == '':
        return ('Sequence appears to be empty.', '', '')

    seqlist = sequence.split('\n')
    
    header = seqlist[0]
    
    if header[0] != '>':
        return ('Header "%s" must start with ">" to be a valid FASTA header.' % header, header, sequence)

    checkheader = header.strip()

    if checkheader == '>':
        return ('Header "%s" is blank; please enter a valid FASTA header.' % header, header, sequence)

    parsed_sequence = ''

    for i in range(1,len(seqlist)):
        parsed_sequence += seqlist[i]

    parsed_sequence = parsed_sequence.replace('\n', '')
    parsed_sequence = parsed_sequence.replace('\t', '')
    parsed_sequence = parsed_sequence.replace(' ', '')
    parsed_sequence = parsed_sequence.replace('\r', '')

    if len(parsed_sequence) < 1:
        return ('Submitted sequence has no characters; sequence must be of non-zero length.', header, parsed_sequence)

    acgtu_count = 0

    for i in range(len(parsed_sequence)):
        thechar = parsed_sequence[i]

        if thechar.islower():
            return ('"%s" (at position %d) is a lower-case characer; characters must be upper-case.' % (thechar, (i+1)), header, parsed_sequence) 

        if thechar not in valid_chars:
            return ('"%s" (at position %d) is not a valid residue code.' % (thechar, (i+1)), header, parsed_sequence)

        if thechar == 'A' or thechar == 'C' or thechar == 'G' or thechar == 'T' or thechar == 'U':
            acgtu_count += 1

    acgtu_proportion = float(acgtu_count) / float(len(parsed_sequence))

    if acgtu_proportion > 0.65:
        return ('This sequence appears to be nucleotides; a protein sequence is required.', header, parsed_sequence)

    return ('', header, parsed_sequence)
# end isValid



def printHeader(print_refresh, working_dir):
    # print header information
    print 'Content-type: text/html'
    print ''
    print '<html>'
    print '<head>'
    print '  <title>INTREPID Inference of Evolutionarily Conserved Residues</title>'
    print '  <link rel="stylesheet" type="text/css" href="%s/intrepid.css"></link>' % web_intrepidhome

    if print_refresh:
        print '  <meta http-equiv="refresh" content="15; url=./run_intrepid.py?workdir=%s"></meta>' % working_dir

    print '</head>'
    print '<body>'
    print '  <table cellspacing=10 border=0 align=center>'
    print '    <tr>'
    print '      <td align=center>'
    print '        <img src="%s/intrepid-banner.jpg"></img>' % web_intrepidhome
    print '      </td>'
    print '    </tr>'
    print '    <tr>'
    print '      <td align=center>'
    print '        <p class="bpg_text">'
    print '        Created by the'
    print '        <a href="http://phylogenomics.berkeley.edu">'
    print '          Berkeley Phylogenomics Group'
    print '        </a>'
    print '        </p>'
    print '      </td>'
    print '    </tr>'
    print '    <tr><td><hr style="color: black; background-color: black; height: 3px;"></td></tr>'
    print '    <tr>'
    print '      <td align=center>'
    print '        <table width=500 border=0>'
    print '          <tr>'
    print '            <td width=60>&nbsp</td>'
    print '            <td align=left>'
    print '              <p class="bpg_text">'

    if print_refresh:
        print '              INTREIPD analysis may take an hour or more; this page will refresh every 15 seconds until results are available.<br><br>'
        print '              You can check results later using this URL: <a href="http://phylogenomics.berkeley.edu/cgi-bin/intrepid/run_intrepid.py?workdir=%s">http://phylogenomics.berkeley.edu/cgi-bin/intrepid/run_intrepid.py?workdir=%s</a><br><br>' % (working_dir, working_dir)
# end printHeader



def getWorkingDir():
    # get working directory
    dir_base = '/home/bpg/intrepid'
    id_file = dir_base + '/next_id.txt'

    existing_dirs = glob.glob('%s/I*' % dir_base)

    mydirnum = len(existing_dirs) + 1
    mydir = dir_base + '/I%d' % mydirnum

    if os.path.exists(mydir):
        return getWorkingDir()

    else:
        os.mkdir(mydir)
        os.system('chmod 777 %s' % mydir)
        return mydir

"""
    dir_bad = True

    while dir_bad:
        handle = open(id_file, 'r')
        my_id = int(handle.readline())
        handle.close()

        handle = open(id_file, 'w')
        next_id = my_id + 1
        handle.write('%d' % next_id)
        handle.close()

        working_dir = dir_base + '/I%d' % my_id

        if os.path.exists(working_dir):
            dir_bad = True

        else:
            os.mkdir(working_dir)
            os.system('chmod 777 %s' % working_dir)
            dir_bad = False

    return working_dir
"""
# end getWorkingDir



def printSequence(sequence):

    myseq = sequence

    myseq = myseq.replace('\r\n', '\n')

    myheader = myseq.split('\n')[0]

    myseq = myseq.replace(myheader, '')

    myseq = myseq.replace('\n', '')
    myseq = myseq.replace('\r', '')
    myseq = myseq.replace('\t', '')
    myseq = myseq.replace(' ', '')

    printseq = myseq

    if len(myseq) > 80:
        printseq = myseq[0:25] + '... ...' + myseq[len(myseq)-25:len(myseq)]

    print 'Executing INTREPID analysis for:</p>'
    print '<p style="font-size: small; font-family: monaco, courier, monospace;">%s<br>%s</p>' % (myheader, printseq)
# end printSequence



def printLog(running):

    log_string = ''

    if running:
        size = 'small'
        log_string = 'INTREPID running...<br>\n'

    else:
        size = 'x-small'

    if os.path.exists('intrepid.log'):
        handle = open('intrepid.log', 'r')

        for line in handle:
            log_string += line

        handle.close()
        
        print '<pre style="font-size: %s; font-family: arial, sans-serif;">%s</pre>' % (size, log_string)

    else:
        print '<pre style="font-size: small; font-family: arial, sans-serif;">INTREPID starting...</pre>'
# end printLog



def printResults(workdir):
    handle = open('intrepid_finished.txt', 'r')

    firstline = handle.readline()
    firstline = firstline.strip()

    print '<p class="bpg_text">%s</p>' % firstline

    if firstline != 'INTREPID finished.':
        print '<p class="bpg_text">INTREPID log:</p>'
        printLog(False)

    else:
        print '<p class="bpg_text">Results are available here: <a href="http://phylogenomics.berkeley.edu/cgi-bin/intrepid/view_intrepid.py?workdir=%s">http://phylogenomics.berkeley.edu/cgi-bin/intrepid/view_intrepid.py?workdir=%s</a></p>' % (workdir, workdir)
        print '<p class="bpg_text">INTREPID log:</p>'
        printLog(False)

    handle.close()
# end printResults



# get form data
form = cgi.FieldStorage()
print_refresh = True

# existing analysis, look it up
if form.has_key('workdir'):
    workdir = form['workdir'].value
    os.chdir(workdir)

    # get sequence

    sequence = ''

    handle = open('seed.fa', 'rU')

    line = handle.readline()

    while line:
        sequence += line
        line = handle.readline()

    handle.close()

    # check if done
    intrepid_finished = False

    if os.path.exists('intrepid_finished.txt'):
        intrepid_finished = True

    if intrepid_finished:
        printHeader(False, workdir)
        printSequence(sequence)
        printResults(workdir)

    else:
        printHeader(True, workdir)
        printSequence(sequence)
        printLog(True)

# new analysis, create it
elif form.has_key('sequence'):
    sequence = form['sequence'].value

    # transform sequence to upper case if necessary
    sequence = sequence.strip()
    sequence = sequence.rstrip()
    sequence = sequence.replace('\r\n', '\n')
    seq_array = sequence.split('\n')

    sequence = seq_array[0] + '\n'

    for i in range(1,len(seq_array)):
        sequence += seq_array[i].upper() + '\n'

    (invalid_seq, myheader, myseq) = isValid(sequence)

    if invalid_seq == '':

        # create new working directory
        working_dir = getWorkingDir()
        os.system('cp exec_intrepid.py %s' % working_dir)
        os.chdir(working_dir)

        # print email address if user supplied one
        if form.has_key('email'):
            email_addr = form['email'].value

            if email_addr != '':
                handle = open('email_addr.txt', 'w')
                print >>handle, email_addr
                handle.close()

        # print seed.fa
        handle = open('seed.fa', 'w')
        handle.write(sequence)
        handle.close()

        # print ID map
        handle = open('id_map.txt', 'w')
        print >>handle, 'QURY_SEQ%s' % myheader
        handle.close()

        # print good_seed.fa
        handle = open('good_seed.fa', 'w')
        print >>handle, '>QURY_SEQ'
        print >>handle, myseq
        handle.close()

        # print and run execution script
        handle = open('exec_intrepid.sh', 'w')

        print >>handle, '#!/bin/sh -l'
        print >>handle, 'source /etc/profile'
        print >>handle, 'export LD_LIBRARY_PATH=/usr/local/modeller9v5/lib/x86_64-intel8:$LD_LIBRARY_PATH'
        print >>handle, 'export DYLD_LIBRARY_PATH=/usr/local/modeller9v5/lib/x86_64-intel8'
        print >>handle, 'export LIB_PATH=/usr/local/modeller9v5/lib/x86_64-intel8'
        print >>handle, 'export PYTHONPATH=/usr/local/modeller9v5/lib/x86_64-intel8:/usr/local/modeller9v5/modlib'
        print >>handle, 'export MODINSTALL9v5=/usr/local/modeller9v5'
        print >>handle, 'cd %s' % working_dir
        print >>handle, '%s/exec_intrepid.py seed.fa' % working_dir

        handle.close()

        os.system('chmod 777 ./exec_intrepid.sh')
        cmd = 'SGE_ROOT=/usr/local/gridengine; export SGE_ROOT; /usr/local/gridengine/bin/lx24-amd64/qsub -cwd -S /bin/sh %s/exec_intrepid.sh &> qsub.submit' % working_dir
        os.system(cmd)

        # print html
        printHeader(True, working_dir)
        printSequence(sequence)
        printLog(True)

    # bad FASTA sequence
    else:
        printHeader(False, '')
        print 'A FASTA-formatted protein sequence is required to run INTREPID.<br>'
        print '<pre style="font-size: small; font-family: monaco, courier, monospace;">%s</pre>' % sequence
        print '<p class="bpg_text">%s</p>' % invalid_seq

# no user-supplied query sequence, register an error
else:
    printHeader(False, '')
    print 'A FASTA-formatted protein sequence is required to run INTREPID.'

# print end of page
print '              </p>'
print '            </td>'
print '            <td width=60>&nbsp</td>'
print '          </tr>'
print '        </table>'
print '      </td>'
print '    </tr>'
print '  </table>'
print '</body>'
print '</html>'

