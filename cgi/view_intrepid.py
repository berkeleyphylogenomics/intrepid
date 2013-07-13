#!/usr/bin/python

import sequence_graph
import cPickle
import os
import sys
import cgi
import cgitb
cgitb.enable()



# install directory of intrepid webserver part
web_intrepidhome = '/intrepid_alpha'
os_intrepidhome = '/var/www/html' + web_intrepidhome



def printStructure(structdir, myscores, cutoff):

    # get web path to structure file
    webstructfile = structdir + '/structure.pdb'

    # print jmol script

    # launch jmol applet
    print '<p class="bpg_text" style="font-size: 12pt;">'
    print '  Evolutionary conservation plotted on 3D structure'
    print '  <a style="font-size: x-small;" href="%s/help/structureview.html" onClick="return popup(this, \'structure viewer\')">more info</a>' % web_intrepidhome
    print '</p>'

    print '<applet name="structviewer" code="JmolApplet" archive="JmolApplet.jar" codebase="/jmol" width="400" height="400" align="center" mayscript="true">'
    print '  <param name="bgcolor" value="black">'
    print '  <param name="progressbar" value="true">'
    print '  <param name="load" value="%s">' % webstructfile
    print '  <param name="script" value="'
    print '    select all; cpk off; wireframe off;'
    print '    cartoon on; color white;'
    
    for i in range(len(myscores)):
        if myscores[i] >= cutoff:
            print '    select %d; wireframe 100; color blue;' % (i+1)

    print '  zoom 125;'
    print '  ">'
    print '</applet>'
# end printStructure



def printSequenceGraph(sequence, scores, cutoff):
 
    print '<p class="bpg_text" style="font-size: 12pt;">Evolutionary conservation plotted against query sequence'
    print '  <a style="font-size: x-small;" href="%s/help/sequenceview.html" onClick="return popup(this, \'sequence graph\')">more info</a>' % web_intrepidhome
    print '</p>'

    # get sequence header
    seqheader = ''
    display_header = ''
    display_description = ''

    handle = open('id_map.txt', 'r')

    for line in handle:
        line = line.strip()
        id = line.split('>')[0]
        header = line.split('>')[1]
        
        if id == 'QURY_SEQ':
            seqheader = header
            break

    if len(seqheader) <= 17:
        display_header = seqheader
        display_description = seqheader

    else:
        display_header = seqheader[0:17] + '...'
        display_description = seqheader

    # parse intrepid results

    init_min = float(sys.maxint)
    init_max = float(-1.0 * sys.maxint)

    sg = sequence_graph.SequenceGraph(init_min, init_max)

    for v in scores:
        sg.addValue(v)

    sg.addSequence(display_header, display_description, '', sequence)

    graphstr = sequence_graph.buildGraph(sg, cutoff, 60)
    print graphstr
# end printSequenceGraph


def printSequence(sequence):

    myseq = sequence
    myseq = myseq.replace('\r\n', '\n')

    myheader = myseq.split('\n')[0]

    myseq = myseq.replace(myheader, '')
                    
    myseq = myseq.replace('\n', '')
    myseq = myseq.replace('\r', '')
    myseq = myseq.replace('\t', '')
    myseq = myseq.replace(' ', '')    
    printseq = myseq[0:25] + '... ...' + myseq[len(myseq)-25:len(myseq)]

    print '<p class="bpg_text">INTREPID results for:</p>'
    print '<p style="font-size: small; font-family: monaco, courier, monospace;">%s<br>%s</p>' % (myheader, printseq)
# end printSequence



# print header information
print 'Content-type: text/html'
print ''
print '<html>'
print '<head>'
print '  <title>INTREPID Inference of Evolutionarily Conserved Residues</title>'
#print '  <link rel="stylesheet" type="text/css" href="/intrepid/ext-2.2/resources/css/ext-all.css"></link>'
print '  <link rel="stylesheet" type="text/css" href="%s/Ext-sliders/common/ext-ux-slidezone.css"></link>' % web_intrepidhome
print '  <link rel="stylesheet" type="text/css" href="%s/intrepid.css"></link>' % web_intrepidhome

print '  <script type="text/javascript" src="%s/ext-2.2/adapter/prototype/prototype.js"></script>' % web_intrepidhome
print '  <script type="text/javascript" src="%s/ext-2.2/adapter/ext/ext-base.js"></script>' % web_intrepidhome
print '  <script type="text/javascript" src="%s/ext-2.2/ext-all.js"></script>' % web_intrepidhome
print '  <script type="text/javascript" src="%s/Ext-sliders/ext2/Ext.ux.SlideZone.js"></script>' % web_intrepidhome
print '</head>'
print '<body>'

print """
        <script type="text/javascript" language="JavaScript">

            Ext.onReady(function() {

                Slider = {}

                var myscorecutoff = 2.0;
                var query = window.location.search.substring(1);
                var parms = query.split('&');

                for(var i=0; i<parms.length; i++) {
                    var pos = parms[i].indexOf('=');
                    
                    if(pos > 0) {
                        var key = parms[i].substring(0,pos);
                        
                        if(key == 'scorecutoff') {
                            myscorecutoff = parseFloat(parms[i].substring(pos+1));
                            break;
                        }
                    }
                }



                Slider.scorecutoffslider = new Ext.ux.SlideZone('scorecutoffslider', {
                    type: 'horizontal',
                    size: 400,
                    sliderWidth: 20,
                    sliderHeight: 20,
                    maxValue: 4.0,
                    minValue: 0.0,
                    sliderSnap: 0.01,
                    sliders: [{ value: myscorecutoff,
                                name: 'scorecutoffval'
                              }]
                });

                Slider.scorecutoffslider.getSlider('scorecutoffval').on('drag',
                    function() {
                        var theval = parseFloat(Slider.scorecutoffslider.getSlider('scorecutoffval').value);
                        var rounded = Math.round(theval * 100) / 100;
                        $('scorecutoffslider_val').innerHTML = rounded;
                    });

                Slider.scorecutoffslider.getSlider('scorecutoffval').on('dragend',
                    function() {
                        var theval = parseFloat(this.value);
                        var rounded = Math.round(theval * 100) / 100;

                        document.resubmitform.scorecutoff.value = rounded;
                        document.resubmitform.submit();
                    });

                var theval = parseFloat(Slider.scorecutoffslider.getSlider('scorecutoffval').value);
                var rounded = Math.round(theval * 100) / 100;
                $('scorecutoffslider_val').innerHTML = rounded;
            
            });

        </script>
"""

print '  <script type="text/javascript">'
print '    function popup(mylink, windowname)'
print '    {'
print '      if (! window.focus) return true;'
print '      var href;'
print '      if (typeof(mylink) == \'string\')'
print '        href=mylink;'
print '      else'
print '        href=mylink.href;'
print '      window.open(href, windowname, \'width=400,height=450,scrollbars=yes\');'
print '      return false;'
print '    }'
print '  </script>'
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
print '            <td valign=top align=left>'
# end printing header



def printInvalidResults():
    print 'A valid results directory is required to view INTREPID results.'
# end printInvalidResults


def printSites(myseq, myscores, cutoff):
    
    meets_cutoff_count = 0

    for i in range(len(myscores)):
        if myscores[i] >= cutoff:
            meets_cutoff_count += 1

    print '<p class="bpg_text">%d sites meet conservation score cutoff.</p>' % meets_cutoff_count

    if meets_cutoff_count < 1:
        pass

    elif meets_cutoff_count <= 20:
        print '<table>'
        print '<tr><th>position</th><th>residue</th><th>score</th></tr>'

        for i in range(len(myscores)):
            if myscores[i] >= cutoff:
                print '  <tr><td style="font-size: 10pt;">%d</td><td style="font-size: 10pt">%s</td><td style="font-size: 10pt">%.4f</td></tr>' % ((i+1), myseq[i], myscores[i])

        print '</table>'

    # print a 2-col table
    else:
        mynewpos = []
        mynewseq = []
        mynewscores = []

        for i in range(len(myscores)):
            if myscores[i] >= cutoff:
                mynewpos.append(i)
                mynewseq.append(myseq[i])
                mynewscores.append(myscores[i])

        icutoff = (len(mynewpos)/2) + 1

        if len(mynewpos) % 2 == 0:
            icutoff = len(mynewpos) / 2

        print '<table>'
        print '<tr><th>position</th><th>residue</th><th>score</th><th width=60>&nbsp</th><th>position</th><th>residue</th><th>score</th></tr>'

        for i in range(icutoff):
            print '<tr><td style="font-size: 10pt;">%d</td><td style="font-size: 10pt;">%s</td><td style="font-size: 10pt;">%.4f</td><td></td>' % ((mynewpos[i]+1), mynewseq[i], mynewscores[i])

            newi = icutoff + i

            if newi < len(mynewpos):
                print '<td style="font-size: 10pt;">%d</td><td style="font-size: 10pt;">%s</td><td style="font-size: 10pt;">%.4f</td></tr>' % ((mynewpos[newi]+1), mynewseq[newi], mynewscores[newi])

            else:
                print '</tr>'

        print '</table>'
# end printSites



def printResults(workingdir):

    # print resubmit form
    print '<form name="resubmitform" method=GET action="./view_intrepid.py">'
    print '<input type=HIDDEN name="workdir" value="%s">' % workingdir
    print '<input type=HIDDEN name="scorecutoff" value="0.0">'
    print '</form>'

    # get sequence

    handle = open('seed.fa', 'rU')

    line = handle.readline()
    sequence = line
    myseq = ''

    line = handle.readline()

    while line:
        sequence += line

        line = line.strip()
        myseq += line

        line = handle.readline()

    handle.close()

    myseq = myseq.replace(' ', '')
    myseq = myseq.replace('\r', '')
    myseq = myseq.replace('\n', '')
    myseq = myseq.replace('\t', '')

    printSequence(sequence)

    # print results link, if results are there
    theparts = workdir.split('/')
    thedirfinal = theparts[len(theparts)-1]

    resdir = web_intrepidhome + '/results/' + thedirfinal
    osresdir = os_intrepidhome + '/results/' + thedirfinal

    if os.path.exists(osresdir + '/intrepid.results.txt'):
        print '<a style="font-size: 10pt;" href="%s/intrepid.results.txt">get raw INTREPID results</a>' % resdir

    # get cutoff for values to highlight
    cutoff = 2.0

    if form.has_key('scorecutoff'):
        try:
            cutoff = float(form['scorecutoff'].value)
        except:
            cutoff = 2.0

    # parse INTREPID output

    myscores = []

    handle = open('intrepid.output.aux', 'r')

    line = handle.readline()
    line = handle.readline()

    ambig = set(['B', 'Z', 'X', 'J', 'O', 'U', '.', '*', '-', '?'])

    index = 0

    while line:
        line = line.strip()
        lst = line.split('|')

        char = lst[3]
        score = float(lst[6])
        
        if myseq[index] in ambig:
            myscores.append(0.0)

        else:
            myscores.append(score)
        
        line = handle.readline()
        index += 1

    handle.close()

    print '<table cellspacing=20>'
    print '  <tr>'
    print '    <td valign=top>'

    # print structure display, if there is a structure
    workparts = workdir.split('/')
    dirfinal = workparts[len(workparts)-1]

    structdir = web_intrepidhome + '/structures/' + dirfinal
    osstructdir = os_intrepidhome + '/structures/' + dirfinal

    if os.path.exists(osstructdir + '/structure.pdb'):
        printStructure(structdir, myscores, cutoff)

    # print slider controls
    print '<table valign=top><tr><td></td></tr>'
    print '<tr><td><p class="bpg_text">Highlighting residues with conservation score at least:</p></td>'
    print '<td id="scorecutoffslider_val"></td></tr></table>'
    print '    <div id="scorecutoffslider"></div>'
    #print '    </td>'
    
    # print sites meeting cutoff
    print '<tr><td>'
    printSites(myseq, myscores, cutoff)
    print '</td></tr>'

    print '</table>' 
    print '</td>'
    
    # print sequence graph
    print '    <td valign=top style="border-left: thin solid black;">'
    printSequenceGraph(myseq, myscores, cutoff)
    print '    </td>'
    
    print '  </tr>'
    print '</table>'
# end printResults



# get form data
form = cgi.FieldStorage()

# look up existing results
if form.has_key('workdir'):
    workdir = form['workdir'].value

    if os.path.exists(workdir):
        os.chdir(workdir)
        printResults(workdir)

    else:
        printInvalidResults()

# no user-supplied results directory, register an error
else:
    printInvalidResults()



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

