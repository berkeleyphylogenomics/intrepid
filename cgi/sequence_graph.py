#!/usr/bin/python

import sys
import cPickle
import cgi

# prnt errors for debugging
import cgitb
cgitb.enable()
        

class SequenceGraph:
    
    def __init__(self, min_val, max_val):
        self.max_value = max_val
        self.min_value = min_val
        self.values = []
        
        self.seqids = []
        self.seqdescs = []
        self.seqURLs = []
        self.seqs = [] 

    def addValue(self, value):
        self.values.append(value)
        
        if value > self.max_value:
            self.max_value = value
            
        if value < self.min_value:
            self.min_value = value

    def addSequence(self, name, description, url, sequence):
        self.seqids.append(name)
        self.seqdescs.append(description)
        self.seqURLs.append(url)
        self.seqs.append(sequence)
        
    def getColumnCount(self):
        return len(self.values)
        
    def getSeqCount(self):
        return len(self.seqids)
        
    def getValue(self, i):
        return self.values[i]
        
    def getSequenceID(self, i):
        return self.seqids[i]
        
    def getSequenceDesc(self, i):
        return self.seqdescs[i]
        
    def getSequenceURL(self, i):
        return self.seqURLs[i]
        
    def getCharacter(self, seq, i):
        return self.seqs[seq][i]
        
        
        
def buildGraph(sg, min_colored, cols_per_line):

    # calculate proportion of graph that's positive and negative
    # based on the max and min values of the sequence graph object
    
    total_height = 100
    
    if sg.min_value >= 0.0:
        positive_height = total_height
        negative_height = 0
        
    else:
        total_value_range = sg.max_value - sg.min_value
        prop_neg = (-1 * sg.min_value) / total_value_range
        
        negative_height = int(total_height * prop_neg)
        positive_height = total_height - negative_height

    # print the html style definitions for the graph

    html = """
<style type="text/css">

.poscol 
{
    background-color: white;
    vertical-align: bottom;
    height: %d;
}

.negcol
{
    background-color:white;
    vertical-align: top;
    height: %d;   
}

.alignment
{
    text-align: center;
    font-family: monaco, courier, monospace;
}

.axis
{
    background-color: black;
}

.seqlabel
{
    text-align: right;
    font-family: arial, sans-serif;
    font-size: small;   
}

.alnlegend
{
    text-align: left;
    font-family: arial, sans-serif;
    font-size: x-small;   
}

.yaxislabel
{
    text-align: right;
    font-family: arial, sans-serif;
    font-size: x-small;   
}

</style>

""" % (positive_height, negative_height)

    # start html graph (it's the table, dummy)
    html += '    <table cellspacing=0 cellpadding=1>\n'
    html += '\n'

    start = 0

    while start < sg.getColumnCount():
        end = start + cols_per_line
        
        if sg.getColumnCount() < end:
            end = sg.getColumnCount()

        html += printGraphRow(sg, min_colored, start, end)
        html += '<tr><td>&nbsp</td></tr>\n'
        start += cols_per_line
    
    html += '</table>\n'
    return html



def printGraphRow(sg, min_colored, start, end):
    # print positive data plot
    html  = '        <tr> <!-- positive Y-axis -->\n'
    html += '\n'
    html += '            <td class="yaxislabel" style="vertical-align: top;">%.2f</td>\n' % sg.max_value
    html += '            <td class="axis"></td>\n'
    html += '\n' 
    html += '            <!-- data plot -->\n'
    
    for i in range(start,end):
        myval = sg.getValue(i)
        html += '            <td class="poscol">'
        
        if myval > 0.0:
            myprop = myval / sg.max_value
            mycolor = 'gray'

            if myval >= min_colored:
                mycolor = 'blue'
            
            myheight = int((myval / sg.max_value) * 100)
            html += '<div title="%s%d [%.4f]" style="height: %d%%; background-color: %s;">&nbsp</div>' % (sg.getCharacter(0,i), (i+1), myval, myheight, mycolor)

        html += '</td>\n'
        
    html += '\n'
    html += '        </tr>\n'
    html += '\n'
        
    # print X-axis
    html += '        <tr> <!-- X-axis -->\n'
    html += '            <td></td>\n'
    html += '            <td class="axis"></td>\n'       
    html += '            <td colspan=%d class="axis"></td>\n' % sg.getColumnCount()
    html += '        </tr>\n'
    html += '\n'
    
    # print negative data plot
    html += '        <tr> <!-- negative Y-axis -->\n'
    html += '\n'
    html += '            <td class="yaxislabel" style="vertical-align: bottom;">%.2f</td>\n' % sg.min_value
    html += '            <td class="axis"></td>\n'
    html += '\n' 
    html += '            <!-- data plot -->\n'

    for i in range(start,end):
        myval = sg.getValue(i)
        html += '            <td class="negcol">'
        
        if myval < 0.0:
            myprop = myval / sg.min_value
            mycolor = 'gray'

            myheight = int((myval / sg.min_value) * 100)
            html += '<div title="%s%d [%.4f]" style="height: %d%%; background-color: %s;">&nbsp</div>' % (sg.getCharacter(0,i), (i+1), myval, myheight, mycolor)

        html += '</td>\n'
        
    html += '\n'
    html += '        </tr>\n'
    html += '\n'

    # print sequences
    html += '        <!-- X-axis sequences -->\n'
    html += '\n'
    
    html += '        <tr>\n'
    
    html += '            <td></td>\n'
    html += '            <td></td>\n'
    
    if end % 10 == 0:
        upper_limit = end-9
        
    else:
        upper_limit = end
    
    j = start+1
    
    while j < upper_limit:
    
        if j == start+1 or j % 10 == 0:
            
            if j < 10:
                html += '            <td class="alnlegend">%d</td>\n' % j
                j += 1
          
            elif j < 100:
                html += '            <td colspan=2 class="alnlegend">%d</td>\n' % j
                j += 2
                
            elif j < 1000:
                html += '            <td colspan=3 class="alnlegend">%d</td>\n' % j
                j += 3
                
            elif j < 10000:
                html += '            <td colspan=4 class="alnlegend">%d</td>\n' % j
                j += 4
                
            elif j < 100000:
                html += '            <td colspan=5 class="alnlegend">%d</td>\n' % j
                j += 5
                
            elif j < 1000000:
                html += '            <td colspan=6 class="alnlegend">%d</td>\n' % j
                j += 6
                
            else:
                html += '            <td colspan=7 class="alnlegend">%d</td>\n' % j
                j += 7
                     
        else:
            html += '            <td class="alnlegend"><br></td>\n'
            j += 1
    
    if end % 10 == 0:
        html += '            <td colspan=10 class="alnlegend" style="text-align: right;">%d</td>\n' % end
        
    html += '        </tr>\n'
    html += '\n'
    
    html += '        <tr>\n'
    html += '            <td></td>\n'
    html += '            <td></td>\n'
    
    for j in range(start+1,end+1):
        if j == start+1 or j % 10 == 0:
            html += '            <td class="alnlegend" style="text-align: center;">|</td>\n'
            
        elif j % 5 == 0:
            html += '            <td class="alnlegend" style="text-align: center;">+</td>\n'

        else:
            html += '            <td class="alnlegend" style="text-align: center;">-</td>\n'

    html += '        </tr>\n'
    html += '\n'
    
    for i in range(sg.getSeqCount()):
        html += '        <tr>\n'
        
        myurl = sg.getSequenceURL(i)
       
        if start == 0:
            if myurl:
                html += '            <td NOWRAP title="%s" class="seqlabel"><a href="%s" target="_blank">%s</a>&nbsp</td>\n' % (sg.getSequenceDesc(i), myurl, sg.getSequenceID(i))
        
            else:
                html += '            <td NOWRAP title="%s" class="seqlabel">%s&nbsp</td>\n' % (sg.getSequenceDesc(i), sg.getSequenceID(i))
        
        else:
            html += '            <td></td>'

        html += '            <td></td>\n'

        sequence_position = start+1

        for j in range(start,end):
            mychar = sg.getCharacter(i,j)

            if mychar == '.' or mychar == '-':
                html += '            <td class="alignment">%s</td>\n' % mychar

            else:
                html += '            <td title="%s%d [%.4f]" class="alignment">%s</td>\n' % (mychar, sequence_position, sg.getValue(j), sg.getCharacter(i,j))
                sequence_position += 1

        html += '        </tr>\n'
        html += '\n'

    return html

    
                
def main():

    if len(sys.argv) < 2:
        print 'Usage: sequence_graph.py inputfile'
        print '  where inputfile is a formatted text file.'
        print '  Here\'s an example of the formatting:'
        print ''
        print '  2 (for 2 sequences, or 1 for 1 sequence)'
        print '  sequence_id_1'
        print '  a 1-line description of sequence_id_1'
        print '  sequence_id_2'
        print '  a 1-line description of sequence_id_2'
        print '  A C 12.3'
        print '  R H 0.0'
        print '  N -  -1.2'
        print '  D W 55.8'
        print ''
        print '  You get the idea.'
        print '  In this case, sequence_id_1 = ARND'
        return 1

    handle = open(sys.argv[1], 'r')

    init_max = float(-1.0 * sys.maxint)
    init_min = float(sys.maxint)
    sg = SequenceGraph(init_min, init_max)
    
    # fill in SequenceGraph object with data
   
    num_seqs = int(handle.readline().strip())

    if num_seqs == 1:

        seqid = handle.readline().strip()
        seqdesc = handle.readline().strip()
        sequrl = ""
        sequence = ""

        line = handle.readline()

        while line:
            line = line.strip()

            residue = line.split(' ')[0]
            value = line.split(' ')[1]

            sg.addValue(float(value))
            sequence += residue

            line = handle.readline()
    
        handle.close()

        sg.addSequence(seqid, seqdesc, sequrl, sequence)

        seqgraphstring = cPickle.dumps(sg)    
        print buildGraph(seqgraphstring)

    else:
        seqid1 = handle.readline().strip()
        seqdesc1 = handle.readline().strip()

        seqid2 = handle.readline().strip()
        seqdesc2 = handle.readline().strip()

        sequrl1 = ""
        sequrl2 = ""

        sequence1 = ""
        sequence2 = ""

        line = handle.readline()

        while line:
            line = line.strip()

            res1 = line.split(' ')[0]
            res2 = line.split(' ')[1]
            value = line.split(' ')[2]

            sg.addValue(float(value))
            sequence1 += res1
            sequence2 += res2

            line = handle.readline()

        handle.close()

        sg.addSequence(seqid1, seqdesc1, sequrl1, sequence1)
        sg.addSequence(seqid2, seqdesc2, sequrl2, sequence2)

        seqgraphstring = cPickle.dumps(sg)
        print buildGraph(seqgraphstring)



if __name__ == '__main__':
    main()

