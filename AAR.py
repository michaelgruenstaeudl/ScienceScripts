#!/usr/bin/env python2
'''Conducting Ancestral Area Reconstructions (AARs) and Visualizing the reconstructions'''
__author__ = "Michael Gruenstaeudl, PhD"
__copyright__ = "Copyright (C) 2014 Michael Gruenstaeudl"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.06.08.1700"
__status__ = "Working"

#########################
### IMPORT OPERATIONS ###
#########################

from termcolor import colored
import argparse, dendropy, itertools, re, subprocess, os, sys
import ExternalProgramComm as EPC
import GeneralFileOperations as GFO
import GeneralStringOperations as GSO

####################
### DEBUG HELPER ###
####################
#print colored("DEBUG "+"results\n"+results, 'magenta')

########################
### GLOBAL VARIABLES ###
########################

executeMesquite = "java -Djava.library.path=lib -Xmx3000M -Xss4m -d64 \
                   -cp /home/michael/binaries/mesquite2.75/Mesquite_Folder/ \
                   mesquite.Mesquite"
executeTG2 = "java -jar /home/michael/binaries/treegraph2/TreeGraph.jar"

customChanges = {'<TextLabel Text="a':'<TextLabel Text="Sp1',
                 '<TextLabel Text="b':'<TextLabel Text="Sp2'}

colordict = {}
#colordict = {'a':'#0000FF',  # #0000FF: dark blue
#             'b':'#FF0000'}  # #FF0000: red

#colordict = {'A':'#0000FF', # Azores (0): dark blue
#             'C':'#FF00FF', # Cape Verde (1): pink
#             'E':'#00FF00', # Mainland (2): green
#             'G':'#FF8300', # La Gomera (3): orange
#             'H':'#00FFFF', # El Hierro (4): light blue
#             'M':'#FF0000', # Madeira (5): red
#             'P':'#FFA07A', # La Palma (6): salmon
#             'R':'#808000', # Gran Canaria (7): olive
#             'T':'#FFFF00'} # Tenerife (8): yellow

###############
### CLASSES ###
###############

class customizingXTG:
    '''class for parsing XTG files: turning raw XTG code into publication-ready XTG code'''
    
    def __init__(self,a,b):
        self.infile = a
        self.flags = b
    
    def go(self):
        outlist = []                                                    
        # Loop through list elements (except first, which is empty)
        for i in GSO.csplit(self.infile,"<")[1:]:

            # If line starts with keyword, do sth.
            if i[0:5] == "<Node":
                # Line width adjustment
                i = GSO.replstr(i,'LineWidth="','"',"0.5")
                # Edge radius adjustment
                i = GSO.replstr(i,'EdgeRadius="','"',"0.9")
                # Text size adjustments
                i = GSO.replstr(i,'TextHeight="','"',"5.0")
                # Text style adjustments
                i = GSO.replstr(i,'TextStyle="','"',"Italics")

            if i[0:5] == "<Bran":
                # Line width adjustment
                i = GSO.replstr(i,'LineWidth="','"',"0.5")

            if i[0:5] == "<Glob":
                # Root presence adjustment 
                if "NOROOT" in self.flags.upper():
                    i = i.replace('ShowRooted="true"','ShowRooted="false"')
                # Brlen adjustment 
                if "BRLX2" in self.flags.upper():
                    brlscale = GSO.exstr(i,'BranchLengthScale="','"')
                    adj_brlscale = str(float(brlscale)*2)
                    i = GSO.replstr(i,'BranchLengthScale="','"',adj_brlscale)
                if "BRLX3" in self.flags.upper():
                    brlscale = GSO.exstr(i,'BranchLengthScale="','"')
                    adj_brlscale = str(float(brlscale)*3)
                    i = GSO.replstr(i,'BranchLengthScale="','"',adj_brlscale)

            outlist.append(i)

        return ''.join(outlist)

###############################################################################

class addingPieCharts:
    '''class for adding pie charts to phylogenetic trees in XTG format'''
    
    def __init__(self,a,b,c,d,e):
        self.infn = a
        self.modelfile = b
        self.flags = c
        self.instring = d
        self.piedata = e
    
    def go(self):
        string_1a = '<PieChartLabel LineColor="#000000" LineWidth="0.3" Width="10.0"\
                     Height="10.0" InternalLines="true" NullLines="false" Id="internals"\
                     Above="true" LineNo="0" LinePos="0"><LabelMargin Left="1.0" Top="1.0"\
                     Right="1.0" Bottom="1.0"></LabelMargin><DataIds>'
        # string_1b is optimal for phylograms
        string_1b = '<PieChartLabel LineColor="#000000" LineWidth="0.2" Width="3.5"\
                     Height="3.5" InternalLines="true" NullLines="false" Id="internals"\
                     Above="true" LineNo="0" LinePos="0"><LabelMargin Left="1.0" \
                     Top="0.0" Right="1.0" Bottom="1.0"></LabelMargin><DataIds>'
        string_2 = '</DataIds></PieChartLabel>'

        infile = self.instring
        outstring, errorlist = "", ["\nWARNING - Adding Pie Charts:"]
      
#   1. Inserting pie chart XTG lines into infile by looping through self.piedata values
        # How to parse node numbers:
        # Snippet: i[0][-3:-1]                                          -> gives the node number
        # Snippet: i[0][i[0].find(" ")+1:i[0].find(";")]                -> extracts the node number, irrespective of its numbers of digits
        # Snippet: i[0][i[0].find(" ")+1:i[0].find(":")].strip()+'"'    -> extract node number, even better

        # Looping through piedata entries
        # --- Start of loop ---
        for i in self.piedata:
            # extracting the node number, irrespective of its numbers of digits
            #node = 'Node UniqueName="Node_'+i[0][-3:-1].strip()+'"' 
            #node = 'Node UniqueName="Node_'+i[0][i[0].find(" ")+1:i[0].find(";")].strip()+'"'
            node = i[0][i[0].find(" ")+1:i[0].find(":")].strip()
            nodeinfo = 'UniqueName="Node_'+node+'"'

            if any(["NaN" in e for e in i]):
                print colored("    Warning: Node "+node+" was not present in reconstruction trees.", 'magenta')
                errorlist.append("  Node "+node+" was not present in reconstruction trees.")

            if nodeinfo not in infile:
                print colored("    Warning: Node "+node+" was not present in plotting tree.", 'magenta')
                errorlist.append("  Node "+node+" was not present in plotting tree.")

            # Remember that '<Node' and 'UniqueName="Node_1"' may not always be adjacent to one another
            if not any(["NaN" in e for e in i]) and "<Node" in infile and nodeinfo in infile:
                # pos = index position immediately before keyword '</Branch>', but after node  
                pos = infile.find('</Branch>',infile.find(nodeinfo))

    # What a typical "i" looks like: ['node 2:', '50', 'a: 25.0', 'b: 25.0'] under MultiTree ML-reconstr.
    # Option 1
                if len(i) > 2 and "." not in i[2]:
                    if "MP" in self.modelfile.upper():
                        print "    - adding pie charts for node: "+node
                    else:
                        print colored("    Warning: Flags and settings don't match for node: "+node,'magenta')

                    #Example element of relvalue_list
                    relvalue_list = GSO.GenerateRelValues(i)
                    outstring = infile[:pos]+string_1a          # Append first infile-section and string_1a
                                                                # to newly-created outstring                        
                    for j in relvalue_list:                     # Loop through relvalue_list and append
                                                                # PieColor definitions to outstring
                        outstring += '<DataId PieColor="'+colordict[j[0]]+'">'+j[0].strip('\n')+'</DataId>'
                    outstring += string_2                       # Append string_2 to outstring
                    
                    for j in relvalue_list:                     # Loop through relvalue_list and append
                                                                # InvisibleData definitions to outstring
		                outstring += '<InvisibleData Id="'+j[0]+'" Text="'+j[2:]+'" IsDecimal="true"></InvisibleData>'

    # What a typical "i" looks like: ['node 2:', '50', 'a: 25.0', 'b: 25.0'] under MultiTree ML-reconstr.
    # Option 2
                if len(i) >= 2 and "." in i[2]:
                    if "ML" in self.modelfile.upper():
                        print "    - adding pie charts for node: "+node
                    else:
                        print colored("    Warning: Flags and settings don't match for node: "+node,'magenta')

                    alist=[]
                    # i[2:], because we no longer need the other info (i[0] is node name, i[1] is number of trees)
                    for j in i[2:]:
                        if "E-" not in j and "0.00" not in j:
                            alist.append(j[:8])
                        else:
                            print colored("    Warning: Strange values for node: "+node,'magenta')                           
                    # Append first infile-section and string_1b to newly-created outstring
                    outstring = infile[:pos]+string_1b                        
                    # Loop through alist and append PieColor definitions to outstring
                    for k in alist:
                        outstring += '<DataId PieColor="'+colordict[k.split(": ")[0]]+'">'+k.split(": ")[0].strip('\n')+'</DataId>'
                    # Append string_2 to outstring
                    outstring += string_2
                    
                    # Loop through relvalue_list and append InvisibleData definitions to outstring
                    for k in alist:
		                outstring += '<InvisibleData Id="'+k.split(": ")[0]+'" Text="'+k.split(": ")[1].strip("*")+'" IsDecimal="true"></InvisibleData>'

    # Option 3            
                if len(i) <= 2:
                    if "MP" in self.modelfile.upper() or "ML" in self.modelfile.upper():
                        print "   - adding pie charts for node: "+node
                    else:
                        print colored("    Warning: Flags and settings don't match for node: "+node,'magenta')
                  
                                                                # Piedata entires have at most two elements,
                    absvalue_list = i[2].split(",")             # Generate absvalue_list by splitting
                                                                # area abbreviations into list
                                                                # Example absvalue_list: ['T','E']

                    outstring = infile[:pos]+string_1b          # Append first infile-section and string_1b
                                                                # to newly-created outstring                        
                    
                    for j in absvalue_list:                     # Loop through absvalue_list and append
                                                                # DataId definitons to outstring
                        outstring += '<DataId PieColor="'+colordict[j]+'">'+j+'</DataId>'
                    outstring += string_2                       # Append string_2 to outstring
                    
                    for j in absvalue_list:                     # Loop through absvalue_list, turn into
                                                                # rel_value and append InvisibleData
                                                                # definitions to outstring
                        rel_value = str(1.0/len(i[2].split(",")))
                        outstring += '<InvisibleData Id="'+j+'" Text="'+rel_value[:3]+'" IsDecimal="true"></InvisibleData>'

                # Close outstring with second infile-section
                outstring += infile[pos:]
                # IMPORTANT STEP: outstring becomes the new infile of the loop
                infile = outstring
        # --- End of loop ---

        # Important line, don't delete
        outlist = GSO.csplit(outstring,"<")[1:]

#   2. Save errorlist to file
        if errorlist >= 2:
            GFO.append(GSO.rmext(self.infn)+"_OmittedNodes.txt",'\n'.join(errorlist))

#   3. Switch position of bootstrap values to below branches, bc they would confilict with pie charts
        if "BS" in self.flags.upper():
            newoutlist = []
            keyword = '<TextLabel Text="'

            for line in outlist:
                if keyword in line and GSO.afind(line,keyword) in [str(s) for s in [1]+range(5,10)]:
                    bsvalue = GSO.exstr(line,'Text="','"')[:-2]
                    line = GSO.replstr(line,'Text="','"',bsvalue)
                    line = line.replace('IsDecimal="true"','IsDecimal="false"')
                    line = line.replace('Above="true"','Above="false"')
                    newoutlist.append(line)
                if keyword in line and GSO.afind(line,keyword) in [str(s) for s in [0]]:
                    line = line.replace('Above="true"','Above="false"')
                    newoutlist.append(line)
                else:
                    newoutlist.append(line)

            # Important line, don't delete
            outlist = newoutlist

#   4. Return outlist                
        return '\n'.join(outlist)

###############################################################################

class addingPieLabels:
    '''class for adding pie labels to phylogenetic trees in XTG format'''
       
    def __init__(self,a,b,c):
        self.infn = a
        self.instring = b
        self.piedata = c
        
    def go(self):
        string_3 = '<LabelMargin Left="1.0" Top="0.0" Right="1.0" Bottom="0.0"></LabelMargin></TextLabel>'
        infile = self.instring
        outstring, errorlist = "", ["\nWARNING - Adding Pie Labels:"]

#   1. Inserting pie chart XTG lines into infile by looping through inlist values
        # --- Start of loop ---
        for i in self.piedata:
            # extracting the node number, irrespective of its numbers of digits
            #node = 'Node UniqueName="Node_'+i[0][-3:-1].strip()+'"' 
            #node = 'Node UniqueName="Node_'+i[0][i[0].find(" ")+1:i[0].find(";")].strip()+'"'
            node = i[0][i[0].find(" ")+1:i[0].find(":")].strip()
            nodeinfo = 'UniqueName="Node_'+node+'"'

            # following line checks in any of the list elements contain "NaN" as part of their strings (e.g. ['node 5:', '0', 'a: NaN', 'b: NaN'])
            if any(["NaN" in e for e in i]):
                print colored("    Warning: Node "+node+" was not present in reconstruction trees.", 'magenta')
                errorlist.append("  Node "+node+" was not present in reconstruction trees.")

            if nodeinfo not in infile:
                print colored("    Warning: Node "+node+" was not present in plotting tree.", 'magenta')
                errorlist.append("  Node "+node+" was not present in plotting tree.")

            if not any(["NaN" in e for e in i]) and nodeinfo in infile:
                # pos = index position immediately before keyword '</Branch>', but after node  
                pos = infile.find('</Branch>',infile.find(nodeinfo))
                outstring = infile[:pos]

# Option 1
                # If an MP-driven reconstr over a multiple trees
                if len(i) > 2 and "." not in i[2]:
                    print "    - adding pie labels for node: "+node
                    alist = GSO.GenerateRelValues(i)

# Option 2
                # If sets of ML-driven reconstr.               
                if len(i) >= 2 and "." in i[2]:                      
                    print "    - adding pie labels for node: "+node
                    alist=[]
                    for j in i[2:]:
                        if "E-" not in j and "0.00" not in j:
                            if "*" in j and "1.0" not in j:
                                alist.append(j.split(": ")[0]+" "+str(( float(j.split(": ")[1].strip("*")) / float(i[1]) ))[:4] + "*")
                            else:
                                alist.append(j.split(": ")[0]+" "+str(( float(j.split(": ")[1].strip("*")) / float(i[1]) ))[:4] )

# Option 3a
                # If an ML-driven reconstr over a single tree
#                if len(i) <= 2:
#                    print "  >   SingleTree MP-reconstr. | Adding pie labels for node: "+node                    
#                    alist = i[1].split(",")

#   2. Actual addition of pie labels         
                for counter,i in enumerate(alist,start=1):
                    pielabel_string = '<TextLabel Text="'+i+'" IsDecimal="false" TextColor="#000000" TextHeight="2.5" TextStyle="" FontFamily="Arial" DecimalFormat="#0.0#####" LocaleLang="en" LocaleCountry="" LocaleVariant="" Id="internals" Above="true" LineNo="'+str(counter)+'" LinePos="0">'
                    outstring += pielabel_string
                    outstring += string_3

                # Close outstring with second infile-section
                outstring += infile[pos:]

                # IMPORTANT STEP: outstring becomes the new infile of the loop
                infile = outstring                                  
        # --- End of loop ---

#   3. Save errorlist to file
        if errorlist >= 2:
            GFO.append(GSO.rmext(self.infn)+"_OmittedNodes.txt",'\n'.join(errorlist))

#   4. Return outstring
        return outstring

###############################################################################

class conversion_NEX2XTG:
    '''class for performing nex2xtg conversion in TreeGraph2 via commandline;
       needs inputfile name <a> and command to start TreeGraph2 <b> as input'''
       
    def __init__(self,a):
        self.infn = a
        
    def go(self):
        # Get current working directory
        cwd = os.getcwd()

#     1 Generating in- and output specs for nex2xtg
        # Note: The os.path.join function constructs pathnames by concatenating 
        # strings, while accounting for Windows pathnames, which require 
        # backslash characters to be escaped.
        infile = os.path.join(cwd, self.infn)
        outfilename = os.path.join(cwd, GSO.rmext(self.infn)+".xtg")
        
#     2. Performing nex2xtg via Treegraph2 and reporting output/error
        command = executeTG2 + " -convert " + infile + " -xtg " + outfilename
        process = subprocess.Popen(command, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        output,error = process.communicate()
        TreeGraph2_ErrorReport(output,error)

#     3. Loading xtg-file into memory
        return GFO.loadR(outfilename)

###############################################################################

class conversion_XTG2IMG:
    '''class for performing XTG2IMG conversion in TreeGraph2 via commandline'''
      
    def __init__(self,a,b):
        self.infn = a
        self.flags = b
        
    def go(self):
        cwd = os.getcwd()  
        
#     1.a. Checking whether xtg-file is present; if yes, generating in- and output specs
        if os.path.exists(os.path.join(cwd, GSO.rmext(self.infn)+".xtg")) == False:
            sys.exit("Problem: Conversion .nex -> .xtg has NOT worked!\nStopping script...")
        if os.path.exists(os.path.join(cwd, GSO.rmext(self.infn)+".xtg")) == True:
            infile = os.path.join(cwd, GSO.rmext(self.infn)+".xtg")
            outfilename1 = os.path.join(cwd, GSO.rmext(self.infn)+".png")
            outfilename2 = os.path.join(cwd, GSO.rmext(self.infn)+".svg")
            
#     1.b. Checking type of tree (phylo- or cladogram) and proceeding accordingly
            if "CLADO" in self.flags.upper():
                resol_command = " -width 600mm -res 120ppi"
            if "PHYLO" in self.flags.upper():
                resol_command = " -phyl -width 600mm -res 120ppi"
                           
#     1.c. Performing xtg2png via Treegraph2 and reporting output/error
            command = executeTG2 + " -image " + infile + " " + outfilename1 + resol_command
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output,error = process.communicate()
            TreeGraph2_ErrorReport(output,error)

#     1.c. Performing xtg2svg via Treegraph2 and reporting output/error
            command = executeTG2 + " -image " + infile + " " + outfilename2 + resol_command
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output,error = process.communicate()
            TreeGraph2_ErrorReport(output,error)

###############################################################################

class TreeGraph2_ErrorReport:
    '''class for reporting errors during TreeGraph2 execution'''
    
    def __init__(self,a,b):
        self.output_stream = a
        self.error_stream = b

    def go():
        if not self.error_stream and not self.output_stream:
            if "java.lang.NullPointerException" in self.output_stream:
                print self.output_stream
                print "Possible problem: Path to cwd must NOT have spaces!"
            else:
                print colored("  > Error: ", 'magenta')+"Unknown TreeGraph2 error."
        if self.error_stream:
            print colored("  > Error: ", 'magenta')+self.error_stream

###############################################################################

class labelingNodesXTG:
    '''class for labeling nodes of phylogenetic tree in XTG format'''
    
    def __init__(self,a):
        self.instring = a
        
    def go(self):
        outlist = []
        alist = GSO.csplit(self.instring,"<")[1:]
        # Loop through list elements (except first, which is empty)
        counter = 2
        for element in alist:
            # If element contains keyword 'UniqueName=', replace characters 
            # within double quotation marks with counter value
            if 'UniqueName="' in element:
                element = GSO.replstr(element,'UniqueName="','"',"Node_"+str(counter))
                counter += 1
            outlist.append(element)

        return ''.join(outlist)

###############################################################################

class humanreadableXTG:
    '''class for adding tabs to make XTG code human-readable'''
    
    def __init__(self,a):
        self.instring = a
        
    def go(self):
        adict = {'<Branch ':'\t<Branch ',
                 '</Branch>':'\t</Branch>',
                 '<LabelMargin ':'\t\t<LabelMargin ',
                 '<TextLabel ':'\t\t<TextLabel ',
                 '<PieChartLabel ':'\n\t\t<PieChartLabel ',
                 '</PieChartLabel>':'\t\t</PieChartLabel>',
                 '<DataId ':'\t\t\t<DataId ',
                 '<DataIds>':'\t\t\t<DataIds>',
                 '</DataIds>':'\t\t\t</DataIds>',
                 '<InvisibleData ':'\t\t\t<InvisibleData ',
                 '\n</DataId>':'</DataId>',
                 '\n</InvisibleData>':'</InvisibleData>',
                 '\n</LabelMargin>':'</LabelMargin>',
                 '\n</LeafMargin>':'</LeafMargin>',
                 '\n</TextLabel>':'</TextLabel>'}
        outstring = self.instring
        for key,value in adict.iteritems():
            outstring = outstring.replace(key,value)

        return outstring


###############
### MODULES ###
###############

def formhelp(results,filename,suffix):  
    GFO.save(GSO.rmext(filename)+suffix, results)
    return GSO.rmext(filename)+suffix

########################################################################

def main(infname, modfname, transl, flags):

# 0. Check if colordict set
    print "  > Checking color dictionary"
    if not colordict:
        print colored("  > Error:", 'magenta'),"No color dictionary specified."
        sys.exit("\n")
    if colordict:
        print "    - current color dictionary: "+colored(str(colordict), 'yellow')
        
# 1. Loading data sets
    indata = GFO.loadR(infname)
    modeldata = GFO.loadR(modfname)
        
# 2. Adjusting taxon names and saving temporary files
#   2.1. Adjusting taxon names (if translation file present)
    if transl:
        print "  > Adjustment of sequence names"
        table = GFO.loadRL(transl)
        # table has title line, hence [1:]
        for line in table[1:]:
            values = line.split(",")
            indata = indata.replace(values[0],values[1])

#   2.2. Checking for number of tree blocks
    if indata.lower().count("begin trees;") < 2:
        print colored("  > Error: ", 'magenta')+"Only a single tree block defined, but minimum is two."
    if indata.lower().count("begin trees;") >= 2:
        pos1 = indata.lower().rfind("begin trees;")
        pos2 = indata.lower().find("end;",pos1)+len("end;")

        # Converting tree from nexus into newick format, because nexus format may contain
        # translation table, which TreeGraph2 cannot parse
        treenex = "#NEXUS\n"+indata[pos1:pos2]
        treehandle = dendropy.Tree.get_from_string(treenex, schema="nexus")
        plottree = treehandle.as_newick_string()

#   2.3. Setting filenames for temporary files and saving to file
    fileprefix = GSO.rmext(infname)+"."+GSO.rmext(modfname)
    tmpinfname = fileprefix+".INPUT"
    treefname =  GSO.rmext(infname)+".plotting.nwk"
    tmpoutname = fileprefix+".OUTPUT"

    # Note: Since Mesquite cannot load filenames containing underscores, a temporary infile 
    # without underscores is generated, but deleted immediately after execution.
    GFO.save(tmpinfname,'\n'.join([indata,modeldata]))
    GFO.save(treefname,plottree)
    print "  > Plotting tree saved to file as: "+treefname

# 3. Character Reconstruction
    print "  > Character Reconstruction in Mesquite"

#   3.1. Execute infile and delete tempinfile
    handle = EPC.InteractWithMesquite(executeMesquite,["cd "+os.getcwd()+"/", "openFile "+tmpinfname, "quit"])
    # ETSUCl needs special command, bc it requests data to be saved and hangs as script; solution: any unknown command breaks that hangup
    #handle = EPC.InteractWithMesquite(executeMesquite,["cd "+os.getcwd()+"/","openFile "+filenamedict["INPUT"].replace("_",""), "BreakHangup", "quit"])
    if not handle:
        print colored("  > Error:", 'magenta'),"No character reconstruction from Mesquite received."
        sys.exit("\n")
    if colordict:
        print "    - character reconstruction from Mesquite received"
    GFO.deleteFile(tmpinfname)

#   3.2. Loading the relevant section of the reconstruction output and saving it to file
    output = GSO.exstr(handle[0],"Reading block: MESQUITE","File reading complete")
    # for ETSUCl only (see 3.3. for explanation)
    #output = GSO.exstr(handle[0],"Original Tree","Module Mesquite")
    GFO.save(tmpoutname,output)

#   3.3. More parsing of the reconstruction output for subsequent piedata generation
    parsing_keywords = ["Trace Character Over Trees","Trace Character History"]
    for keyword in parsing_keywords:
        if keyword in output:
            output = output[output.find(keyword):]
    if "File reading complete" in output:
        output = output[output.find(keyword):]   

# 4. Parsing the reconstruction output
    linesepList = output[output.find("\nnode"):].split("\n")
    # removing all empty elements of linesepList
    linesepList = filter(None, linesepList)
    # generating piedata
    piedata = [element.split("  ") for element in linesepList if element]
    if not piedata:
        print colored("  > Error:", 'magenta'),"General parsing of reconstruction data unsuccessful."
        sys.exit("\n")
    if colordict:
        print "    - reconstruction data parsed"

#   4.a. Special parsing for ML over multiple trees
    ### section was unintentionally deleted

#   4.b. Special parsing for ML over multiple trees
    if "ML" in modfname and "OverMultTrees" in modfname:
        print "  > Special data parsing under scheme 'MultiTree ML-reconstr.'"
        parsedpiedata = []                
        for i in piedata:
            e1 = i[0]
            e2 = GSO.exstr(i[1],'Node in ',' trees.')
            e3 = i[2].split("each: ")[1].split("; ")
            parsedpiedata.append([e1] + [e2] + e3)
        piedata = parsedpiedata
        if not piedata:
            print colored("  > Error:", 'magenta'),"Special data parsing unsuccessful. Possible issue: Malformed NEXUS file."
            sys.exit("\n")

# 5. Conversion TRE->PNG including all associated steps (nex->xtg, 
#    labelling nodes, adding pie labels and charts, etc.)
    print "  > Visualization in TreeGraph2"

#   5.1. Conversion from .nex to .xtg format
    print "  > Step 1: Conversion .nex -> .xtg"        
    results = conversion_NEX2XTG(treefname).go()

#   5.2. Label internal nodes
    print "  > Step 2: Labeling nodes of tree"
    results = labelingNodesXTG(results).go()

#   5.3. Format tree
    print "  > Step 3: Customizing XTG file"
    results = customizingXTG(results,flags).go()
    if piedata:
#   5.4. Add Pie labels
        # Indentation important, because pie data if statement above
        print "  > Step 4: Adding pie labels to tree"
        results = addingPieLabels(infname, results, piedata).go()
#   5.5. Add Pie data
        print "  > Step 5: Adding pie charts to tree"
        results = addingPieCharts(infname, modfname, flags, results, piedata).go()
    else:
        print colored("  >  Warning:",'magenta'),"Skipping Steps 4 and 5: No pie data available."

#   5.6. Improving visualization
    print "  > Step 6: Improving visualization"
    outlist = []
    # following section can't be in customizingXTG(), because lines starting with
    # "\t\t<Text" don't exist there yet; they are added in 
    for i in GSO.csplit(results,"<")[1:]:
        if i.strip("\t")[0:5] == "<Text":
            # Text color adjustment
            i = GSO.replstr(i,'TextColor="','"',"#808080")
            # Check if bootstrap values present; if so, strip comma and digits
            if "BS" in flags.upper():
                bsvalue = GSO.exstr(i,'Text="','"')[:-2]                     
                i = GSO.replstr(i,'Text="','"',bsvalue)
                i = i.replace('IsDecimal="true"','IsDecimal="false"')
        outlist.append(i)
    results = ''.join(outlist)    

    if customChanges:
        for key,value in customChanges.iteritems():
            results = results.replace(key,value)

#   5.7. Formatting .xtg file
    print "  > Step 7: Making XTG code human-readable"
    results = humanreadableXTG(results).go()
    tmpfn = GSO.rmext(infname)+".mod.xtg"
    GFO.save(tmpfn, results)

#   5.8. Conversion from .xtg to .png format
    print "  > Step 8: Conversion .xtg  -> .png"
    conversion_XTG2IMG(tmpfn, flags).go()

###############
### EXECUTE ###
###############

print ""
print colored("  Script name: "+sys.argv[0], 'cyan')
print colored("  Author: "+__author__, 'cyan')
print colored("  Version: "+__version__, 'cyan')
print colored("  (Note: Mesquite cannot process filenames with underscores.)", 'yellow')
print ""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Performing AAR; 2014 Michael Gruenstaeudl')
    parser.add_argument('-i','--infile', help='/path_to_working_dir/infile.nex', required=True)
    parser.add_argument('-m','--model', help='/path_to_working_dir/model.nex', required=True)
    parser.add_argument('-t','--translation', help='/path_to_working_dir/translation_file.csv', required=False)
    parser.add_argument('-f','--flags', help='BS(bootstrap),\
                        CLADO(cladogram),\
                        PHYLO(phylogram),\
                        NOROOT(dont display root),\
                        BRLX2(multiply branch length scale by 2),\
                        BRLX3(multiply branch length scale by 3)', required=True)
    args = parser.parse_args()

main(args.infile, args.model, args.translation, args.flags)

print ""
print colored("  Done.", 'cyan')
print ""
