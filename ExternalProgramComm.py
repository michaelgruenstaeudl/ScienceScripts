#!/usr/bin/env python2
''' External Program Communication'''
__author__ = "Michael Gruenstaeudl"
__maintainer__ = "Michael Gruenstaeudl"
__email__ = "michael.gruenstaeudl@utexas.edu"
__version__ = "2012.08.30.1500"
__status__ = "Working"

#########################
### IMPORT OPERATIONS ###
#########################

import subprocess, os, commands, time

###############
### CLASSES ###
###############

class interactwithPaup():
    '''class for interacting with the shell-version of PAUP;
    needs inputfilename <a> and commandlist <b> as input.'''

    def __init__(self,a,b):
        self.inputfilename = a
        self.commandlist = b
        
    def executeCommands(self):
        
        process = subprocess.Popen("paup.exe", stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        process.stdin.write("log file="+str(self.inputfilename)+"_log.txt;"+"\n")
        process.stdin.write("execute "+str(self.inputfilename)+".nex;"+"\n")
        for i in self.commandlist:
            process.stdin.write(str(i)+"\n")        
        process.stdin.write("log stop;"+"\n")
        process.stdin.write("quit;"+"\n")
        output,error = process.communicate()
        #DEBUGLINE  print output, error


class interactwithMesquite():
    '''class for interacting with the headless version of Mesquite;
    needs command to start Mesquite <a> and commandlist <b> as input.'''

    def __init__(self,a,b):
        self.command_to_start_mesquite = a
        self.commandlist = b
    
    def executeCommands(self):
        #DEBUGLINE print self.command_to_start_mesquite
        #DEBUGLINE print self.commandlist
    	process = subprocess.Popen(self.command_to_start_mesquite, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
    	for i in self.commandlist:
    		process.stdin.write(i+"\n")
    		time.sleep(15)  # 20 is usually safer
    	output,error = process.communicate()
    	return [output,error]


class getbashanswer():
    '''class for obtaining bash answers;
    needs appropriate command <a> as input'''

    def __init__(self,a):
        self.command = a
        
    def go(self):
        process = subprocess.Popen(self.command, shell=True, stdout = subprocess.PIPE)
        output = process.communicate()[0]
        #DEBUGLINE print output
        return output.split("\n")[0]
        
        
###################
### DEFINITIONS ###
###################        

def InteractWithPAUP(a,b):
    return interactwithPaup(a,b).executeCommands()

def InteractWithMesquite(a,b):
    return interactwithMesquite(a,b).executeCommands()

def getBashAnswer():
    return getbashanswer("dir").go()

