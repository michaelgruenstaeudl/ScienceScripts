#!/usr/bin/env python2
''' External Program Communication'''
__author__ = "Michael Gruenstaeudl, PhD"
__maintainer__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.12.03.1900"
__status__ = "Working"

# IMPORT OPERATIONS

import commands
import os
import subprocess
import time


# CLASSES

class InteractWithPAUP():
    '''class for interacting with the shell-version of PAUP;
    needs inFn <a> and cmdlist <b> as input.'''

    def __init__(self,a,b):
        self.inFn = a
        self.cmdlist = b
        
    def executeCommands(self):        
        p = subprocess.Popen("paup.exe",
                             stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE)
        p.stdin.write("log file="+str(self.inFn)+"_log.txt;"+"\n")
        p.stdin.write("execute "+str(self.inFn)+".nex;"+"\n")
        for i in self.cmdlist:
            process.stdin.write(str(i)+"\n")        
        p.stdin.write("log stop;"+"\n")
        p.stdin.write("quit;"+"\n")
        output, error = p.communicate()


class getbashanswer():
    '''class for obtaining bash answers;
    needs appropriate command <a> as input'''

    def __init__(self,a):
        self.cmd = a
        
    def go(self):
        p = subprocess.Popen(self.cmd,
                             stdout = subprocess.PIPE,
                             shell=True)
        output = p.communicate()[0]
        return output.split("\n")[0]
        
        
# DEFINITIONS

def CmdsPAUP(a,b):
# formerly: InteractWithPAUP
    return InteractWithPAUP(a,b).executeCommands()

def getBashAnswer():
    return getbashanswer("dir").go()
