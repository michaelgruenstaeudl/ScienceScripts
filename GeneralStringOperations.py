#!/usr/bin/env python2
'''General String Operations'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.06.07.2300"
__status__ = "Working"

#########################
### IMPORT OPERATIONS ###
#########################
# none

###############
### CLASSES ###
###############

class excisingStrings:
    '''class for excising strings via keyword delimitation;
       needs instring <a>, keyword1 <b>, and keyword2 <c> as input'''
       
    def __init__(self,a,b,c):
        self.instring = a           # string
        self.keyword1 = b           # string
        self.keyword2 = c           # string
        
    def go(self):
        pos1 = afterFind(self.instring,self.keyword1).go()
        pos2 = self.instring.find(self.keyword2,pos1)
        return self.instring[pos1:pos2]
        
class replacingStrings:
    '''class for replacing strings via keyword delimitation;
       needs instring <a>, keyword1 <b>, keyword2 <c>, and replacement <d> as input'''
       
    def __init__(self,a,b,c,d):
        self.instring = a           # string
        self.keyword1 = b           # string
        self.keyword2 = c           # string
        self.replacement = d        # string
        
    def go(self):
        pos1 = afterFind(self.instring,self.keyword1).go()
        pos2 = self.instring.find(self.keyword2,pos1)
        if self.instring[pos1:pos2]:
            output = self.instring.replace(self.instring[pos1:pos2],self.replacement)
        # If string to be replaced were empty:
        if not self.instring[pos1:pos2]:
            output = self.instring[:pos1]+self.replacement+self.instring[pos2:]
        return output

class generatingRelativeValues:
    '''class for generating relative values from a list of absolute values;
       needs inlist <a> as input'''
       
    def __init__(self,a):
        self.inlist = a
        
    def go(self):
        alist = self.inlist[2][self.inlist[2].find("each:")+5:].split(";")  
        summe = sum(float(x.strip()[3:]) for x in alist)
        return "\n".join(str(x)+(str(float(x.strip()[:3])/summe)[:4]) for x in alist)

class clearSplit:
    '''class for splitting string without removing delimiter;
       needs instring <a> and delimiter <b> as input'''
    def __init__(self,a,b):
        self.instring = a
        self.delimiter = b
        
    def go(self):
        return [self.delimiter+element for element in self.instring.split(self.delimiter)]

class afterFind:
    '''class for returning index number of end of keyword in instring;
       needs instring <a> and keyword <b> as input '''
    def __init__(self,a,b):
        self.instring = a
        self.keyword = b
                
    def go(self):
        pos1 = self.instring.find(self.keyword)+len(self.keyword)
        return pos1


class removeExtension:
    '''class for returning instring without extension (i.e. ".txt" or ".trees";
       needs instring <a> and delimiter <b> as input '''
    def __init__(self,a,b):
        self.instring = a
        self.delimiter = b
                
    def go(self):
        return self.instring[:self.instring.rfind(self.delimiter)]
        

###################
### DEFINITIONS ###
###################

def GenerateRelValues(inlist):
    return generatingRelativeValues(inlist).go()
    
def afind(a,b):
    return afterFind(a,b).go()

def csplit(a,b):
    return clearSplit(a,b).go()

def exstr(instring,keyword1,keyword2):
    return excisingStrings(instring,keyword1,keyword2).go()

def replstr(instring,keyword1,keyword2,replacement):
    return replacingStrings(instring,keyword1,keyword2,replacement).go()

def rmext(a):
    return removeExtension(a,".").go()
