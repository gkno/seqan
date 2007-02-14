import os
import copy
import string
import sys

################################################################################

FUNCS = {}
CLASSES = {}
TYPEDEFS = {}


################################################################################

def main(project_path):
    pos1 = project_path.rfind('/')
    if (pos1 < 0):
        exit ('ERROR: wrong argument "' + project_path + '"');
        
    project = project_path[pos1+1:]
    
    for root, dirs, files in os.walk(project_path):
        for file in files:
            if file == forwardFilename(project):
                continue
            path = os.path.join(root, file)
            if testFileType(path):
                parseFile(path)
        if 'CVS' in dirs:
            dirs.remove('CVS')
        if '.svn' in dirs:
            dirs.remove('.svn')
    
    if FUNCS != {}:        
        outAll(project_path, project)

################################################################################

def forwardFilename(project):
    return project + "_generated_forwards.h"

################################################################################

def testFileType(filename):
    pos = filename.rfind(".")
    if (pos >= 0): ext = filename[pos+1:]
    else: ext = ""

    return ext in ["c", "C", "cpp", "CPP", "c++", "C++", "h", "H", "hpp", "HPP", "h++", "H++"]
  
################################################################################
    
def parseFile(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    
    for line in lines:
        if (line.find("SEQAN_NO_GENERATED_FORWARDS") >= 0):
            print "-",
            return
    
    print ".",

    sigs = preprocess(lines, filename);
    
    createEntries(sigs)

################################################################################
# returns True if line ends with '\'

def isMultiLine(line):
    if (len(line) == 0): return False
    return line[len(line)-1] == '\\'
    
################################################################################
# removes comments, linebreaks, and macro definitions returns a string

def preprocess(lines, filename):
    ret = []
    
    inComment = False
    inString = False
    inDefine = False
    curlyCount = 0
    namespaces = []
    lineNumber = 0
    
    str = ""
    
    for line in lines:
        lineNumber += 1;

        line = line.strip()
    
        #remove some characters to make parsing simpler
        line = line.replace("'}'", "");
        line = line.replace("'{'", "");
        line = line.replace("';'", "");
        line = line.replace("'\"'", "");
        line = line.replace('\\"', "");
        

        #skip multiline defines
        if inDefine or ((not inComment) and (not inString) and (line.find("#define") >= 0)):
            inDefine = isMultiLine(line)
            continue
            
        if len(line) <= 0: 
            continue

        #skip preprocessor lines
        if (not inComment) and (not inString) and (line[0]=="#"):
            continue

        while len(line) > 0 :
            if inComment:
                pos1 = line.find("*/")
                if pos1 >= 0:
                    line = line[pos1 + 2:]
                    inComment = False
                else:
                    break
                    
            elif inString:
                pos1 = line.find('"')
                if pos1 >= 0:
                    line = line[pos1 + 1:]
                    inString = False
                else: 
                    break
                    
            else:
                pos1 = line.find("/*")
                pos2 = line.find("//")
                pos3 = line.find('"')
                pos4 = line.find('{')
                pos5 = line.find('}')
                pos6 = line.find(';')
                
                if (pos1 >= 0) and ((pos2 < 0) or (pos1 < pos2)) and ((pos3 < 0) or (pos1 < pos3)) and ((pos4 < 0) or (pos1 < pos4)) and ((pos5 < 0) or (pos1 < pos5)) and ((pos6 < 0) or (pos1 < pos6)):
                    if curlyCount == 0:
                        str += ' ' + line[0:pos1].strip()
                    line = line[pos1 + 2:]
                    inComment = True
                    
                elif (pos2 >= 0) and ((pos3 < 0) or (pos2 < pos3)) and ((pos4 < 0) or (pos2 < pos4)) and ((pos5 < 0) or (pos2 < pos5)) and ((pos6 < 0) or (pos2 < pos6)):
                    if curlyCount == 0:
                        str += ' ' + line[0:pos2].strip()
                    break
                    
                elif (pos3 >= 0) and ((pos4 < 0) or (pos3 < pos4)) and ((pos5 < 0) or (pos3 < pos5)) and ((pos6 < 0) or (pos3 < pos6)):
                    if curlyCount == 0:
                        str += ' ' + line[0:pos3].strip()
                    line = line[pos3 + 1:]
                    inString = True
                    
                elif (pos4 >= 0) and ((pos5 < 0) or (pos4 < pos5)) and ((pos6 < 0) or (pos4 < pos6)):
                    if curlyCount == 0:
                        entry = process2(str + ' ' + line[:pos4])
                        nam = isNamespace(entry)
                        if nam != "":
                            namespaces += [nam]
                        else:             
                            ret = ret + [[filename, lineNumber, entry, namespaces + []]]
                            curlyCount = 1
                        str = ""
                    else:
                        curlyCount += 1

                    line = line[pos4 + 1:]
                    
                elif (pos5 >= 0) and ((pos6 < 0) or (pos5 < pos6)):
                    line = line[pos5 + 1:]
                    if curlyCount > 0:
                        curlyCount -= 1
                    elif len(namespaces) > 0:
                        namespaces = namespaces[:len(namespaces)-1]
                    else:
                        print "ERROR in" , filename , "(", lineNumber, "): Too many }"
                        
                elif (pos6 >= 0):
                    if curlyCount == 0:
                        entry = process2(str + ' ' + line[:pos6])
                        if isStuctDeclaration(entry) or isTypedef(entry):
                            ret = ret + [[filename, lineNumber, entry, namespaces + []]]
                            
                    str = ""
                    line = line[pos6 + 1:]
                        
                else:
                    if curlyCount == 0:
                        str += ' ' + line
                    break
    
    return ret

################################################################################
# shrinks to the interesing part of the signature

def process2(str):
    str = str.replace('\t', ' ')
    str = str.replace('  ', ' ')
    
    pos1 = str.rfind(';')
    if pos1 >= 0:
        str = str[pos1+1:]
        
    str = removeBases(str)
        
    str = str.replace('template<', 'template <')
    str = str.replace('template < ', 'template <')
    str = str.replace(' operator ', ' operator')
        
    str = str.strip()
    return str;


################################################################################
# determines whether a signature is a namespace
# returns the namespace name or "" if it is no namespace

def isNamespace(str):
    pos1 = str.rfind(' ')
    if pos1 < 0: return ""
    while (pos1 > 0) and (str[pos1] == ' '): pos1 -= 1
    if (pos1 >= 8) and (str[pos1-8:pos1+1] == 'namespace'):
        return str[pos1 + 1:].strip()
    else:
        return ""


################################################################################
# determines whether a signature is a struct declaration

def isStuctDeclaration(str):
    str = removeBases(str)
    
    pos1 = str.rfind(' ')
    if pos1 < 0: return False
    while (pos1 > 0) and (str[pos1] == ' '): pos1 -= 1
    return ((pos1 >= 5) and (str[pos1-5:pos1+1] == 'struct')) or ((pos1 >= 4) and (str[pos1-4:pos1+1] == 'class'))
   
################################################################################
# determines whether a signature is typedef

def isTypedef(str):
    str = str.strip()
    return (str[:7] == "typedef") and (str.find('::') < 0)

################################################################################
# removes list of base class 

def removeBases(str):
    pos1 = -2
    while pos1 < len(str):
        pos1 = str.find(':', pos1 + 2)
        if pos1 < 0: return str
        if (pos1 == len(str)-1) or (str[pos1+1]!=':'): return str[:pos1]
        
    return str.strip()
    
################################################################################
# gets a list of [filename, lineNumber, signature, namespaces] tupels
# analyse signature and adds entries in FUNCS and CLASSES

def createEntries(sigs):
    for data in sigs:
        filename = data[0]
        lineNumber = data[1]
        sig = data[2]
        namespaces = data[3]
        
        entry = makeEntry(filename, lineNumber, sig)

        name = getTypedefName(sig)
        if name != '':
            addEntry(TYPEDEFS, name, entry, namespaces)
        else:
            name = getStructName(sig)
            if name != '':
                addEntry(CLASSES, name, deleteDefaultArguments(entry, '<', '>'), namespaces)
            else:
                name = getFuncName(sig)
                if name != '':
                    addEntry(FUNCS, name, deleteDefaultArguments(entry, '(', ');'), namespaces)


################################################################################
# deletes all default arguments from argument lists
# use delim = '>' for template argument lists and ')' for function argument lists

def deleteDefaultArguments(str, start_delim, stop_delim):
    ret = ""

    start = str.find(start_delim);
    if start >= 0: 
        ret = str[:start]
        str = str[start:]

    while str != "":
        pos1 = str.find("=")
        if pos1 < 0:
            ret += str
            break
        ret += str[:pos1]
        pos2 = str.find(",", pos1)
        pos3 = str.rfind(stop_delim, pos1)
        if ((pos2 > pos3) and (pos3 >= 0)) or (pos2 < 0): 
            pos2 = pos3
        if pos2 >= 0:
            str = str[pos2:]
        else:
            break
            
    return ret

################################################################################
# returns the string that is inserted into the header

def makeEntry(filename, lineNumber, sig):
    text = sig + ";       \t// \"" + filename + "\"(" + str(lineNumber) + ")"
    return text

################################################################################
# returns the key the functions and structs are sorted for

def getSortKey(name, namespaces):
    return str(namespaces) + name

################################################################################
# adds a signature to FUNCS or CLASSES

def addEntry(arr, name, entry, namespaces):
    key = getSortKey(name, namespaces)
    if not arr.has_key(key):
        arr[key] = []
    arr[key] += [[name, entry, namespaces]]
    

################################################################################
# get the function name from a signature or '', if signature is not a function

def getFuncName(sig):
    sig = sig.strip()
    pos1 = sig.rfind('(')
    if pos1 < 0:
        return ""
    else:
        pos1 -= 1
        while (pos1 >= 0) and (sig[pos1] == ' '):
            pos1 -= 1
            
        pos2 = sig.rfind(' ', 0, pos1)
        return sig[pos2 + 1: pos1 + 1]

################################################################################
# get the class name from a signature or '', if signature is not a class

def getStructName(sig):
    sig = sig.strip()

    pos1 = sig.rfind(' ')
    if pos1 < 0: return ""
    
    while (pos1 > 0) and (sig[pos1] == ' '): pos1 -= 1
    
    if ((pos1 >= 5) and (sig[pos1-5:pos1+1] == 'struct')) or ((pos1 >= 4) and (sig[pos1-4:pos1+1] == 'class')):
        name = sig[pos1 + 1:].strip()
        if (name.find('<') >= 0): return ""
        else: return name
       
    return ""


################################################################################
# get the typedef name from a signature or '', if signature is not a typedef

def getTypedefName(sig):
    sig = sig.strip()
    if not isTypedef(sig): return ""
    pos1 = sig.rfind(' ')
    if pos1 < 0: return ""
    return sig[pos1+1:].strip()


################################################################################
# main Function for output of forward header

def outAll(path, project):
    
    header_switch = "SEQAN_HEADER_" + project + "_GENERATED_FORWARDS_H"

    str = ""
    str += "#ifndef " + header_switch.upper() + " \n"
    str += "#define " + header_switch.upper() + " \n\n"
    
    str += "//////////////////////////////////////////////////////////////////////////////\n"
    str += "// NOTE: This file is automatically generated by build_forwards.py\n"
    str += "//       Do not edit this file manually!\n"
    str += "//////////////////////////////////////////////////////////////////////////////\n\n\n"

    str += "//////////////////////////////////////////////////////////////////////////////\n"
    str += "// CLASSES\n"
    str += outList(CLASSES)
    
    str += "\n//////////////////////////////////////////////////////////////////////////////\n"
    str += "// TYPEDEFS\n"
    str += outList(TYPEDEFS)

    str += "\n//////////////////////////////////////////////////////////////////////////////\n"
    str += "// FUNCTIONS\n"
    str += outList(FUNCS)
    
    str += "#endif\n"

    filename = os.path.join(path, forwardFilename(project))
    fl = file(filename, "w")
    fl.write(str + "\n")
    fl.close()

################################################################################
# 

def outList(lst):
    keys = lst.keys()
    keys.sort()
    
    namespaces = []
    str = ""
    
    for key in keys:
        first_entry = lst[key][0]
        new_namespaces = first_entry[2]
        str += outChangeNamespaces(namespaces, new_namespaces)
        namespaces = new_namespaces
            
        str += "//____________________________________________________________________________\n"
        str += "// " + first_entry[0] + "\n\n"
        
        for entry in lst[key]:
            str += entry[1] + "\n"
            
    str += outChangeNamespaces(namespaces, [])
    
    return str
    

################################################################################
# if old_namespaces != new_namespaces, close old namespace and open new one

def outChangeNamespaces(old_namespaces, new_namespaces):
    str = ""
    if old_namespaces != new_namespaces:
        if len(old_namespaces) > 0: 
            str += "\n"
        
        while len(old_namespaces) > 0:
            str += "} //namespace " + old_namespaces[len(old_namespaces)-1] + "\n"
            old_namespaces = old_namespaces[:len(old_namespaces)-1]
            
        if len(new_namespaces) > 0:
            str += "//////////////////////////////////////////////////////////////////////////////\n\n"
            
        while len(new_namespaces) > 0:
            str += "namespace " + new_namespaces[0] + " {\n"
            new_namespaces = new_namespaces[1:]
            
    return str + "\n"

################################################################################
# Start: parse arguments and call main

#main("../projects/library/seqan")


if len(sys.argv) < 2: 
    exit ('too few arguments');

#if (os.path.exists("V:/seqan2/" + sys.argv[1])):
print "create forwards for", sys.argv[1]
main (sys.argv[1]);

