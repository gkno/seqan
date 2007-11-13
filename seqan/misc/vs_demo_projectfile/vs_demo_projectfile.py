import sys
import string
import os

BUILD_FLAGS_7 = '/I ".." /I "../library" /I "../../platforms/windows" /Op- /EHsc /D "DEBUG" /D "WIN32" /Zi /GR /W2 /Zc:wchar_t'

################################################################################

def addFile(name):
    global BUILD_FLAGS
    global vcproj1
    global vcproj2
    
    print ".",
    
    vcproj1 += '\t\t<Configuration Name="' + name + '|Win32" IntermediateDirectory="temp" ConfigurationType="0">\n'
    vcproj1 += '\t\t\t<Tool Name="VCNMakeTool"\n\t\t\t\tBuildCommandLine="&quot;$(VCInstallDir)bin\\cl.exe&quot; &quot;$(ConfigurationName).cpp&quot; ' + BUILD_FLAGS + '" Output="&quot;$(ConfigurationName).exe&quot;"/>\n'
    vcproj1 += '\t\t</Configuration>\n'

    vcproj2 += '\t\t<File RelativePath=".\\' + name + '.cpp"></File>\n'
    

################################################################################

def scanDemos(search_path):
    global vcproj1
    global vcproj2

    vcproj1 = ""
    vcproj2 = ""
    
    for root, dirs, files in os.walk(search_path):
        if 'CVS' in dirs: dirs.remove('CVS')
        if '.svn' in dirs: dirs.remove('.svn')
        for file in files:
            pos = file.rfind(".")
            if pos < 0: continue
            if file[pos+1:] in ["cpp", "CPP"]:
                addFile(file[:pos])

################################################################################

def createProject(search_path, version):
    global vcproj1
    global vcproj2

    f = open(version + ".vcproj")
    t = f.read()
    f.close()
    
    t = t.replace("####PLATZHALTER1####", vcproj1)
    t = t.replace("####PLATZHALTER2####", vcproj2)
    f = open(os.path.join(search_path, version + ".vcproj"), "w")
    f.write(t)
    f.close()
    

    f = open(version + ".sln")
    t = f.read()
    f.close()
    f = open(os.path.join(search_path, version + ".sln"), "w")
    f.write(t)
    f.close()

   
################################################################################
    

def main():
    global BUILD_FLAGS
    
    print "This script creates Visual Studio files (.vcproj and .sln) for SeqAn demos"
    print "The created files are stored in the demos folder"
    
    BUILD_FLAGS = BUILD_FLAGS.replace('"', '&quot;');
    
    scanDemos("..\\..\\projects\\demos")
    
    createProject("..\\..\\projects\\demos", "Seqan_7)
    createProject("..\\..\\projects\\demos", "Seqan_8)
    
    print "Files created."
    
    
main()