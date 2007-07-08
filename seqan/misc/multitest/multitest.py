import subprocess
import sys
import string
import os
import getpass
from stat import *

################################################################################

PATH_MODULES = '../../projects/tests'
PATH_OUTPUT = ''
FILE_VIEWER = 'index.html'
FILE_OUTPUT = 'data.js'
FILE_OUTPUT_OLD = 'data_old.js'
PATH_RESULTFILES = 'results/'

LOAD_MAX = 2

# [plaform, server]
SERVERS = [
    ['gcc', 'mouse'], 
    ['gcc', 'frosch'], 
    ['gcc', 'kiefer'], 
    ['gcc', 'fichte'], 
    ['gcc', 'aplysia'], 
#    ['gcc', 'nawab'], 
#    ['gcc', 'bauch'], 
    ['windows', 'localhost']
]
    
# [plaform, name, compiler]
COMPILERS = [
    ['gcc',     'g++-3.0',   'g++-3.0'], 
    ['gcc',     'g++-3.2',   'g++-3.2'], 
    ['gcc',     'g++-3.3',   'g++-3.3'], 
    ['gcc',     'g++-3.4',   'g++-3.4'], 
    ['gcc',     'g++-4.1.1',  '/import/testing/bin/g++'],
#    ['windows', 'vc++-2003', 'cl']
]
    
MODES = ['release build']
#MODES = ['release build', 'release runonly']

################################################################################

################################################################################

MODULES = []
ALL_STATES = ''
USERNAME = 'unknown'
PASSWORD = ''

################################################################################

def main():
    init()
    createPage()
    showPage()
    while not finishedAll():
        scheduleTasks()
        if pollTasks(): createPage()
    pollTasks()
    createPage()


################################################################################

def init():
    global SERVERS, COMPILERS, MODES, PATH_MODULES, LOAD_MAX
    global PATH_OUTPUT, PATH_RESULTFILES
    global MODULES, TASKS, LOAD
    global USERNAME, PASSWORD
    
    USERNAME = getpass.getuser()
    PASSWORD = getpass.getpass()
    
    
    MODULES = []
    for f in os.listdir(PATH_MODULES):
        if (f == 'CVS') or (f == '.svn'): continue
        if f[0] == '_': continue
        p = PATH_MODULES + "/" + f
        m = os.stat(p)[ST_MODE]
        if S_ISDIR(m):
            MODULES.append(f)
            
    MODULES.sort()
    
    TASKS = []
    for mode in MODES:
        compnum = 0
        for [platform, compilername, compiler] in COMPILERS:
            compnum += 1
            for module in MODULES:
                task = {
                    'platform': platform,
                    'compilername': compilername,
                    'compiler': compiler,
                    'mode': mode,
                    'module': module,
                    'state': 'pending',
                    'filename': module + '_' + platform + '_' + compilername + '_' + mode,
                    'proc' : None                    
                }
                TASKS.append(task)

    servers = []
    i = 0
    for [platform, servername] in SERVERS:
        server = {
            'i': i,
            'platform': platform,
            'server': servername,
            'load': 0,
            'loadmax': LOAD_MAX,
            'count': 0
        }
        i += 1
        servers.append(server)
    SERVERS = servers
    
    path_resultfiles = PATH_OUTPUT + PATH_RESULTFILES
    if not os.access(path_resultfiles, os.F_OK):
        os.mkdir(path_resultfiles)

   
################################################################################
    
def javaScriptEscape(s):
    s = s.replace("'", "\\'")
    s = s.replace("\n", "'\n+ '")
    s = "DATA = '" + s + "';"
    return s;
    

def createPage():
    global SERVERS, COMPILERS, MODES, MODULES, TASKS
    global PATH_OUTPUT, FILE_OUTPUT, FILE_OUTPUT_OLD

    s = ''
    i = 0
    for mode in MODES:
        s += '<div class=mode>' + mode + '</div>\n'
        s += '<table cellspacing=0 cellpadding=0 border=0>\n'
        s += '  <tr><th>&nbsp;</th>'
        for module in MODULES:
            s += '<th>' + module + '</th>'
        s += '  </tr>\n'
        
        for [platform, compilername, compiler] in COMPILERS:
            s += '  <tr>'
            s += '<td><nobr>' + platform + ': ' + compilername + '</nobr></td>'
            for module in MODULES:
                task = TASKS[i]
                link = '<a href="' + PATH_RESULTFILES + task['filename'] + '.txt">'
                i += 1
                state = task['state']
                if state == 'pending': c = '<td class=pending align=middle>&nbsp;</td>'
                if state == 'started': 
                    c = '<td class=started align=middle>(' + task['server']['server'] + ')</td>'
                if state == 'error': c = '<td class=error align=middle>' + link + 'err</a></td>'
                if state == 'ok': c = '<td class=ok align=middle>' + link + 'ok</a></td>'
                s += c
                
            s += '</tr>\n'
            
        s += '</table>\n'
    
    s += '<div class=mode>Server loads</div>\n'
    s += '<table cellspacing=0 cellpadding=0 border=0>\n'
    for server in SERVERS:
        s += '<tr><td class=server>' + server['server'] + '</td><td>active=' + str(server['load']) + '</td><td>finished=' + str(server['count']) + '</td></tr>'
    s += '</table>'
    
    f = open(PATH_OUTPUT + FILE_OUTPUT)
    t = f.read()
    f.close()
    
    f = open(PATH_OUTPUT + FILE_OUTPUT_OLD, "wb")
    f.write(t)
    f.close()

    f = open(PATH_OUTPUT + FILE_OUTPUT, "wb")
    f.write(javaScriptEscape(s))
    f.close()

################################################################################

def showPage():
    global FILE_VIEWER, PATH_OUTPUT
    c = [PATH_OUTPUT + FILE_VIEWER]
    subprocess.Popen(c, shell=True)

################################################################################
    
def finishedAll():
    global TASKS
    
    for task in TASKS:
        if (task['state'] == 'pending') or (task['state'] == 'started'): return False
    return True

################################################################################

def scheduleTasks():
    global TASKS, SERVERS
    
    for server in SERVERS:
        if (server['load'] < server['loadmax']):
            startTaskForServer(server)
        
################################################################################
       
def pollTasks():
    global TASKS
    global ALL_STATES
    
    old_all_states = ALL_STATES
    ALL_STATES = ''

    for task in TASKS:
        ALL_STATES += task['state']
        
        if task['state'] != 'started': continue
        
        p = task['proc']
        if p == None: continue
        
        if (p.poll() != None): 
            onTaskStopped(task)
        
    return (ALL_STATES == old_all_states)


def onTaskStopped(task):
    proc = task['proc']
    out = proc.stdout.read()
    err = proc.stderr.read()
    
    is_error = (proc.poll() != 0)
    
    if is_error: task['state'] = 'error'
    else: task['state'] = 'ok'
    
    f = open(PATH_OUTPUT + PATH_RESULTFILES + task['filename'] + '.txt', "wb")
    f.write('module:   ' + task['module'] + '\n')
    f.write('platform: ' + task['platform'] + '\n')
    f.write('compiler: ' + task['compiler'] + '\n')
    f.write('mode:     ' + task['mode'] + '\n')
    f.write('server:   ' + task['server']['server'] + '\n')
    f.write('result:   ' + task['state'] + '\n')
    f.write('------------------------------------------------------------------\n')
    f.write('----------------------------- STDOUT -----------------------------\n')
    f.write('------------------------------------------------------------------\n')
    f.write(out)
    f.write('\n')
    if len(err) > 1:
        f.write('------------------------------------------------------------------\n')
        f.write('----------------------------- STDERR -----------------------------\n')
        f.write('------------------------------------------------------------------\n')
        f.write(err)
        f.write('\n')
    f.close()
    
    server = task['server']
    server['load'] -= 1
    server['count'] += 1
    
    print task['filename'] + ":",
    if is_error: print "ERROR"
    else: print "OK"

################################################################################

def startTaskForServer(server):
    global TASKS
    
    for task in TASKS:
        if (task['state'] == 'pending') and (task['platform'] == server['platform']):
            task['proc'] = startTask(task, server)
            server['load'] += 1
            task['state'] = 'started'
            task['server'] = server #server['i']
            return
            
    #no more tasks for that platform
    server['loadmax'] = server['load']
    
    
def startTask(task, server):
    platform = task['platform']
    if platform == 'gcc': return startTaskGCC(task, server)
    if platform == 'windows': return startTaskWindows(task, server)
    
################################################################################

def startTaskGCC(task, server):
    global USERNAME, PASSWORD
    
    c = ["plink.exe", "-batch", 
            "-ssh", "-l", USERNAME, "-pw", PASSWORD, server['server'], 
                "make", 
                    "-s", 
                    "-C\"~/MyDocuments/VELOP/seqan2/\"", 
                    "Platform=" + task['platform'],
                    "Project=" + task['module'],
                    "Mode=" + task['mode'],
                    "Compiler=" + task['compiler'],
                    "BuildFolder=" + task['platform'] + "/" + task['compilername']
         ]
#    print c;
    return subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        

def startTaskWindows(task, server):
    c = ["dir"]
    return subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

################################################################################


main()

