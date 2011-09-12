from distutils.core import setup
import sys
import py2exe
from glob import glob
import os

# ISSUE: import classes (via includes) instead of copying it (via data_files)

sys.path.append("C:\\Program Files\\CMake 2.8\\bin")
sys.path.append(os.path.normpath(os.getcwd() + "/../bin"))
setup(options={
    "py2exe": { "includes": ["json", "ftplib", "platform", "glob", "shutil", "webbrowser"], "skip_archive": True } },
    data_files=[("Microsoft.VC90.CRT", glob(r'.\Microsoft.VC90.CRT\*.*')), ("classes", glob(r'..\bin\classes\*.*'))],
    console=["../bin/seqan_instrumentation.py"])
