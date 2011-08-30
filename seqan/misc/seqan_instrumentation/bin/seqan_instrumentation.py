#!/usr/bin/env python

import sys
from bkahlert import DiffCollector

def main():
	# sys.argv[1]: event name
	# sys.argv[2]: cmake binary dir (e.g. [SEQAN]/build/Debug)
	# sys.argv[3]: seqan src directory
	diffCollector = DiffCollector(sys.argv[2], sys.argv[3])
        
	# cmake run
	if(sys.argv[1] == "cmake"):
		diffCollector.prepare()
		print("[NOTE] Your ID is " + diffCollector.getID().get())
		print("[NOTE] Your documented data are saved in file\n       " + diffCollector.getStatsFile() + ".")
		diffCollector.build()
		return 0

        # build (only called no matter how many target need to be rebuild)
	if(sys.argv[1] == "build"):
		diffCollector.prepare()
		diffCollector.build()
		return 0

	# pre build called before each target build
	# sys.argv[4]: target name
	if(sys.argv[1] == "pre_build"):
		return 0

	# post build called after each target build
	# sys.argv[4]: target name
	if(sys.argv[1] == "post_build"):
		return 0

	return 0

if __name__ == '__main__':
	sys.exit(main())

