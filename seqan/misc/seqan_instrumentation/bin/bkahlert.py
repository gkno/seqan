#!/usr/bin/env python

import os
import sys
import threading
import binascii
from classes import *
sys.stdout = flushfile.Flushfile(sys.stdout)

class DiffCollector(object):
	def __init__(self, cMake_binary_dir, src_dir):
		self.cMake_binary_dir = cMake_binary_dir
		self.src_dir = src_dir
		
		self.dirs		= dirs.Dirs(self.src_dir, [ ".svn", "build", "util", "misc", "docs"""", "core", "extras""" ], [ "*.o", "Thumbs.db", ".DS_Store", "CMakeCache.txt" ])
		self.bin_dir		= self.dirs.get_abs_dir_path("misc/seqan_instrumentation/bin")
		self.last_revision_dir	= self.dirs.get_abs_dir_path("misc/seqan_instrumentation/last_revision_copy")
		self.userdata_dir	= self.dirs.get_abs_dir_path("misc/seqan_instrumentation/userdata")

		self.id_file		= self.dirs.get_abs_file_path("misc/seqan_instrumentation/userdata/id.txt")
		self.id			= id.ID(self.id_file)

		self.stats_file		= self.dirs.get_abs_file_path("misc/seqan_instrumentation/userdata/" + self.id.get() + "_stats.txt")
		self.stats		= stats.Stats(cMake_binary_dir, self.stats_file)
		self.stats.save("id", self.id.get())
		
	def getID(self):
		return self.id

	def getStatsFile(self):
		return self.stats_file

	def getStats(self):
		return self.stats

	def create_diff(self, differ):
		print("[1] Creating revision...");	
		next_diff_filename = differ.get_next_diff_filename()
		if os.path.isdir(self.last_revision_dir) and len(os.listdir(self.last_revision_dir)) > 0:
			differ.save_diff(next_diff_filename)
			print("[1] Revision created in " + next_diff_filename[len(self.src_dir):] + ".")
		else:
			differ.save_diff_reset(next_diff_filename)
			print("[1] Revision reset in " + next_diff_filename[len(self.src_dir):] + " due to non existing workspace comparison copy.")

	def upload_diff(self, differ):
		print("[2] Uploading diff files...");
		differ.sync_ftp("ftp.bkahlert.com", "u49635740-promo", binascii.unhexlify("684a753337685f5f5f6d"))
		print("[2] Diff files uploaded.")

	def update_comparision_copy(self, syncer):
		print("[3] Syncing comparison copy...");	
		syncer.make_comparison_copy(self.src_dir, self.last_revision_dir, self.dirs)
		print("[3] Comparison copy synced.")
	
	def prepare(self):
		if(not os.path.exists(self.last_revision_dir)):
			os.makedirs(self.last_revision_dir)
		if(not os.path.exists(self.userdata_dir)):
			os.makedirs(self.userdata_dir)
	
	def build(self):
		differ = diff.Diff(bin_dir=self.bin_dir, results_dir=self.userdata_dir, root_dir=self.src_dir, from_dir_local=self.last_revision_dir[len(self.src_dir):], to_dir_local=".", id=self.id, excluded_resources=self.dirs)
		syncer = sync.Sync(bin_dir=self.bin_dir)

		t1 = threading.Thread(target=self.create_diff, args=(differ,))
		t1.start()
		t1.join()

		t2 = threading.Thread(target=self.upload_diff, args=(differ,))
		t2.start()

		t3 = threading.Thread(target=self.update_comparision_copy, args=(syncer,))
		t3.start()
	
		t2.join()
		t3.join()
