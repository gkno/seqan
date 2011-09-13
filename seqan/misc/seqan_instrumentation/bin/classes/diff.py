import os
import subprocess
import glob
from datetime import datetime
from ftplib import FTP
import socket

class Diff(object):
	MAX_REVISION_LENGTH = 8
        
	def __init__(self, bin_dir, results_dir, root_dir, from_dir_local, to_dir_local, id, excluded_resources):
		self.results_dir = results_dir
		self.root_dir = root_dir
		self.from_dir_local = from_dir_local
		self.to_dir_local = to_dir_local
		self.id = id
		self.excluded_resources = excluded_resources
		self.bin_dir = bin_dir

	def get_diff_filename(self, revision, datetime):
		return self.results_dir + "/" + self.id.get() + "_r" + str.rjust(str(revision), self.MAX_REVISION_LENGTH, "0") + "_" + datetime.strftime("%Y-%m-%dT%H-%M-%S") + ".diff"
	
	def get_diff_filenames(self):
		return glob.glob(self.results_dir + "/" + self.id.get() + "_r*_*.diff")

	def get_next_diff_filename(self):
		r = len(self.get_diff_filenames())
		if(r >= 2**self.MAX_REVISION_LENGTH):
			raise Exception("Maximum number of " + (2**self.MAX_REVISION_LENGTH) + " exceeded!")
		return self.get_diff_filename(r, datetime.today())

	def save_diff_reset(self, filename):
		f = open(filename, "w")
		f.write("RESET\n")
		f.close()

	def save_diff(self, filename):
		if(os.name == "nt"):
			p = subprocess.Popen([self.bin_dir + "/diffutils/bin/diff.exe", "-u", "-r", "-N"] + self.excluded_resources.get_exclude_list_diff() + ["." + self.from_dir_local, self.to_dir_local], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.root_dir)
		else:
			p = subprocess.Popen(["diff", "-u", "-r", "-N"] + self.excluded_resources.get_exclude_list_diff() + ["." + self.from_dir_local, self.to_dir_local], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.root_dir)	
		f = open(filename, "wb")
		output = p.communicate()[0]
		f.write(output)
		f.close()

	def upload(self, ftp, file):
		ftp.storbinary("STOR " + os.path.basename(file), open(file, "rb"), 1024)

	def sync_ftp(self, server, username, password):
		try:
			ftp = FTP(server, timeout=5)
			ftp.login(username, password)
		except socket.timeout:
			return

		# get remote file list
		remote_files = []
		try:
			remote_files = ftp.nlst()
		except:
			pass

		# upload all files that are not present on the ftp server
		for diff_filename in glob.glob(self.results_dir + "/" + self.id.get() + "_*"):
			if(os.path.basename(diff_filename) not in remote_files):
				self.upload(ftp, diff_filename)	
		ftp.quit()
