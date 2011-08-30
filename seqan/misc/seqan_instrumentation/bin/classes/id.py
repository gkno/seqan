import os
import random
import string

"""Generates a permanent ID."""
class ID(object):
	"""
	@param cache_file is the file the generated ID is stored in
	"""
	def __init__(self, cache_file):
		self.length = 16
		self.cache_file = cache_file

	def get(self):
		if(os.path.isfile(self.cache_file)):
			f = open(self.cache_file, "r")
			id = f.readlines()[0]
			f.close
			if(len(id) == self.length): return id
		
		id = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(self.length))

		self.create_dir()
		
		f = open(self.cache_file, "w")
		f.write(id)
		f.close
		return id

	def __repr__(self):
		return get()

	def create_dir(self):
		try:
			os.makedirs(os.path.dirname(self.cache_file))
		except:
			pass
