import os
import glob
import json
import socket
import sys
import platform
import shutil
from classes.configobj import ConfigObj

class Stats(object):
	CMakeCachePathKey = "CMakeCachePath"
        
	def __init__(self, cMake_binary_dir, src_dir, stats_cache):
		self.cMake_binary_dir = cMake_binary_dir
		self.src_dir = src_dir
		self.stats_cache = stats_cache
		self.stats = {}
		try:
			if os.path.isfile(self.stats_cache):
				self.load()
		except:
			pass

		self.get_os()
		self.get_machine()
		self.get_devenv()
		self.get_ip()

	def get(self, key):
		return self.stats[key]

	def set(self, key, value):
		if not key == None:
			self.stats[key] = value

	def load(self):
		f = open(self.stats_cache, "r")
		self.stats = json.load(f)
		f.close()

	def save(self, key=None, value=None):
                self.set(key, value)
                
		f = open(self.stats_cache, "w")
		json.dump(self.stats, f, indent = 4)
		f.close()

		cMakeCache_path = self.get_cMakeCache_path()
		if not cMakeCache_path is None and os.path.isfile(cMakeCache_path):
			shutil.copyfile(cMakeCache_path, self.get_cMakeCache_copy_path())
		
	def needs_refresh(self, stat):
		return not stat in self.stats or self.stats[stat] is None or len(self.stats[stat]) == 0
	
	def get_os(self):
		if(self.needs_refresh("os")):
			self.stats["os"] = { "sys.platform": sys.platform, "platform": platform.platform(), "system": platform.system(), "release": platform.release(), "version": platform.version(), "node": platform.node() }
			self.save()
		return self.stats["os"]

	def get_machine(self):
		if(self.needs_refresh("machine")):
			self.stats["machine"] = { "machine": platform.machine(), "processor": platform.processor(), "architecture": [ platform.architecture(), "64" if sys.maxsize > 2**32 else "no64" ]}
			self.save()
		return self.stats["machine"]

	def get_devenv(self):
		if(self.needs_refresh("devenv")):
			f = open(self.get_cMakeCache_path(), "r")
			lines = f.readlines()
			f.close()
			lines_without_slash_comments = map(lambda x: "" if len(x) >= 2 and x[0:2] == "//" else x, lines)

			config = ConfigObj(lines_without_slash_comments)
			
			self.stats["devenv"] = {
				"CMAKE_GENERATOR": config["CMAKE_GENERATOR:INTERNAL"] if "CMAKE_GENERATOR:INTERNAL" in config else None,
				"CMAKE_BUILD_TYPE": config["CMAKE_BUILD_TYPE:STRING"] if "CMAKE_BUILD_TYPE:STRING" in config else None,
				"CMAKE_C_COMPILER": config["CMAKE_GENERATOR:INTERNAL"] if "CMAKE_GENERATOR:INTERNAL" in config else None,
				"CMAKE_C_FLAGS": config["CMAKE_C_FLAGS:STRING"] if "CMAKE_C_FLAGS:STRING" in config else None,
				"CMAKE_CXX_COMPILER": config["CMAKE_CXX_COMPILER:FILEPATH"] if "CMAKE_CXX_COMPILER:FILEPATH" in config else None,
				"CMAKE_CXX_FLAGS": config["CMAKE_CXX_FLAGS:STRING"] if "CMAKE_CXX_FLAGS:STRING" in config else None,
				"CMAKE_LINKER": config["CMAKE_LINKER:FILEPATH"] if "CMAKE_LINKER:FILEPATH" in config else None,
				"CMAKE_MODULE_LINKER_FLAGS": config["CMAKE_MODULE_LINKER_FLAGS:STRING"] if "CMAKE_MODULE_LINKER_FLAGS:STRING" in config else None,
				"CMAKE_MAKE_PROGRAM": config["CMAKE_MAKE_PROGRAM:FILEPATH"] if "CMAKE_MAKE_PROGRAM:FILEPATH" in config else None
			}
			self.save()
		return self.stats["devenv"]		
	    
	def get_ip(self):
		if(self.needs_refresh("ip")):
			self.stats["ip"] = []
			for addrinfo in socket.getaddrinfo(socket.gethostname(), None, 0, 0, socket.SOL_TCP):
				self.stats["ip"].append(addrinfo[4][0])
			self.save()
		return self.stats["ip"]

	def get_cMakeCache_path(self):
        	cMakeCache_paths = map(os.path.normpath, glob.glob(self.cMake_binary_dir + "/../*/CMakeCache.txt"))
        	return None if len(cMakeCache_paths) == 0 else cMakeCache_paths[0]

	def get_cMakeCache_copy_path(self):
		(l, m, r) = self.stats_cache.rpartition(".")
		l = l + "_CMakeCache"
        	return os.path.normpath("".join([l, m, r]))
