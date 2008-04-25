# The following line is for the SVN repository
# You can remove it, if you use the release/snapshot
CPPFLAGS += -I../../../../platforms/gcc

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

CPPFLAGS += -I../..
CPPFLAGS += -O3
CPPFLAGS += -pedantic
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

pizzachili: pizzachili.cpp ../../lib/libpizzachili.a
	$(CXX) $(CPPFLAGS) pizzachili.cpp -o $@ -L../../lib -lpizzachili

../../lib/libpizzachili.a:
	$(MAKE) -C ../../extra/pizzachili
	
clean:
	rm -f $(basename $(wildcard *.cpp))
