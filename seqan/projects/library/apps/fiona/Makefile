SEQAN = ../..

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

CPPFLAGS += -I$(SEQAN)
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS += -O3 -g
CPPFLAGS += -pedantic


suffix: $(basename $(wildcard Fiona.cpp))

	
clean:
	rm -f $(basename $(wildcard FIONA.cpp))
	
