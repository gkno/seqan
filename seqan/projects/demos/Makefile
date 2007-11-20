SEQAN_DIR=..

# The following 4 lines are for the SVN repository
# You can remove them, if you use the snapshot
ifeq ($(shell test -d ../library && echo Y),Y)
  SEQAN_DIR=../library
  CPPFLAGS+=-I../../platforms/gcc
endif

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

CPPFLAGS+=-I$(SEQAN_DIR)
#CPPFLAGS+=-O3

all: $(basename $(wildcard *.cpp))
index: $(basename $(wildcard index_*.cpp))
graph: $(basename $(wildcard graph_*.cpp))
modifier: $(basename $(wildcard modifier_*.cpp))
file: $(basename $(wildcard file_*.cpp))
clean:
	rm -f $(basename $(wildcard *.cpp))
