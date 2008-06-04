SEQAN = ..

# The following line is for the SVN repository
# You can remove it, if you use the release/snapshot
CPPFLAGS += -I$(SEQAN)/../../platforms/gcc

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

CPPFLAGS += -I$(SEQAN)
CPPFLAGS += -O9
CPPFLAGS += -pedantic
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

TARGETS = razers/razers seqan_tcoffee/seqan_tcoffee

all: $(TARGETS)
clean:
	rm -f $(TARGETS)
