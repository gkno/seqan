SEQAN_BASE = ..

# The following line is for the SVN repository
# You can remove it, if you use the release/snapshot
CPPFLAGS += -I$(SEQAN_BASE)/../../platforms/gcc

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

CPPFLAGS += -I$(SEQAN_BASE)
CPPFLAGS += -O9
CPPFLAGS += -pedantic
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

TARGETS = razers/razers seqan_tcoffee/seqan_tcoffee

all: check_seqan_base $(TARGETS)
razers: check_seqan_base razers/razers
seqan_tcoffee: check_seqan_base seqan_tcoffee/seqan_tcoffee

check_seqan_base:
	@if [ ! -d "$(SEQAN_BASE)/seqan" ]; then \
		echo "The directory $(SEQAN_BASE)/seqan could not be found!"; \
		exit 1; \
	fi

clean:
	rm -f $(TARGETS)
