SEQAN_BASE = ..

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

#CXXFLAGS += -m32
#CXXFLAGS += -m64
CXXFLAGS += -I$(SEQAN_BASE)
CXXFLAGS += -O9
#CXXFLAGS += -O0 -g -pg
CXXFLAGS += -pedantic -W -Wall
CXXFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 ${CXXEXTRAFLAGS}

TARGETS = razers/razers seqan_tcoffee/seqan_tcoffee seqcons/seqcons razers/paramChooser

all: check_seqan_base $(TARGETS)
razers: check_seqan_base razers/razers 
paramChooser: check_seqan_base razers/paramChooser 

check_seqan_base:
	@if [ ! -d "$(SEQAN_BASE)/seqan" ]; then \
		echo "The directory $(SEQAN_BASE)/seqan could not be found!"; \
		exit 1; \
	fi

clean:
	rm -f $(TARGETS)
