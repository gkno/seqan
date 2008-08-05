SEQAN_BASE = ..

# The following line is for the SVN repository
# You can remove it, if you use the release/snapshot
CXXFLAGS += -I$(SEQAN_BASE)/../../platforms/gcc

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

#CXXFLAGS += -m32
#CXXFLAGS += -m64
CXXFLAGS += -I$(SEQAN_BASE)
CXXFLAGS += -O9
CXXFLAGS += -pedantic
CXXFLAGS += -W -Wall
CXXFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

#QGRAM_LENS = 7 8 9 10 11 12 13 14
QGRAM_LENS = 8 11 12 13
RAZER_TARGETS := $(patsubst %,razers/razers.q%,$(QGRAM_LENS))
TARGETS = razers/razers $(RAZER_TARGETS) seqan_tcoffee/seqan_tcoffee

all: check_seqan_base $(TARGETS)
razers: check_seqan_base razers/razers 

razers/razers.q%: check_seqan_base razers/razers.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -DQGRAM_LEN=$(patsubst razers/razers.q%,%,$@) razers/razers.cpp -o $@

check_seqan_base:
	@if [ ! -d "$(SEQAN_BASE)/seqan" ]; then \
		echo "The directory $(SEQAN_BASE)/seqan could not be found!"; \
		exit 1; \
	fi

clean:
	rm -f $(TARGETS)
