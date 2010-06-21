SEQAN_BASE = ..

# CXX = /home/takifugu2/data/teaching/mriese/bin/myGcc2/bin/g++

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

######################## Name ##############################################################

# 1.4.2 1.5.2 1.6
VERSION = 1

# tanne sequoia
SERVER = mac

FULL_NAME = razerMC-${VERSION}-${SERVER}
# FULL_NAME = razersWindow-${SERVER}

######################## My Flags ##########################################################
# Compile window version
# CXXFLAGS += -D RAZERS_WINDOW

# Compile parallel version choose one of the following versions: A, B, C
CXXFLAGS += -D RAZERS_PARALLEL_READS -D SEQAN_PARALLEL

# Compile version A (1.4.2): with which waits after each window
# CXXFLAGS += -D RAZERS_PARALLEL_READS_WAIT_AFTER_WINDOW
# VERSION = 1.4.2

# Compile verstion B (1.5.2): with flexible verification blocks
# CXXFLAGS += -D RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
# VERSION = 1.5.2

# Debug output for that option
# CXXFLAGS += -D FLEX_TIMER

# Sort and Partition hits with accuracy
# CXXFLAGS += -D RAZER_USE_ACCURACY

# Take time for sorting in Flex version
# CXXFLAGS += -D RAZER_SORT_TIMER

# Compile version C (1.6): that runs independent and does work stealing in the end
CXXFLAGS += -D RAZERS_PARALLEL_READS_INDEPENDENT
VERSION = 1.6.3

#################### General Flags ##########################################################
# Optimization
CXXFLAGS += -O3

# OpenMP
CXXFLAGS += -fopenmp

# Parallel sort
#CXXFLAGS += -D_GLIBCXX_PARALLEL

# Debug
# CXXFLAGS += -g -save-temps

# Library
CXXFLAGS += -I$(SEQAN_BASE)

# Seqan Flags
CXXFLAGS += -pedantic -W -Wall -Wno-variadic-macros
CXXFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 ${CXXEXTRAFLAGS}

razerMC:	check_seqan_base razers3/razers
	cp razers3/razers ~/masterKram/bin/${FULL_NAME}

test:	CXXFLAGS += -D SEQAN_ENABLE_TESTING
test:	check_seqan_base razers3/test_razers3/test_razers3
	./razers3/test_razers3/test_razers3

check_seqan_base:
	@if [ ! -d "$(SEQAN_BASE)/seqan" ]; then \
		echo "The directory $(SEQAN_BASE)/seqan could not be found!"; \
		exit 1; \
	fi

clean:
	rm -f $(TARGETS) $(TARGETS:=.o)
