DEBUG ?= 0
PARALLEL ?= 0

# Submodules
PWD = $(shell pwd)
SEQTK_ROOT = ${PWD}/src/htslib/htslib
KLIB_ROOT = ${PWD}/src/klib

# Flags
CXX=g++
GCC=gcc
CXXFLAGS += -isystem ${SEQTK_ROOT} -isystem ${KLIB_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas
LDFLAGS += -lz

# Additional flags for release/debug
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O9 -DNDEBUG
	#LDFLAGS += --static
endif

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SAMSOURCES = $(wildcard src/samtools/*.c) $(wildcard src/samtools/*.h)
BWASOURCES = $(wildcard src/bwa/*.c) $(wildcard src/bwa/*.h)
DELLYSOURCES = $(wildcard src/delly/src/*.h) $(wildcard src/delly/src/*.cpp)
BCSPLITSOURCES = $(wildcard src/*.h) $(wildcard src/*.c)

# Targets
TARGETS = .htslib .samtools .bwa .delly src/bcsplit 

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && cd ../../ && touch .htslib

.samtools: $(SAMSOURCES)
	cd src/samtools && make && cd ../../ && touch .samtools

.bwa: $(BWASOURCES)
	cd src/bwa && make && cd ../../ && touch .bwa

.delly: ${DELLYSOURCES}
	cd src/delly && make && cd ../../ && touch .delly

src/bcsplit: ${BCSPLITSOURCES}
	$(GCC) $(CXXFLAGS) $@.c -o $@ $(LDFLAGS)

clean:
	cd src/htslib && make clean
	cd src/samtools && make clean
	cd src/bwa && make clean
	cd src/delly && make clean
	rm -f $(TARGETS) $(TARGETS:=.o) .htslib .samtools .bwa .delly
