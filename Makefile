# ==================================
# define our source and object files
# ==================================

SOURCES=Spanner.cpp SpanDet.cpp RunControlParameterFile.cpp Function-Generic.cpp Function-Sequence.cpp Histo.cpp MosaikAlignment.cpp PairedData.cpp headerSpan.cpp cluster.cpp steps.cpp DepthCnvDet.cpp BamReader.cpp SHA1.cpp
OBJECTS=$(SOURCES:.cpp=.o)

CSOURCES=fastlz.c bgzf.c
COBJECTS=$(CSOURCES:.c=.o)

# ================
# compiler options
# ================
RE2DIR =/home/stewardg/Projects/source/re2
BAMTOOLS_ROOT=/home/stewardg/Projects/source/bamtools

CPPFLAGS=-Wall -O2 -march=nocona -std=c++0x -I$(RE2DIR) -I$(BAMTOOLS_ROOT)/include

#CPPFLAGS=-Wall -g
LDFLAGS=-Wl,-s -static
#LDFLAGS=-Wl
PROGRAM=Spanner
LIBS=-L$(RE2DIR)/obj/lib -L$(BAMTOOLS_ROOT)/lib -lz -lre2 -lbamtools -lpthread

all: $(PROGRAM)
 
$(PROGRAM): $(OBJECTS) $(COBJECTS)
	@echo "- linking" $(PROGRAM)
	@$(CXX) $(LDFLAGS) $(FLAGS) -o $@ $^ $(LIBS) 

.PHONY: clean

clean:
	rm -f *.o $(PROGRAM) *~
