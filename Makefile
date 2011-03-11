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

CPPFLAGS=-Wall -O2 -march=nocona -std=c++0x
#CPPFLAGS=-Wall -g
LDFLAGS=-Wl,-s -static
#LDFLAGS=-Wl
PROGRAM=Spanner
LIBS=-lz -lboost_regex -lpthread

all: $(PROGRAM)
 
$(PROGRAM): $(OBJECTS) $(COBJECTS)
	@echo "- linking" $(PROGRAM)
	@$(CXX) $(LDFLAGS) $(FLAGS) -o $@ $^ $(LIBS) 

.PHONY: clean

clean:
	rm -f *.o $(PROGRAM) *~

