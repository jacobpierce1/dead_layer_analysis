ifeq ($(ROOTSYS),)
$(error ROOTSYS is not define, do ssetup root)
endif
ROOTCFLAGS := $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS := $(shell $(ROOTSYS)/bin/root-config --libs)
CXXFLAGS = -c -O2 $(ROOTCFLAGS)
LIBS = 

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: all clean

targets = Singles_Hist_Maker

all: $(targets)

Singles_Hist_Maker: Singles_Hist_Maker.o
	$(CXX) $^ -o $@ $(LIBS) $(ROOTLIBS)

clean:
	rm -f $(targets) *.o
