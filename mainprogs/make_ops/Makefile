CHROMA=/home/jbulava/usqcd_soft/install/chroma_scalar_single
CONFIG=$(CHROMA)/bin/chroma-config
CXX=$(shell $(CONFIG) --cxx)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) #-I.
LDFLAGS=$(shell $(CONFIG) --ldflags)
LIBS=$(shell $(CONFIG) --libs)

HDRS=

OBJS= make_ops.o

make_ops: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf make_ops $(OBJS) *~
