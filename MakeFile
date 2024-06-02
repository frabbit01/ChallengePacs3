CXX      = mpic++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -fopenmp -O3 -Wall -pedantic -I$(PACS_ROOT)/include -I$(PACS_ROOT)/src/Parallel/MPI/Matrix
LDFLAGS ?= -L$(PACS_ROOT)/lib
LIBS  ?= -lmpi -lmuparser 
LINK.o := $(LINK.cc)

EXEC = main
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o) 

DEPEND = make.dep

.PHONY = all $(EXEC) clean distclean $(DEPEND)

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(DEPEND)
	$(RM) *.o
doc:
	doxygen $(DOXYFILE)
	
distclean: clean
	$(RM) $(EXEC)
	$(RM) *.csv *.out *.bak *~
$(DEPEND): $(SRCS)
	@$(RM) $(DEPEND)
	@for file in $(SRCS); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

-include $(DEPEND)