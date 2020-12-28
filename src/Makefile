PROGRAM = CombPlaneToSVG


OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++11
CC		= g++ -std=c++11 -Wall -g -O2

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL = -I../basic_common_cpp/src_basic -I../basic_common_cpp/src_number -I../basic_common_cpp/src_matrix -I../basic_common_cpp/src_comb -I../basic_common_cpp/src_graph

CFLAGS = -I$(BOOST_INCDIR) -I$(EIGEN_PATH) -I$(GMP_INCDIR) -I$(PERMLIB_INCDIR)

LIBS = $(GMP_CXX_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
