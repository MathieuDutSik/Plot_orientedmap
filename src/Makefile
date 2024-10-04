PROGRAM = CombPlaneToSVG

OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++17 -Wall -g -O2

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL = -I../basic_common_cpp/src_basic -I../basic_common_cpp/src_number -I../basic_common_cpp/src_matrix -I../basic_common_cpp/src_comb -I../basic_common_cpp/src_graph -I../permutalib/src

CFLAGS = -I$(BOOST_INCDIR) -I$(EIGEN_PATH)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
