#___________________________________________________________________________
#
#  Makefile
#
#  author: 	Ivica Stevanovic	OFCOM
# ___________________________________________________________________________



PRGMAIN=test
PRGBIN=test

CC=g++ -O3


CFLAGS =   -m64 #-arch x86_64 

INCDIR =  -I.

INCLUDES := $(wildcard *.hpp)

all: $(PRGBIN)
	    
OBJECTS = P452.cpp \
          P452.test.cpp
	 

%.h:;

%.o: %.cpp $(INCLUDES)
	$(CC) $(CFLAGS) $(INCDIR) $(PFLAGS) -c $< -o $@ 

LCOM = -L/usr/local/lib -lboost_unit_test_framework

$(PRGBIN): $(OBJECTS) $(PRGMAIN).o
	( $(CC) $(LFLAGS) $(OBJECTS) $(PRGMAIN).o $(LIBDIR) $(LIBS) $(LCOM) \
	      -o $(PRGBIN)   \
	) 


clean:
	rm -f *.o core
	
quick: 
	touch test.cpp	


