#___________________________________________________________________________
#
#  Makefile
#
#  author: 	Ivica Stevanovic		OFCOM
# ___________________________________________________________________________



PRGMAIN=test
PRGBIN=test

CC=g++

CFLAGS =   -O3 -m64 #-arch x86_64 

INCDIR =  -I.

INCLUDES := $(wildcard *.hpp)

all: $(PRGBIN)
	    
OBJECTS = P452.cpp \
          P452.test.cpp
	 

%.h:;

%.o: %.cpp $(INCLUDES)
	$(CC) $(CFLAGS) $(INCDIR) $(PFLAGS) -c $< -o $@ 


LCOM	= /usr/local/lib/libboost_unit_test_framework.a


$(PRGBIN): $(OBJECTS) $(PRGMAIN).o
	( $(CC) $(LFLAGS) $(OBJECTS) $(PRGMAIN).o $(LIBDIR) $(LIBS) $(LCOM) \
	      -o $(PRGBIN)   \
	) 


clean:
	rm -f *.o core
	
quick: 
	touch test.cpp	


