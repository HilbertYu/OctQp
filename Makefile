###############################################################################
# compiler setting
###############################################################################
CC       = gcc
CXX      = g++
CFLAGS   = -g -Wall
CXXFLAGS = $(CFLAGS) -std=c++0x
LIBS     = -lm -loctave -loctinterp
INCPATH  = -I/usr/include/octave-4.0.0/ -I./ -I/usr/include/octave-4.0.0/octave
DIR     = $(shell pwd)

###############################################################################
# source files setting
###############################################################################
C_SOURCES   = $(shell find . -name "*.c")
CXX_SOURCES = $(shell find . -name "*.cpp")
C_OBJS      = $(patsubst %.c,%.o,$(wildcard $(C_SOURCES)))
CXX_OBJS    = $(patsubst %.cpp,%.o,$(wildcard $(CXX_SOURCES)))
OBJS        = $(C_OBJS) $(CXX_OBJS)
EXEC      = $(shell basename $(DIR))

###############################################################################
.PHONY : clean clean_all
###############################################################################
all: $(EXEC)

%.o:%.c
	$(CC) -c $(CFLAGS) $(INCPATH) $< -o $@
%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@ $(LIBS)

$(EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $(EXEC) $(LIBS)

###############################################################################
clean:
	@rm -vfr $(OBJS)
clean_all: clean
	@rm -vfr $(EXEC)
###############################################################################
