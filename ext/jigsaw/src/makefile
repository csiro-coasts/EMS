#
# Basic JIGSAW static makefile
#
# NB: Not using any configure supplied tokens as only [newer] gcc
#     compilers are supported
#

CC=g++
ifdef DEBUG
 CFLAGS=-std=c++11 -pedantic -Wall -fPIC -g -D__lib_jigsaw -Ilibcpp
else
 CFLAGS=-std=c++11 -pedantic -Wall -O3 -flto -fPIC -DNDEBUG -D__lib_jigsaw -Ilibcpp
endif

# JIGSAW has only one source file
OBJS := jigsaw.o

JIGSAW_LIB := ../lib/libjigsaw64r.so
JIGSAW_EXE := ../bin/jigsaw64r

ALL_TARGETS := $(OBJS) $(JIGSAW_LIB) $(JIGSAW_EXE)
all : $(ALL_TARGETS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCS) -c $*.cpp -o $*.o

$(JIGSAW_LIB) : $(OBJS)
	$(CC) -shared -o $@ $^

$(JIGSAW_EXE) : jigsaw.cpp
	$(CC) -o $@ -static-libstdc++ -std=c++11 -pedantic -Wall -O3 -flto -DNDEBUG -Ilibcpp $^

clean:
	rm -f *.o *~ $(JIGSAW_LIB)

