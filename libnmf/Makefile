# Makefile for libnmf
# 
# make libnmf.a - create static library
# make libnmf.so - create dynamic library

CCFLAGS = -Wall -Wextra -Iinclude/ -g -fPIC
LIBS = -llapack -larpack

SOURCE = $(wildcard *.c)
OBJ = $(SOURCE:%.c=%.o)

.INTERMEDIATE: $(OBJ)

all: libnmf.a libnmf.so

libnmf.a: $(OBJ)
	ar rcs $@ $^

libnmf.so: $(OBJ)
	$(CC) -shared -Wl,-soname,$@ $(LIBS) -o $@ $^

nmf_mu.so: nmf_mu.o
	$(CC) -shared -Wl,-soname,$@ $(LIBS) -o $@ $^

%.o: %.c include/%.h
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -f libnnmf.so libnnmf.a $(OBJ)

