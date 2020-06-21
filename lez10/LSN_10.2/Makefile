CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

EXECUTABLES := main2.exe
all : $(EXECUTABLES)

%.exe : %.o ../lab/ParRandGen/random.o
	$(CC) ../lab/ParRandGen/random.o $< -o $@
%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)
random.o : ../lab/ParRandGen/random.cpp ../lab/ParRandGen/random.h
	$(CC) -c ../lab/ParRandGen/random.cpp -o ../lab/ParRandGen/random.o $(CFLAGS)
clean :
	rm *.o *.exe seed.out
	

