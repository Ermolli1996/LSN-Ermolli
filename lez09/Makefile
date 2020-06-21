CC = g++
CFLAGS = -Wall -O3 --std=c++11

EXECUTABLES := main.exe 
all : $(EXECUTABLES)

%.exe : %.o ../lez00/ParRandGen/random.o
	$(CC) ../lez00/ParRandGen/random.o $< -o $@
%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)
random.o : ../lez00/ParRandGen/random.cpp ../lez00/ParRandGen/random.h
	$(CC) -c ../lez00/ParRandGen/random.cpp -o ../lez00/ParRandGen/random.o $(CFLAGS)
clean :
	rm *.o *.exe seed.out
	

