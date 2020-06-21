CC = g++
CFLAGS = -Wall -O3 --std=c++11

EXECUTABLES := main.exe main2.exe main3.exe
all : $(EXECUTABLES)

%.exe : %.o ../lez00/ParRandGen/random.o sigma.o
	$(CC) ../lez00/ParRandGen/random.o sigma.o $< -o $@
%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)
random.o : ../lez00/ParRandGen/random.cpp ../lez00/ParRandGen/random.h
	$(CC) -c ../lez00/ParRandGen/random.cpp -o ../lez00/ParRandGen/random.o $(CFLAGS)
sigma.o : sigma.cpp sigma.h
	$(CC) -c sigma.cpp -o sigma.o $(CFLAGS)
clean :
	rm *.o main.exe main2.exe main3.exe seed.out
	

