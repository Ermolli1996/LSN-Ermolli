CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

EXECUTABLES := main2.exe
all : $(EXECUTABLES)

%.exe : %.o ./ParRandGen/random.o
	$(CC) ./ParRandGen/random.o $< -o $@
%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)
random.o : ./ParRandGen/random.cpp ./ParRandGen/random.h
	$(CC) -c ./ParRandGen/random.cpp -o ./ParRandGen/random.o $(CFLAGS)
clean :
	rm *.o *.exe seed.out
	

