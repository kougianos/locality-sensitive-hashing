OBJS 	= main.o ATD.o LSHfunctions.o LSH.o
SOURCE	= main.c ATD.c LSHfunctions.c LSH.c
HEADER  = ATD.h LSHfunctions.h LSH.h
OUT  	= lsh
CC	= gcc
CFLAGS   = -g -c 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c

ATD.o: ATD.c
	$(CC) $(CFLAGS) ATD.c
	
LSHfunctions.o: LSHfunctions.c
	$(CC) $(CFLAGS) LSHfunctions.c

LSH.o: LSH.c
	$(CC) $(CFLAGS) LSH.c

clean:
	rm -f $(OBJS) $(OUT)

count:
	wc $(SOURCE) $(HEADER)