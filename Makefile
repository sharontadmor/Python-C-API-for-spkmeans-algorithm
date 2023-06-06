CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm

spkmeans: spkmeans_main.o spkmeans.o utils.o
	$(CC) -o spkmeans spkmeans_main.o spkmeans.o utils.o $(CFLAGS)

spkmeans_main.o: spkmeans_main.c spkmeans.h
	$(CC) -c spkmeans_main.c $(CFLAGS)

spkmeans.o: spkmeans.c spkmeans.h
	$(CC) -c spkmeans.c $(CFLAGS)

utils.o: utils.c spkmeans.h
	$(CC) -c utils.c $(CFLAGS)

all: spkmeans

spkmeans_asan: spkmeans_main_asan.o spkmeans_asan.o utils_asan.o
	$(CC) -o $@ spkmeans_main_asan.o spkmeans_asan.o utils_asan.o $(CFLAGS) -fsanitize=address -static-libasan -g

spkmeans_main_asan.o: spkmeans_main.c spkmeans.h
	$(CC) -o $@ -c spkmeans_main.c $(CFLAGS) -fsanitize=address -static-libasan -g

spkmeans_asan.o: spkmeans.c spkmeans.h
	$(CC) -o $@ -c spkmeans.c $(CFLAGS) -fsanitize=address -static-libasan -g

utils_asan.o: utils.c spkmeans.h
	$(CC) -o $@ -c utils.c $(CFLAGS) -fsanitize=address -static-libasan -g

clean: 
	rm -f spkmeans spkmeans_main.o spkmeans.o utils.o spkmeans_asan spkmeans_main_asan.o spkmeans_asan.o utils_asan.o
