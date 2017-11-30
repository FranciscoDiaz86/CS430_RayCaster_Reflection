all: raycast.c
	gcc raycast2.c -o raycast -lm

clean:
	rm -rf raycast *~
