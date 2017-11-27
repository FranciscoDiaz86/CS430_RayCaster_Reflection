all: raycast.c
	gcc raycast.c -o raycast -lm

clean:
	rm -rf raycast *~
