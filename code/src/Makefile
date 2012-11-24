CC = mpicc
CC_FLAGS = -lm
SOURCE = main.c edgeview.c pgmio.c reconstruct.c parallel_reconstruct.c
HEADER = edgeview.h reconstruct.h parallel_reconstruct.h
EXECUTABLE_FOLDER = ../bin
EXECUTABLE = casestudy

$(EXECUTABLE): $(SOURCE) Makefile $(HEADER)
	$(CC) $(CC_FLAGS) $(SOURCE) -o $(EXECUTABLE_FOLDER)/$@

run: $(EXECUTABLE)
	mpiexec -n 2 $(EXECUTABLE_FOLDER)/$(EXECUTABLE)

clean:
	-@ rm *.o  
	-@ rm $(EXECUTABLE_FOLDER)/$(EXECUTABLE)
