CC = mpicc
CC_FLAGS = -lm
SOURCE = main.c edgeview.c pgmio.c reconstruct.c parallel_reconstruct.c
HEADER = edgeview.h reconstruct.h parallel_reconstruct.h
EXECUTABLE_FOLDER = ../bin
EXECUTABLE = image_reconstruction_2d

$(EXECUTABLE): $(SOURCE) Makefile $(HEADER) | $(EXECUTABLE_FOLDER) ../output_files ../output_files/logs
	$(CC) $(CC_FLAGS) $(SOURCE) -o $(EXECUTABLE_FOLDER)/$@

run: $(EXECUTABLE)
	mpiexec -n 8 $(EXECUTABLE_FOLDER)/$(EXECUTABLE)

clean:
	-@ rm *.o  
	-@ rm $(EXECUTABLE_FOLDER)/$(EXECUTABLE)
