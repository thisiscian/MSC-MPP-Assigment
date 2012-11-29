##
##	The following variables are free to be modified
##  Folder paths should be given relative to the "src" folder
##  

EXECUTABLE_FOLDER = ../bin
PEXECUTABLE = parallel
SEXECUTABLE = serial
INPUT_FOLDER = ../input_files
OUTPUT_FOLDER = ../output_files
REQUIRED_FOLDERS = $(OUTPUT_FOLDER) $(OUTPUT_FOLDER)/logs $(INPUT_FOLDER) $(EXECUTABLE_FOLDER)

NTHREADS = 2
LIMIT_TYPE = 0
LIMIT = 100

##### These input files are included with the problem

INPUT_FILE=edge192x128.pgm
#INPUT_FILE=edge256x192.pgm 
#INPUT_FILE=edge512x384.pgm  
#INPUT_FILE=edge768x768.pgm
#INPUT_FILE=edge1024x1408.pgm 

CC = mpicc
CC_FLAGS = -lm
PSOURCE = main_parallel.c edgeview.c pgmio.c reconstruct.c parallel_reconstruct.c
SSOURCE = main_serial.c edgeview.c pgmio.c reconstruct.c
PHEADER = edgeview.h reconstruct.h parallel_reconstruct.h
SHEADER = edgeview.h reconstruct.h 

main:
	-@ mkdir $(REQUIRED_FOLDERS)
	make $(SEXECUTABLE)
	make $(PEXECUTABLE)

$(PEXECUTABLE): $(PSOURCE) Makefile $(PHEADER) | $(EXECUTABLE_FOLDER) $(REQUIRED_FOLDERS)
	$(CC) $(CC_FLAGS) $(PSOURCE) -o $(EXECUTABLE_FOLDER)/$@

$(SEXECUTABLE): $(SSOURCE) Makefile $(SHEADER) | $(EXECUTABLE_FOLDER) $(REQUIRED_FOLDERS)
	$(CC) $(CC_FLAGS) $(SSOURCE) -o $(EXECUTABLE_FOLDER)/$@

run: $(PEXECUTABLE)
	make srun
	make prun

srun: $(SEXECUTABLE)
	./$(EXECUTABLE_FOLDER)/$(SEXECUTABLE) $(LIMIT_TYPE) $(LIMIT) $(INPUT_FOLDER) $(INPUT_FILE) $(OUTPUT_FOLDER)

	
prun: $(PEXECUTABLE)
	mpiexec -n $(NTHREADS) $(EXECUTABLE_FOLDER)/$(PEXECUTABLE) $(LIMIT_TYPE) $(LIMIT) $(INPUT_FOLDER) $(INPUT_FILE) $(OUTPUT_FOLDER)


clean:
	-@ rm $(EXECUTABLE_FOLDER)/$(PEXECUTABLE)
	-@ rm $(EXECUTABLE_FOLDER)/$(SEXECUTABLE)
	-@ rm *.o  
	-@ rm $(OUTPUT_FOLDER)/*
	-@ rm $(OUTPUT_FOLDER)/logs/*
