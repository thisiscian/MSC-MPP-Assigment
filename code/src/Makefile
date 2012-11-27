CC = mpicc
CC_FLAGS = -lm
PSOURCE = main_parallel.c edgeview.c pgmio.c reconstruct.c parallel_reconstruct.c
SSOURCE = main_serial.c edgeview.c pgmio.c reconstruct.c
PHEADER = edgeview.h reconstruct.h parallel_reconstruct.h
SHEADER = edgeview.h reconstruct.h 
EXECUTABLE_FOLDER = ../bin
PEXECUTABLE = parallel
SEXECUTABLE = serial
ITERATIONS = 100
INPUT_FOLDER = ../input_files
OUTPUT_FOLDER = ../output_files
INPUT_FILE = edge192x128.pgm
REQUIRED_FOLDERS = $(OUTPUT_FOLDER) $(OUTPUT_FOLDER)/logs
NTHREADS = 12

main:
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
	./$(EXECUTABLE_FOLDER)/$(SEXECUTABLE) $(ITERATIONS) $(INPUT_FOLDER) $(INPUT_FILE) $(OUTPUT_FOLDER)

	
prun: $(PEXECUTABLE)
	mpiexec -n $(NTHREADS) $(EXECUTABLE_FOLDER)/$(PEXECUTABLE) $(ITERATIONS) $(INPUT_FOLDER) $(INPUT_FILE) $(OUTPUT_FOLDER)


clean:
	-@ rm $(EXECUTABLE_FOLDER)/$(PEXECUTABLE)
	-@ rm $(EXECUTABLE_FOLDER)/$(SEXECUTABLE)
	-@ rm *.o  
	-@ rm $(OUTPUT_FOLDER)/*
	-@ rm $(OUTPUT_FOLDER)/logs/*
