objects = main.o parallel.o 

all: $(objects)
		nvcc -arch=sm_20 $(objects) -o app

%.o: %.cpp
		g++ -m64 -O3 -Wall -g -c -o $@
%.o: %.cu
		nvcc -O3 -m64 cu -arch=compute_35 -I. -dc $< -o $@

clean:
		rm -f *.o app
