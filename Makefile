all: clean compile generate run_serial run_parallel

clean:
	rm -f cout_serial.txt cout.txt cerr.txt cerr_serial.txt
compile:
	gcc -fopenmp rand_test_generator.c -o rand_test_generator.o -lm
	g++ -g -fopenmp main.cpp netlist.cpp batch.cpp Grid_Graph.cpp -o VLSI_Router.o
generate:
	./rand_test_generator.o 4
run_serial:
	./VLSI_Router.o sample_test.txt sampleout.txt 15 1 1 2> cerr_serial.txt 1> cout_serial.txt
run_parallel:
	./VLSI_Router.o sample_test.txt sampleout.txt 15 1 10 2> cerr.txt 1> cout.txt

default:
	rm -f cout_serial.txt cout.txt cerr.txt cerr_serial.txt

visualise: 
	python3 visualize.py sampleout.txt

evaluate:
	python3 Evaluator.py sample_test.txt sampleout.txt
