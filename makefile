all:
	g++ -std=c++14  main.cpp -o main.out
clean:
	rm -f *.out
run: clean all
	./main.out
