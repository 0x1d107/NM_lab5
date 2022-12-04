all: task1 task2 task3 task4
task1: task1.o generate.o
	$(CXX) -o $@ $^
task2: task2.o generate.o
	$(CXX) -o $@ $^
task3: task3.o generate.o
	$(CXX) -o $@ $^
task4: task4.o generate.o
	$(CXX) -o $@ $^
task5: task5.o generate.o
	$(CXX) -o $@ $^

%.o:%.cpp
	$(CXX) -c -o $@ $^
clean:
	rm -f task? *.o generate.o
.PHONY: clean
