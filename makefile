ALL:
	mpicc src/main.c ../util/Lab4_IO.c -o main -lm

main:
	mpicc src/main.c util/Lab4_IO.c -o main -lm

smain:
	gcc src/smain.c util/Lab4_IO.c -o smain -lm

dtrim:
	gcc util/datatrim.c -o datatrim

tester:
	gcc util/serialtester.c util/Lab4_IO.c -o serialtester -lm
