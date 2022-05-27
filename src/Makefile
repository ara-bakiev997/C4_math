CC=gcc
CFLAGS=-std=c11 -Wall -Werror -Wextra -g
SOURCE=s21*.c
OBJECTS=$(SOURCE:.c=.o)
GCOVFLAGS=-L. --coverage
HTML=lcov -t "test" -o rep.info -c -d .
OS=$(shell uname)

ifeq ($(OS), Darwin)
	CHECK_H=-lcheck
else
	CHECK_H=-lcheck -pthread -lrt -lm
endif

all: gcov_report

s21_math.a: clean $(OBJECTS)
	ar rcs s21_math.a $(OBJECTS)

test: s21_math.a
	$(CC) $(CFLAGS) math_test.c $(SOURCE) -o test $(CHECK_H) $(GCOVFLAGS)
	./test

gcov_report: test
	$(HTML)
	genhtml -o report rep.info

check:
	cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config $(SOURCE) *.h
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extensions=c $(SOURCE) *.h
	rm -rf CPPLINT.cfg
	make test
ifeq ($(OS), Darwin)
	leaks --atExit -- test
else
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./test
endif

clean:
	rm -rf *.o *.a *.so *.gcda *.gcno *.gch rep.info *.html *.css test report *.txt