#!gmake

.PHONY: all
all: qsim-parse-input qsim-parse-tandem

.PHONY: allall
all: qsim-parse-input qsim-parse-tandem qsim-parse-input-debug qsim-parse-tandem-debug

qsim-parse-input: qsim_parse_input.cpp
	g++ -O3 -o $@ $^

qsim-parse-input-debug: qsim_parse_input.cpp
	g++ -g -O0 -o $@ $^

qsim-parse-tandem: qsim_parse_tandem.cpp
	g++ -O3 -o $@ $^

qsim-parse-tandem-debug: qsim_parse_tandem.cpp
	g++ -g -O0 -o $@ $^

.PHONY: clean
clean:
	rm -f qsim-parse-input* qsim-parse-tandem* *.csv
