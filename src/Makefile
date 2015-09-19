#!gmake

.PHONY: all
all: qsim-parse-input

.PHONY: allall
all: qsim-parse-input qsim-parse-input-debug

PARSE_DEPS = qsim_parse_input.cpp simplesim.cpp input_model.cpp ranlib.cpp rnglib.cpp fasta.cpp

qsim-parse-input: $(PARSE_DEPS)
	g++ -O3 -o $@ $^

qsim-parse-input-debug: $(PARSE_DEPS)
	g++ -g -O0 -o $@ $^

.PHONY: clean
clean:
	rm -f qsim-parse-input* *.csv
