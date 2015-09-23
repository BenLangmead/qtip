#!gmake

.PHONY: all
all: qsim-parse

.PHONY: allall
all: qsim-parse qsim-parse-debug

PARSE_DEPS = qsim_parse.cpp simplesim.cpp input_model.cpp ranlib.cpp rnglib.cpp fasta.cpp

qsim-parse: $(PARSE_DEPS)
	g++ -O3 -o $@ $^

qsim-parse-debug: $(PARSE_DEPS)
	g++ -g -O0 -o $@ $^

.PHONY: clean
clean:
	rm -f qsim-parse* *.csv
