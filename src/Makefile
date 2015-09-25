#!gmake

.PHONY: all
all: qsim-parse qsim-rewrite

.PHONY: allall
allall: all qsim-parse-debug qsim-rewrite-debug

PARSE_DEPS = qsim_parse.cpp simplesim.cpp input_model.cpp ranlib.cpp rnglib.cpp fasta.cpp

REWRITE_DEPS = qsim_rewrite.cpp

qsim-parse: $(PARSE_DEPS)
	g++ -O3 -o $@ $^

qsim-parse-debug: $(PARSE_DEPS)
	g++ -g -O0 -o $@ $^

qsim-rewrite: $(REWRITE_DEPS)
	g++ -O3 -o $@ $^

qsim-rewrite-debug: $(REWRITE_DEPS)
	g++ -g -O0 -o $@ $^

.PHONY: clean
clean:
	rm -rf *.dSYM
	rm -f qsim-parse* *.csv qsim-rewrite*
