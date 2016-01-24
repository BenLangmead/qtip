#!gmake

.PHONY: all
all: qsim-parse qsim-rewrite VERSION

.PHONY: allall
allall: all qsim-parse-debug qsim-rewrite-debug

PARSE_DEPS = qsim_parse.cpp simplesim.cpp input_model.cpp ranlib.cpp rnglib.cpp fasta.cpp

REWRITE_DEPS = qsim_rewrite.cpp

# http://programmers.stackexchange.com/questions/141973/how-do-you-achieve-a-numeric-versioning-scheme-with-git
# git tag -a v1.0.4 -m 'Version 1.0.4'
# git push --tags
VERSION:
	git describe --tags --long > $@

qsim-parse: $(PARSE_DEPS)
	g++ -O3 -o $@ $^

# note, on some JHU systems I have to use -gdwarf-3
qsim-parse-debug: $(PARSE_DEPS)
	g++ -g -O0 -o $@ $^

qsim-rewrite: $(REWRITE_DEPS)
	g++ -O3 -o $@ $^

# note, on some JHU systems I have to use -gdwarf-3
qsim-rewrite-debug: $(REWRITE_DEPS)
	g++ -g -O0 -o $@ $^

.PHONY: clean
clean:
	rm -rf *.dSYM
	rm -f qsim-parse* *.csv qsim-rewrite*
	rm -f VERSION
