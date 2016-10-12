#!gmake

TOOL=qtip

.PHONY: all
all: $(TOOL)-parse $(TOOL)-rewrite VERSION

.PHONY: allall
allall: all $(TOOL)-parse-debug $(TOOL)-rewrite-debug

PARSE_DEPS = $(TOOL)_parse.cpp simplesim.cpp input_model.cpp ranlib.cpp rnglib.cpp fasta.cpp

REWRITE_DEPS = $(TOOL)_rewrite.cpp

# http://programmers.stackexchange.com/questions/141973/how-do-you-achieve-a-numeric-versioning-scheme-with-git
# git tag -a v1.4.1 -m 'Version 1.4.1'
# git push --tags
VERSION:
	git describe --tags --long > $@

$(TOOL)-parse: $(PARSE_DEPS)
	g++ -O3 -o $@ $^

# note, on some JHU systems I have to use -gdwarf-3
$(TOOL)-parse-debug: $(PARSE_DEPS)
	g++ -g -O0 -o $@ $^

$(TOOL)-rewrite: $(REWRITE_DEPS)
	g++ -O3 -o $@ $^

# note, on some JHU systems I have to use -gdwarf-3
$(TOOL)-rewrite-debug: $(REWRITE_DEPS)
	g++ -g -O0 -o $@ $^

.PHONY: clean
clean:
	rm -rf *.dSYM
	rm -f $(TOOL)-parse* *.csv $(TOOL)-rewrite*
	rm -f VERSION
	rm -f *.pyc
