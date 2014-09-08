#
# bowtie2
#

GIT=$(shell which git)

bowtie2: bowtie2-git/bowtie2-align
	ln -s -f bowtie2-git/bowtie2 bowtie2
	ln -s -f bowtie2-git/bowtie2-align bowtie2-align

bowtie2-git/bowtie2-align: bowtie2-git
	cd bowtie2-git && make

bowtie2-git: $(GIT)
	$(GIT) clone git@github.com:BenLangmead/bowtie2.git
	mv bowtie2 $@

.PHONY: clean
clean:
	rm -f bowtie2 bowtie2-align
	rm -rf bowtie2-git
