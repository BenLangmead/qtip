//
//  simplesim.cpp
//  qsim
//
//  Created by Ben Langmead on 9/12/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#include <stdlib.h>
#include <vector>
#include <inttypes.h>
#include "simplesim.h"
#include "fasta.h"
#include "rnglib.hpp"

using namespace std;

/**
 * Return draw from binomial distribution with given n, p.
 */
static inline size_t draw_binomial(size_t n, float p) {
	return (size_t)ignbin((int)n, p);
}

static const char * startswith = "!!ts!!";
static const char * sep = "!!ts-sep!!";

/**
 * Mutate given simulated read in-place.
 */
void SimulatedRead::mutate(const char *seq) {
	const size_t newsz = strlen(qual_);
	while(newsz >= seq_buf_len_) {
		double_seq_buf();
	}
	size_t rdoff = 0, rfoff = 0;
	const size_t edlen = strlen(edit_xscript_);
	for(size_t i = 0; i < edlen; i++) {
		if(edit_xscript_[i] == '=') {
			seq_buf_[rdoff++] = seq[rfoff++];
		} else if(edit_xscript_[i] == 'X') {
			do {
				seq_buf_[rdoff] = "ACGT"[(int)(r4_uni_01() * 4)];
			} while(seq_buf_[rdoff] == seq[rfoff]);
			rdoff++;
			rfoff++;
		} else if(edit_xscript_[i] == 'I') {
			seq_buf_[rdoff++] = "ACGT"[(int)(r4_uni_01() * 4)];
		} else if(edit_xscript_[i] == 'D') {
			rfoff++;
		} else {
			throw 1;
		}
	}
	assert(rdoff == newsz);
	seq_buf_[rdoff] = '\0';
	assert(newsz == strlen(seq_buf_));
}

/**
 * Mapping from ASCII characters for ambiguous nucleotides into masks:
 */
char asc2dnacomp[] = {
	/*   0 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/*  16 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/*  32 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,'-',  0,  0,
	/*  48 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/*  64 */ 0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
	/*    A   B   C   D           G   H           K       M   N */
	/*  80 */ 0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
	/*        R   S   T       V   W       Y */
	/*  96 */ 0,'T','V','G','H',  0,  0,'C','D',  0,  0,'M',  0,'K','N',  0,
	/*   a   b   c   d           g   h           k       m   n */
	/* 112 */ 0,  0,'Y','S','A',  0,'B','W',  0,'R',  0,  0,  0,  0,  0,  0,
	/*        r   s   t       v   w       y */
	/* 128 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/* 144 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/* 160 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/* 176 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/* 192 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/* 208 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/* 224 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	/* 240 */ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

/**
 * Write a simulated read to an output file.
 */
void SimulatedRead::write(FILE *fh, const char *typ) {
	size_t len = strlen(qual_);
	fprintf(fh, "@%s%s%s%s%c%s%" PRIuPTR "%s%d%s%s\n",
			startswith, sep,
			refid_, sep,
			fw_ ? '+' : '-', sep,
			(uintptr_t)refoff_, sep,
			score_, sep,
			typ);
	if(fw_) {
		fputs(seq_buf_, fh);
	} else {
		for(size_t i = 0; i < len; i++) {
			fputc(asc2dnacomp[(int)seq_buf_[len-i-1]], fh);
		}
	}
	fputs("\n+\n", fh);
	if(fw_) {
		fputs(qual_, fh);
	} else {
		for(size_t i = 0; i < len; i++) {
			fputc(qual_[len-i-1], fh);
		}
	}
	fputc('\n', fh);
}

/**
 * Simulate a batch of reads
 */
void StreamingSimulator::simulate_batch(
	float fraction,
	size_t min_u,
	size_t min_c,
	size_t min_d,
	size_t min_b)
{
	size_t nc = 0, nd = 0, nu = 0, nb = 0;
	if(!model_u_.empty()) {
		nu = std::max((size_t)(fraction * model_u_.num_added()), min_u);
	}
	if(!model_b_.empty()) {
		nb = std::max((size_t)(fraction * model_b_.num_added()), min_b);
	}
	if(!model_c_.empty()) {
		nc = std::max((size_t)(fraction * model_c_.num_added()), min_c);
	}
	if(!model_d_.empty()) {
		nd = std::max((size_t)(fraction * model_d_.num_added()), min_d);
	}
	assert(nu + nb + nc + nd > 0);
	
	std::string refid, refid_full;
	size_t refoff = 0, retsz = 0;
	SimulatedRead rd1, rd2;
	int hist[256];
	const int max_attempts = 10;
	while(true) {
		const char * buf = fa_.next(refid, refid_full, refoff, retsz);
		if(buf == NULL && fa_.done()) {
			break;
		}
		memset(hist, 0, sizeof(int) * 256);
		for(size_t i = 0; i < retsz; i++) {
			hist[(int)buf[i]]++;
		}
		if(hist['N'] > (int)(0.9 * retsz)) {
			continue; // mostly Ns
		}
		
		// Maybe N content should affect choice for n*_chances
		
		//
		// Unpaired & bad-end
		//
		
		size_t nu_chances = retsz - model_u_.avg_len() + 1;
		size_t nb_chances = retsz - model_b_.avg_len() + 1;
		size_t nu_samp = draw_binomial(nu, float(nu_chances) / tot_fasta_len_);
		size_t nb_samp = draw_binomial(nb, float(nb_chances) / tot_fasta_len_);
		for(size_t i = 0; i < nu_samp + nb_samp; i++) {
			bool unp = i < nu_samp;
			int attempts = 0;
			do {
				if(attempts > max_attempts) {
					break;
				}
				attempts++;
				const TemplateUnpaired &t = (unp ? model_u_.draw() : model_b_.draw());
				size_t nslots = retsz - olap_;
				size_t off = (size_t)(r4_uni_01() * nslots);
				assert(off < nslots);
				for(size_t j = off; j < off+t.len_; j++) {
					if(buf[j] != 'A' && buf[j] != 'C' &&
					   buf[j] != 'G' && buf[j] != 'T')
					{
						// have to worry about the case where the distribution
						// of Ns in the chunk is such that this loop iterates
						// many times or forever
						continue;
					}
				}
				rd1.init(
					buf + off,
					t.qual_,
					t.edit_xscript_,
					t.fw_flag_ == 'T',
					t.best_score_,
					refid.c_str(),
					refoff + off);
				const char *lab = unp ? "unp" : (t.mate_flag_ == '1' ?
												 "bad_end_mate1" :
												 "bad_end_mate2");
				rd1.write(unp ? fh_u_ : fh_b_, lab);
			} while(false);
		}
		
		//
		// Concordant & discordant
		//
		
		size_t nc_chances = retsz - model_c_.avg_len() + 1;
		size_t nd_chances = retsz - model_d_.avg_len() + 1;
		size_t nc_samp = draw_binomial(nc, float(nc_chances) / tot_fasta_len_);
		size_t nd_samp = draw_binomial(nd, float(nd_chances) / tot_fasta_len_);
		for(size_t i = 0; i < nc_samp + nd_samp; i++) {
			bool conc = i < nc_samp;
			int attempts = 0;
			do {
				if(attempts > max_attempts) {
					break;
				}
				attempts++;
				const TemplatePaired &t = conc ? model_c_.draw() : model_d_.draw();
				size_t nslots = retsz - olap_;
				size_t off = (size_t)(r4_uni_01() * nslots);
				assert(off < nslots);
				size_t off_1, off_2;
				if(t.upstream1_) {
					off_1 = off;
					off_2 = off + t.fraglen_ - t.len_2_;
				} else {
					off_2 = off;
					off_1 = off + t.fraglen_ - t.len_1_;
				}
				for(size_t j = off_1; j < off_1+t.len_1_; j++) {
					if(buf[j] != 'A' && buf[j] != 'C' &&
					   buf[j] != 'G' && buf[j] != 'T')
					{
						// have to worry about the case where the distribution
						// of Ns in the chunk is such that this loop iterates
						// many times or forever
						continue;
					}
				}
				for(size_t j = off_2; j < off_2+t.len_2_; j++) {
					if(buf[j] != 'A' && buf[j] != 'C' &&
					   buf[j] != 'G' && buf[j] != 'T')
					{
						// have to worry about the case where the distribution
						// of Ns in the chunk is such that this loop iterates
						// many times or forever
						continue;
					}
				}
				rd1.init(buf + off_1,
						 t.qual_1_,
						 t.edit_xscript_1_,
						 t.fw_flag_1_ == 'T',
						 t.score_1_,
						 refid.c_str(),
						 refoff + off_1);
				rd2.init(buf + off_2,
						 t.qual_2_,
						 t.edit_xscript_2_,
						 t.fw_flag_2_ == 'T',
						 t.score_2_,
						 refid.c_str(),
						 refoff + off_2);
				const char *lab = conc ? "conc" : "disc";
				rd1.write(conc ? fh_c_1_ : fh_d_1_, lab);
				rd2.write(conc ? fh_c_2_ : fh_d_2_, lab);
			} while(false);
		}
	}
}

#ifdef SIMPLESIM_MAIN

#include <fstream>
#include <iostream>

static void test1() {
	SimulatedRead rd;
	const char *ref = "ACGT";
	const char *qual = "ABCD";
	const char *edit_xscript = "====";
	rd.init(ref, qual, edit_xscript,true, 0, "r1", 0);
	assert(strcmp(rd.mutated_seq(), "ACGT") == 0);
	assert(strcmp(rd.qual(), "ABCD") == 0);
	assert(strcmp(rd.edit_xscript(), "====") == 0);
}

static void test2() {
	SimulatedRead rd;
	const char *ref = "AACC";
	const char *qual = "ABCD";
	const char *edit_xscript = "====";
	rd.init(ref, qual, edit_xscript,false, 0, "r1", 0);
	assert(strcmp(rd.mutated_seq(), "AACC") == 0);
	assert(strcmp(rd.qual(), "ABCD") == 0);
	assert(strcmp(rd.edit_xscript(), "====") == 0);
}

static void test3() {
	SimulatedRead rd;
	const char *ref = "ACGT";
	const char *qual = "ABCD";
	const char *edit_xscript = "=X==";
	rd.init(ref, qual, edit_xscript,true, 0, "r1", 0);
	assert(strcmp(rd.mutated_seq(), "ACGT") != 0);
	assert(rd.mutated_seq()[0] == 'A');
	assert(rd.mutated_seq()[1] != 'C');
	assert(rd.mutated_seq()[2] == 'G');
	assert(rd.mutated_seq()[3] == 'T');
	assert(strcmp(rd.qual(), "ABCD") == 0);
	assert(strcmp(rd.edit_xscript(), "=X==") == 0);
}

static void test4() {
	SimulatedRead rd;
	const char *ref = "ACGT";
	const char *qual = "ABC";
	const char *edit_xscript = "=D==";
	rd.init(ref, qual, edit_xscript,true, 0, "r1", 0);
	assert(strcmp(rd.mutated_seq(), "AGT") == 0);
	assert(strcmp(rd.qual(), "ABC") == 0);
	assert(strcmp(rd.edit_xscript(), "=D==") == 0);
}

static void test5() {
	SimulatedRead rd;
	const char *ref = "AGT";
	const char *qual = "ABCD";
	const char *edit_xscript = "=I==";
	rd.init(ref, qual, edit_xscript,true, 0, "r1", 0);
	assert(rd.mutated_seq()[0] == 'A');
	assert(rd.mutated_seq()[2] == 'G');
	assert(rd.mutated_seq()[3] == 'T');
	assert(strcmp(rd.qual(), "ABCD") == 0);
	assert(strcmp(rd.edit_xscript(), "=I==") == 0);
}

/**
 * A test where we read the FASTQ output and confirm it makes sense.
 */
static void test6() {
	SimulatedRead rd;
	const char *ref = "ACGT";
	const char *qual = "ABCD";
	const char *edit_xscript = "====";
	const char *fn = ".test6.tmp";
	{
		FILE *fh = fopen(fn, "wb");
		rd.init(ref, qual, edit_xscript, true, 0, "r1", 0);
		rd.write(fh, "hello");
		fclose(fh);
	}
	ifstream t(fn);
	string str((istreambuf_iterator<char>(t)), istreambuf_iterator<char>());
	assert(str.rfind("\nACGT\n+\nABCD\n") != string::npos);
	remove(fn);
}

static void test7() {
	SimulatedRead rd;
	const char *ref = "AAACC";
	const char *qual = "EDCBA";
	const char *edit_xscript = "=====";
	const char *fn = ".test7.tmp";
	{
		FILE *fh = fopen(fn, "wb");
		rd.init(ref, qual, edit_xscript, false, 0, "r1", 0);
		rd.write(fh, "hello");
		fclose(fh);
	}
	ifstream t(fn);
	string str((istreambuf_iterator<char>(t)), istreambuf_iterator<char>());
	assert(str.rfind("\nGGTTT\n+\nABCDE\n") != string::npos);
	remove(fn);
}

int main(void) {
	test1();
	test2();
	test3();
	test4();
	test5();
	test6();
	test7();
	cerr << "ALL TESTS PASSED" << endl;
}
#endif


