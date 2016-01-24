//
//  simplesim.cpp
//  qsim
//
//  Created by Ben Langmead on 9/12/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cctype>
#include <math.h>
#include "simplesim.h"
#include "fasta.h"
#include "rnglib.hpp"
#include "edit_xscript.h"

using namespace std;

/**
 * Return draw from binomial distribution with given n, p.
 */
static inline size_t draw_binomial(size_t n, float p) {
	return (size_t)ignbin((int)n, p);
}

/**
 * Mutate given simulated read in-place.
 */
void SimulatedRead::mutate(const char *seq) {
	const size_t newsz = strlen(qual_);
	while(newsz+1 >= seq_buf_len_) {
		double_seq_buf();
	}
	size_t rdoff = 0, rfoff = 0;
	const size_t edlen = strlen(edit_xscript_);
	for(size_t i = 0; i < edlen; i++) {
		if(edit_xscript_[i] == '=') {
			assert(seq[rfoff] != '\0');
			assert(isalpha(seq[rfoff]));
			seq_buf_[rdoff++] = seq[rfoff++];
		} else if(edit_xscript_[i] == 'X') {
			assert(seq[rfoff] != '\0');
			assert(isalpha(seq[rfoff]));
			do {
				seq_buf_[rdoff] = draw_base();
			} while(seq_buf_[rdoff] == seq[rfoff]);
			rdoff++;
			rfoff++;
		} else if(edit_xscript_[i] == 'I') {
			seq_buf_[rdoff++] = draw_base();
		} else if(edit_xscript_[i] == 'D') {
			rfoff++;
		} else if(edit_xscript_[i] == 'S') {
			seq_buf_[rdoff++] = draw_base();
			rfoff++;
		} else {
			throw 1;
		}
	}
	assert(rdoff == newsz);
	seq_buf_[rdoff] = '\0';
	if(newsz != strlen(seq_buf_)) {
		fprintf(stderr, "rdoff:%u, newsz:%u, strlen(qual_):%u, strlen(seq_buf_):%u, edit_xscript_:%s\n",
		        (unsigned)rdoff, (unsigned)newsz, (unsigned)strlen(qual_),
		        (unsigned)strlen(seq_buf_), edit_xscript_);
	}
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
	fprintf(fh, "@%s%s%s%s%c%s%llu%s%d%s%s\n",
			sim_startswith, sim_sep,
			refid_, sim_sep,
			fw_ ? '+' : '-', sim_sep,
			(unsigned long long)refoff_, sim_sep,
			score_, sim_sep,
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
 * Write pair of simulated reads to parallel FASTQ files.
 */
void SimulatedRead::write_pair(
	const SimulatedRead& rd1,
	const SimulatedRead& rd2,
	FILE *fh1,
	FILE *fh2,
	const char *typ)
{
	FILE *fhs[2] = {fh1, fh2};
	for(size_t i = 0; i < 2; i++) {
		fprintf(fhs[i], "@%s%s%s%s%c%s%llu%s%d%s%s%s%c%s%llu%s%d%s%s\n",
				sim_startswith, sim_sep,
				rd1.refid_, sim_sep,
				rd1.fw_ ? '+' : '-', sim_sep, // got different fws
				(unsigned long long)rd1.refoff_, sim_sep,
				rd1.score_, sim_sep,
				rd2.refid_, sim_sep,
				rd2.fw_ ? '+' : '-', sim_sep,
				(unsigned long long)rd2.refoff_, sim_sep,
				rd2.score_, sim_sep,
				typ);
		const SimulatedRead &rd = ((i == 0) ? rd1 : rd2);
		size_t len = strlen(rd.qual_);
		if(rd.fw_) {
			fputs(rd.seq_buf_, fhs[i]);
		} else {
			for(size_t j = 0; j < len; j++) {
				fputc(asc2dnacomp[(int)rd.seq_buf_[len-j-1]], fhs[i]);
			}
		}
		fputs("\n+\n", fhs[i]);
		if(rd.fw_) {
			fputs(rd.qual_, fhs[i]);
		} else {
			for(size_t j = 0; j < len; j++) {
				fputc(rd.qual_[len-j-1], fhs[i]);
			}
		}
		fputc('\n', fhs[i]);
	}
}

/**
 * Apply the specified function: max(minimum, factor * f(x)), where f might be
 * x or sqrt(x).
 */
size_t apply_function(float fraction, int function, size_t mn, size_t n) {
    if(n == 0) {
        return 0;
    }
    double nn = (double)n;
    if(function == FUNC_SQRT) {
        nn = sqrt(nn);
    }
    return std::max((size_t)(fraction * nn), mn);
}

/**
 * Simulate a batch of reads
 */
void StreamingSimulator::simulate_batch(
	float fraction,
	int function,
	size_t min_u,
	size_t min_c,
	size_t min_d,
	size_t min_b)
{
	size_t nc = 0, nd = 0, nu = 0, nb = 0;
	size_t n_wrote_c = 0, n_wrote_d = 0, n_wrote_u = 0, n_wrote_b = 0;
	nu = apply_function(fraction, function, min_u, model_u_.num_added());
	nb = apply_function(fraction, function, min_b, model_b_.num_added());
	nc = apply_function(fraction, function, min_c, model_c_.num_added());
	nd = apply_function(fraction, function, min_d, model_d_.num_added());
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
		if(retsz < olap_) {
			continue;
		}
		const size_t nchances = retsz - olap_ + 1;
		const float binom_p = min(((float)nchances) * 1.1f / tot_fasta_len_, 0.999f);
		memset(hist, 0, sizeof(int) * 256);
		for(size_t i = 0; i < retsz; i++) {
			hist[(int)buf[i]]++;
		}
		if(hist['N'] > (int)(0.9 * retsz)) {
			continue; // mostly Ns
		}
		
		// Maybe N content should affect choice for n*_chances
		
		//
		// Unpaired
		//
		
		size_t nu_samp = draw_binomial(nu, binom_p);
		for(size_t i = 0; i < nu_samp; i++) {
			int attempts = 0;
			do {
				if(attempts > max_attempts) {
					break;
				}
				attempts++;
				const TemplateUnpaired &t = model_u_.draw();
				size_t nslots = retsz - olap_;
				assert(nslots > 0);
				size_t off = std::min((size_t)(r4_uni_01() * nslots), nslots-1);
				assert(off < nslots);
				const size_t rflen = t.reflen();
				for(size_t j = off; j < off + rflen; j++) {
					const int b = buf[j];
					if(b != 'A' && b != 'C' && b != 'G' && b != 'T') {
						continue; // uses 1 attempt
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
				n_wrote_u++;
				rd1.write(fh_u_, "u");
			} while(false);
		}
		
		//
		// Bad-end
		//

		size_t nb_samp = draw_binomial(nb, binom_p);
		for(size_t i = 0; i < nb_samp; i++) {
			int attempts = 0;
			do {
				if(attempts > max_attempts) {
					break;
				}
				attempts++;
				const TemplateUnpaired &t = model_b_.draw();
				bool mate1 = t.mate_flag_ == '1';
				size_t nslots = retsz - olap_;
				size_t off = std::min((size_t)(r4_uni_01() * nslots), nslots-1);
				assert(off < nslots);
				const size_t rflen = t.reflen();
				for(size_t j = off; j < off + rflen; j++) {
					const int b = buf[j];
					if(b != 'A' && b != 'C' && b != 'G' && b != 'T') {
						continue; // uses 1 attempt
					}
				}
				if(mate1) {
					rd1.init(
						buf + off,
						t.qual_,
						t.edit_xscript_,
						t.fw_flag_ == 'T',
						t.best_score_,
						refid.c_str(),
						refoff + off);
					rd2.init_random(
						t.opp_len_,
						t.fw_flag_ == 'T', // doesn't matter much, but need them for name
						t.best_score_, // doesn't matter much, but need them for name
						refid.c_str(), // doesn't matter much, but need them for name
						refoff + off); // doesn't matter much, but need them for name
				} else {
					rd2.init(
						buf + off,
						t.qual_,
						t.edit_xscript_,
						t.fw_flag_ == 'T',
						t.best_score_,
						refid.c_str(),
						refoff + off);
					rd1.init_random(
						t.opp_len_,
						t.fw_flag_ == 'T', // doesn't matter much, but need them for name
						t.best_score_, // doesn't matter much, but need them for name
						refid.c_str(), // doesn't matter much, but need them for name
						refoff + off); // doesn't matter much, but need them for name
				}
				n_wrote_b++;
				const char *lab = mate1 ? "b1" : "b2";
				SimulatedRead::write_pair(rd1, rd2, fh_b_1_, fh_b_2_, lab);
			} while(false);
		}

		//
		// Concordant & discordant
		//
		
		size_t nc_samp = draw_binomial(nc, binom_p);
		size_t nd_samp = draw_binomial(nd, binom_p);
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
				size_t off = std::min((size_t)(r4_uni_01() * nslots), nslots-1);
				assert(off < nslots);
				size_t off_1, off_2;
				const size_t rflen_1 = edit_xscript_to_rflen(t.edit_xscript_1_);
				const size_t rflen_2 = edit_xscript_to_rflen(t.edit_xscript_2_);
				if(t.upstream1_) {
					off_1 = off;
					off_2 = off + std::max(t.fraglen_, rflen_2) - rflen_2;
				} else {
					off_2 = off;
					off_1 = off + std::max(t.fraglen_, rflen_1) - rflen_1;
				}
				for(size_t j = off_1; j < off_1 + rflen_1; j++) {
					const int b = buf[j];
					if(b != 'A' && b != 'C' && b != 'G' && b != 'T') {
						continue; // uses 1 attempt
					}
				}
				for(size_t j = off_2; j < off_2 + rflen_2; j++) {
					const int b = buf[j];
					if(b != 'A' && b != 'C' && b != 'G' && b != 'T') {
						continue; // uses 1 attempt
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
				if(conc) { n_wrote_c++; } else { n_wrote_d++; }
				const char *lab = conc ? "c" : "d";
				SimulatedRead::write_pair(rd1, rd2,
										  conc ? fh_c_1_ : fh_d_1_,
										  conc ? fh_c_2_ : fh_d_2_,
										  lab);
			} while(false);
		}
	}
	cerr << "    Wrote " << n_wrote_u << " unpaired tandem reads "
	     << "(target=" << nu << ")" << endl;
	cerr << "    Wrote " << n_wrote_b << " bad-end tandem reads "
	     << "(target=" << nb << ")" << endl;
	cerr << "    Wrote " << n_wrote_c << " concordant tandem pairs "
	     << "(target=" << nc << ")" << endl;
	cerr << "    Wrote " << n_wrote_d << " discordant tandem pairs "
	     << "(target=" << nd << ")" << endl;
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


