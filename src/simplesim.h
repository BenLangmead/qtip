//
//  simplesim.h
//  qsim
//
//  Created by Ben Langmead on 9/12/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#ifndef __qsim__simplesim__
#define __qsim__simplesim__

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include "fasta.h"
#include "input_model.h"
#include "ranlib.hpp"

#define SIM_STARTSWITH_LITERAL "qsim!"
#define SIM_SEPARATOR_LITERAL ':'

static const char * sim_startswith = SIM_STARTSWITH_LITERAL;
static const char sim_sep = SIM_SEPARATOR_LITERAL;

static inline char draw_base() {
    return "ACGT"[std::min((int)(r4_uni_01() * 4), 3)];
}

/**
 * The thing is really built in two stages.  First we get seq (but it's pointer
 * into the reference buffer), fw, qual (but it's a pointer into the template
 * object), edit_xscript (also pointer to template), score, refid (from string
 * variable in simulator), and refoff (from size_t in simulator).
 *
 * Then we mutate the reference sequence into the appropriate read sequence
 * using the edit transcript.  When we do this, we have to allocate a new
 */
class SimulatedRead {

public:
	
	SimulatedRead() :
		fw_(),
		qual_(NULL),
		edit_xscript_(NULL),
		score_(),
		refid_(NULL),
		refoff_()
	{
		seq_buf_len_ = 64;
		seq_buf_ = new char[seq_buf_len_];
		qual_buf_len_ = 64;
		qual_buf_ = new char[qual_buf_len_];
	}
	
	~SimulatedRead() {
		if(seq_buf_ != NULL) {
			delete[] seq_buf_;
			seq_buf_ = NULL;
		}
		if(qual_buf_ != NULL) {
			delete[] qual_buf_;
			qual_buf_ = NULL;
		}
	}
	
	void init(
		const char *seq,
		char *qual,
		char *edit_xscript,
		bool fw,
		int score,
		const char *refid,
		size_t refoff)
	{
		qual_ = qual;
		edit_xscript_ = edit_xscript;
		fw_ = fw;
		score_ = score;
		refid_ = refid;
		refoff_ = refoff;
		mutate(seq);
	}
	
	/**
	 * Initialize with random bases and sequence string.
	 */
	void init_random(
		size_t len,
		bool fw,
		int score,
		const char *refid,
		size_t refoff)
	{
		assert(len > 0);
		while(len+1 >= qual_buf_len_) {
			double_qual_buf();
		}
		while(len+1 >= seq_buf_len_) {
			double_seq_buf();
		}
		fw_ = fw;
		score_ = score;
		refid_ = refid;
		refoff_ = refoff;
		char *seq_cur = seq_buf_;
		char *qual_cur = qual_buf_;
		for(size_t i = 0; i < len; i++) {
			int c = draw_base();
			int q = 'I';
			*seq_cur++ = c;
			*qual_cur++ = q;
		}
		*seq_cur = *qual_cur = '\0';
		qual_ = qual_buf_;
	}
	
	/**
	 * Write unpaired simulated read to a FASTQ file.
	 */
	void write(FILE *fh, const char *typ);
	
	/**
	 * Write pair of simulated reads to parallel FASTQ files.
	 */
	static void write_pair(
		const SimulatedRead& rd1,
		const SimulatedRead& rd2,
		FILE *fh1,
		FILE *fh2,
		const char *typ);
	
	/**
	 * Return read sequence, as mutated from reference.
	 */
	const char *mutated_seq() const {
		return seq_buf_;
	}
	
	/**
	 * Return quality string.
	 */
	const char *qual() const {
		return qual_;
	}

	/**
	 * Return edit transcript.
	 */
	const char *edit_xscript() const {
		return edit_xscript_;
	}

protected:
	
	/**
	 * Mutate simulated read in-place using edit transcript and reference
	 * string.
	 *
	 * Note: edit transcript is always presented as though 5' end of read is on the left?
	 */
	void mutate(const char *seq);

	/**
	 * Double the size of the sequence buffer.
	 */
	void double_seq_buf() {
		delete[] seq_buf_;
		seq_buf_len_ *= 2;
		seq_buf_ = new char[seq_buf_len_];
	}

	/**
	 * Double the size of the quality buffer.
	 */
	void double_qual_buf() {
		delete[] qual_buf_;
		qual_buf_len_ *= 2;
		qual_buf_ = new char[qual_buf_len_];
	}

	bool fw_;
	char *qual_;
	const char *edit_xscript_;
	int score_;
	const char *refid_;
	size_t refoff_;
	
	size_t seq_buf_len_;
	char *seq_buf_;

	size_t qual_buf_len_;
	char *qual_buf_;
};

enum {
    FUNC_LINEAR = 1,
    FUNC_SQRT = 2
};

/**
 * What do we need from the dists?
 * 1. Average read/fragment lengths for all 4 classes
 * 2. Maximum read/fragment length across all 4 classes
 * 3. # alignments encountered of each of the 4 types
 * 4. A way of randomly drawing templates of the 4 types
 */
class StreamingSimulator {
	
public:
	
	StreamingSimulator(
		const std::vector<std::string>& fns, // FASTAs to simulate reads from
		size_t chunksz, // rolling FASTA buffer length
		const InputModelUnpaired& model_u,
		const InputModelUnpaired& model_b,
		const InputModelPaired& model_c,
		const InputModelPaired& model_d,
		FILE *fh_u,
		FILE *fh_b_1,
		FILE *fh_b_2,
		FILE *fh_c_1,
		FILE *fh_c_2,
		FILE *fh_d_1,
		FILE *fh_d_2) :
		olap_(std::max(model_u.max_len(),
			  std::max(model_b.max_len(),
			  std::max(model_c.max_len(), model_d.max_len())))),
		fa_(fns, chunksz, olap_),
		model_u_(model_u),
		model_b_(model_b),
		model_c_(model_c),
		model_d_(model_d),
		fh_u_(fh_u),
		fh_b_1_(fh_b_1),
		fh_b_2_(fh_b_2),
		fh_c_1_(fh_c_1),
		fh_c_2_(fh_c_2),
		fh_d_1_(fh_d_1),
		fh_d_2_(fh_d_2)
	{
		tot_fasta_len_ = estimate_fasta_length(fns);
	}
	
	/**
	 * Simulate a batch of reads
	 */
	void simulate_batch(
		float fraction,
		int function,
		size_t min_u,
		size_t min_c,
		size_t min_d,
		size_t min_b);
	
	/**
	 * Return the estimated number of bases in all the FASTA files, based on
	 * the file sizes.
	 */
	size_t num_estimated_bases() const {
		return tot_fasta_len_;
	}

protected:
	
	/**
	 * Return size of file in bytes.
	 */
	static size_t filesize(const std::string& filename) {
		std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
		return (size_t)in.tellg();
	}
	
	/**
	 * Given list of FASTA file names, estimate total size.
	 */
	static size_t estimate_fasta_length(const std::vector<std::string>& fns) {
		size_t tot = 0;
		for(size_t i = 0; i < fns.size(); i++) {
			tot += filesize(fns[i]);
		}
		return tot;
	}
	
	size_t olap_;
	FastaChunkwiseParser fa_;
	size_t tot_fasta_len_;
	const InputModelUnpaired& model_u_;
	const InputModelUnpaired& model_b_;
	const InputModelPaired& model_c_;
	const InputModelPaired& model_d_;
	FILE *fh_u_;
	FILE *fh_b_1_;
	FILE *fh_b_2_;
	FILE *fh_c_1_;
	FILE *fh_c_2_;
	FILE *fh_d_1_;
	FILE *fh_d_2_;
};

#endif /* defined(__qsim__simplesim__) */
