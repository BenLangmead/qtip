//
//  fasta.h
//  qsim
//
//  Created by Ben Langmead on 9/4/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#ifndef __qsim__fasta__
#define __qsim__fasta__

#include <stdio.h>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cassert>
#include <limits>

/**
 * Class that facilitates iterating through overlapping stretches of all the
 * entries in one or more (multi-)FASTA files.
 */
class FastaChunkwiseParser {
public:
	FastaChunkwiseParser(
		const std::vector<std::string>& fns,
		size_t chunksz,
		size_t olap) :
		fns_(fns),
		fni_(0),
		foff_(0),
		refoff_(0),
		fh_(NULL),
		buf_(NULL),
		bufcur_(0),
		chunksz_(chunksz),
		olap_(olap),
		pushback_(std::numeric_limits<int>::min())
	{
		assert(chunksz > olap);
		buf_ = new char[chunksz];
	}
	
	~FastaChunkwiseParser() {
		delete[] buf_;
	}
	
	/**
	 * Reset back to first chunk of first FASTA.
	 */
	void reset() {
		fni_ = foff_ = bufcur_ = refoff_ = 0;
		pushback_ = std::numeric_limits<int>::min();
		if(fh_ != NULL) {
			fclose(fh_); fh_ = NULL;
		}
	}
	
	/**
	 * Return true if we've iterated through all chunks of all FASTAs.
	 */
	inline bool done() const {
		return fni_ >= fns_.size();
	}

		/**
	 * Returns pointer to sequence.
	 */
	const char *next(
		std::string& refid, // out: name of reference buffer is from
		std::string& refid_full, // out: name of reference buffer is from
		size_t& refoff, // out: reference offset of first character in buffer
		size_t& retsz); // out: number of characters in returned buffer
	
protected:

	inline void pushback(int c) {
		assert(pushback_ == std::numeric_limits<int>::min());
		pushback_ = c;
	}
	
	inline int get(FILE *fh) {
		if(pushback_ != std::numeric_limits<int>::min()) {
			int ret = pushback_;
			pushback_ = std::numeric_limits<int>::min();
			return ret;
		}
		return fgetc(fh);
	}
	
	const static size_t FASTA_BUFSZ = 65536;

	const std::vector<std::string> fns_;
	size_t fni_; // offset into list of files
	size_t foff_; // offset into current FASTA file
	size_t refoff_; // offset into current ref
	FILE *fh_;
	char *buf_;
	size_t bufcur_; // current number of chars in buf_
	const size_t chunksz_;
	const size_t olap_;
	char fabuf_[FASTA_BUFSZ];
	int pushback_;
};

#endif /* defined(__qsim__fasta__) */
