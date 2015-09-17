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
#include "fasta.h"
#include "input_model.h"
#include "ranlib.hpp"

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
		size_t olap, // amt of overlap between consecutive FASTA windows
		const InputModelUnpaired& model_u,
		const InputModelUnpaired& model_b,
		const InputModelPaired& model_c,
		const InputModelPaired& model_d) :
		fa_(fns, chunksz, olap),
		model_u_(model_u),
		model_b_(model_b),
		model_c_(model_c),
		model_d_(model_d)
	{
		tot_fasta_len_ = estimate_fasta_length(fns);
	}
	
	/**
	 * Simulate a batch of reads
	 */
	void simulate_batch(
		float fraction,
		size_t min_u,
		size_t min_c,
		size_t min_d,
		size_t min_b);

protected:
	
	/**
	 * Given
	 */
	static size_t estimate_fasta_length(const std::vector<std::string>& fns) {
		return 0; //TODO
	}
	
	FastaChunkwiseParser fa_;
	size_t tot_fasta_len_;
	const InputModelUnpaired& model_u_;
	const InputModelUnpaired& model_b_;
	const InputModelPaired& model_c_;
	const InputModelPaired& model_d_;
};

#endif /* defined(__qsim__simplesim__) */
