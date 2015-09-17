//
//  simplesim.cpp
//  qsim
//
//  Created by Ben Langmead on 9/12/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#include "simplesim.h"
#include "fasta.h"
#include <vector>

/**
 * Return draw from binomial distribution with given n, p.
 */
static inline size_t draw_binomial(size_t n, float p) {
	return (size_t)ignbin((int)n, p);
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
	while(true) {
		const char * buf = fa_.next(refid, refid_full, refoff, retsz);
		if(buf == NULL && fa_.done()) {
			break;
		}
		size_t nu_chances = retsz - model_u_.avg_len() + 1;
		size_t nb_chances = retsz - model_b_.avg_len() + 1;
		size_t nc_chances = retsz - model_c_.avg_len() + 1;
		size_t nd_chances = retsz - model_d_.avg_len() + 1;
		size_t nu_samp = draw_binomial(nu, float(nu_chances) / tot_fasta_len_);
		size_t nb_samp = draw_binomial(nb, float(nb_chances) / tot_fasta_len_);
		size_t nc_samp = draw_binomial(nc, float(nc_chances) / tot_fasta_len_);
		size_t nd_samp = draw_binomial(nd, float(nd_chances) / tot_fasta_len_);
	}
}
