//
//  input_model.h
//  qsim
//
//  Created by Ben Langmead on 9/15/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#ifndef __qsim__input_model__
#define __qsim__input_model__

#include <stdio.h>
#include <algorithm>
#include "ds.h"
#include "template.h"

/**
 * Encapsulates the input model so we can simulate reads similar to the input
 * reads/alignments.
 */
class InputModelUnpaired {
	
public:
	
	InputModelUnpaired(const EList<TemplateUnpaired>& ts) : ts_(ts) {
		fraglen_avg_ = 0.0f;
		fraglen_max_ = 0;
		for(size_t i = 0; i < ts.size(); i++) {
			fraglen_avg_ += ((float)ts[i].len_ / ts.size());
			fraglen_max_ = std::max(fraglen_max_, (size_t)ts[i].len_);
		}
	}
	
	/**
	 *
	 */
	const TemplateUnpaired& draw() const {
		assert(!empty());
		return ts_[0];
	}
	
	/**
	 * Return true iff no templates were added.
	 */
	bool empty() const {
		return ts_.empty();
	}
	
	/**
	 * Return the number of unpaired input models encountered.
	 * TODO: if we subsample prior to initializing ts_, need to update this
	 * function.
	 */
	size_t num_added() const {
		return ts_.size();
	}

	/**
	 * Return average length of all reads.
	 */
	float avg_len() const {
		return fraglen_avg_;
	}

	/**
	 * Return maximum length of any unpaired template.
	 */
	float max_len() const {
		return fraglen_max_;
	}

protected:
	
	const EList<TemplateUnpaired>& ts_;
	float fraglen_avg_;
	size_t n_;
	size_t fraglen_max_;
};

class InputModelPaired {
	
public:
	
	InputModelPaired(const EList<TemplatePaired>& ts) : ts_(ts) {
		fraglen_avg_ = 0.0f;
		for(size_t i = 0; i < ts.size(); i++) {
			fraglen_avg_ += ((float)ts[i].fraglen_ / ts.size());
			fraglen_max_ = std::max(fraglen_max_, ts[i].fraglen_);
		}
	}
	
	/**
	 *
	 */
	const TemplatePaired& draw() const {
		assert(!empty());
		return ts_[0];
	}

	/**
	 * Return true iff no templates were added.
	 */
	bool empty() const {
		return ts_.empty();
	}

	/**
	 * Return the number of unpaired input models encountered.
	 * TODO: if we subsample prior to initializing ts_, need to update this
	 * function.
	 */
	size_t num_added() const {
		return ts_.size();
	}
	
	/**
	 * Return average length of all fragments.
	 */
	float avg_len() const {
		return fraglen_avg_;
	}

	/**
	 * Return maximum length of any unpaired template.
	 */
	float max_len() const {
		return fraglen_max_;
	}

protected:
	
	const EList<TemplatePaired>& ts_;
	float fraglen_avg_;
	size_t n_;
	size_t fraglen_max_;
};

#endif /* defined(__qsim__input_model__) */
