#ifndef __qsim__align__
#define __qsim__align__

#include <stdio.h>
#include <cassert>
#include "ds.h"

class EditDistanceSolver {
	
public:
	
	EditDistanceSolver(size_t max_seq_size) :
		x_(NULL), xlen_(0), y_(NULL), ylen_(0), sz_(0), filled_(false)
	{
		max_seq_size++;
		mat_.reserveExact(max_seq_size * max_seq_size);
	}
	
	/**
	 * Initialize with new X and Y.
	 */
	void init(const char *x, size_t xlen, const char *y, size_t ylen) {
		sz_ = (xlen + 1) * (ylen + 1);
		mat_.reserveExact(sz_);
		filled_ = false;
	}
	
	/**
	 * Fill matrix.
	 */
	void fill() {
		{
			int dist = 0;
			for(int i = 0; i < sz_; i += (ylen_+1)) {
				mat_[i] = dist++;
			}
		}
		for(int j = 0; j <= ylen_; j++) {
			mat_[j] = j;
		}
		for(size_t i = 1; i <= xlen_; i++) {
			for(size_t j = 1; j <= ylen_; j++) {
				
			}
		}
		filled_ = true;
	}
	
	/**
	 *
	 */
	void edit_transcript() {
		assert(filled_);
	}
	
protected:
	
	const char *x_;
	size_t xlen_;
	const char *y_;
	size_t ylen_;
	size_t sz_;
	EList<int> mat_;
	bool filled_;
};

#endif /* defined(__qsim__align__) */
