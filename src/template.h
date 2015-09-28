//
//  template.h
//  qsim
//
//  Created by Ben Langmead on 9/13/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#ifndef qsim_template_h
#define qsim_template_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "edit_xscript.h"

/*
 * Encapsulates:
 *
 * 1. Best score
 * 2. FW flag (T or F)
 * 3. Quality string
 * 4. Read length
 * 5. Mate flag (0, 1 or 2)
 * 6. Opposite mate read length
 * 7. Edit transcript
 */
struct TemplateUnpaired {
	
	TemplateUnpaired() :
		best_score_(0),
		fw_flag_(0),
		len_(0),
		mate_flag_(0),
		opp_len_(0),
		qual_(NULL),
		edit_xscript_(NULL) { }
	
	TemplateUnpaired(
		int best_score,
		int len,
		char fw_flag,
		char mate_flag,
		int opp_len,
		const char *qual,
		const char *edit_xscript)
	{
		init(best_score, len, fw_flag, mate_flag, opp_len, qual, edit_xscript);
	}
	
	void init(
		int best_score,
		int len,
		char fw_flag,
		char mate_flag,
		int opp_len,
		const char *qual,
		const char *edit_xscript)
	{
		best_score_ = best_score;
		len_ = len;
		fw_flag_ = fw_flag;
		mate_flag_ = mate_flag;
		opp_len_ = opp_len;
		assert(qual != NULL);
		if(qual != NULL) {
			qual_ = strdup(qual);
			assert(qual_ != NULL);
		}
		assert(edit_xscript != NULL);
		if(edit_xscript != NULL) {
			edit_xscript_ = strdup(edit_xscript);
			assert(edit_xscript_ != NULL);
		}
	}
	
	~TemplateUnpaired() {
		if(qual_ != NULL) {
			free(qual_);
			qual_ = NULL;
		}
		if(edit_xscript_ != NULL) {
			free(edit_xscript_);
			edit_xscript_ = NULL;
		}
	}
	
	/**
	 * Return number of reference characters involved in alignment, as
	 * indicated by edit transcript.
	 */
	size_t reflen() const {
		assert(edit_xscript_ != NULL);
		return edit_xscript_to_rflen(edit_xscript_);
	}
	
	int best_score_;
	char fw_flag_;
	int len_;
	char mate_flag_;
	int opp_len_;
	
	char *qual_;
	char *edit_xscript_;
};

/*
 * Paired-end columns:
 * 1. Sum of best scores of both mates
 * 2. Mate 1 FW flag (T or F)
 * 3. Mate 1 quality string
 * 4. Mate 1 best score
 * 5. Mate 1 read length
 * 6. Mate 1 edit transcript
 * 7 .Mate 2 FW flag (T or F)
 * 8. Mate 2 quality string
 * 9. Mate 2 best score
 * 10. Mate 2 read length
 * 11. Mate 2 edit transcript
 * 12. Mate-1-upstream flag (T or F)
 * 13. Fragment length
 */
struct TemplatePaired {
	
	TemplatePaired() :
		score_12_(0),
		score_1_(0),
		len_1_(0),
		fw_flag_1_(0),
		qual_1_(NULL),
		edit_xscript_1_(NULL),
		score_2_(0),
		len_2_(0),
		fw_flag_2_(0),
		qual_2_(NULL),
		edit_xscript_2_(NULL),
		upstream1_(false),
		fraglen_(0) { }
	
	TemplatePaired(
		int score_12,
		int score_1,
		int len_1,
		char fw_flag_1,
		const char *qual_1,
		const char *edit_xscript_1,
		int score_2,
		int len_2,
		char fw_flag_2,
		const char *qual_2,
		const char *edit_xscript_2,
		bool upstream1,
		size_t fraglen)
	{
		init(score_12,
			 score_1, len_1, fw_flag_1, qual_1, edit_xscript_1,
			 score_2, len_2, fw_flag_2, qual_2, edit_xscript_2,
			 upstream1, fraglen);
	}
	
	void init(
		int score_12,
		int score_1,
		int len_1,
		char fw_flag_1,
		const char *qual_1,
		const char *edit_xscript_1,
		int score_2,
		int len_2,
		char fw_flag_2,
		const char *qual_2,
		const char *edit_xscript_2,
		bool upstream1,
		size_t fraglen)
	{
		score_12_ = score_12;
		score_1_ = score_1;
		len_1_ = len_1;
		fw_flag_1_ = fw_flag_1;
		score_2_ = score_2;
		len_2_ = len_2;
		fw_flag_2_ = fw_flag_2;
		upstream1_  = upstream1;
		fraglen_  = fraglen;
		if(qual_1 != NULL) {
			qual_1_ = strdup(qual_1);
		}
		if(qual_2 != NULL) {
			qual_2_ = strdup(qual_2);
		}
		if(edit_xscript_1 != NULL) {
			edit_xscript_1_ = strdup(edit_xscript_1);
		}
		if(edit_xscript_2 != NULL) {
			edit_xscript_2_ = strdup(edit_xscript_2);
		}
	}
	
	~TemplatePaired() {
		if(qual_1_ != NULL) {
			free(qual_1_);
			qual_1_ = NULL;
		}
		if(qual_2_ != NULL) {
			free(qual_2_);
			qual_2_ = NULL;
		}
		if(edit_xscript_1_ != NULL) {
			free(edit_xscript_1_);
			edit_xscript_1_ = NULL;
		}
		if(edit_xscript_2_ != NULL) {
			free(edit_xscript_2_);
			edit_xscript_2_ = NULL;
		}
	}
	
	int score_12_;
	int score_1_;
	int len_1_;
	char fw_flag_1_;
	char *qual_1_;
	char *edit_xscript_1_;
	int score_2_;
	int len_2_;
	char fw_flag_2_;
	char *qual_2_;
	char *edit_xscript_2_;
	bool upstream1_;
	size_t fraglen_;
};

#endif
