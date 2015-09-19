//
//  edit_xscript.h
//  qsim
//
//  Created by Ben Langmead on 9/19/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#ifndef qsim_edit_xscript_h
#define qsim_edit_xscript_h

/**
 * Given edit transcript, get total number of reference characters involved.
 */
static inline size_t edit_xscript_to_rflen(const char *edit_xscript) {
	const char *cur = edit_xscript;
	size_t rflen = 0;
	while(*cur != '\0') {
		if(*cur == 'S' || *cur == '=' || *cur == 'X' || *cur == 'D') {
			rflen++;
		}
		cur++;
	}
	return rflen;
}

#endif
