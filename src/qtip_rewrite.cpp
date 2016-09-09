//
//  qtip_rewrite.cpp
//  qtip
//
//  Created by Ben Langmead on 9/25/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <cassert>
#include "qtip_rewrite.h"

using namespace std;

/**
 * Re-write a SAM file to include our new MAPQ prediction.
 */

bool write_orig_mapq = false;
const char *orig_mapq_flag = "Zm:i";

bool write_precise_mapq = false;
const char *precise_mapq_flag = "Zp:Z";

bool keep_ztz = false;

const static size_t BUFSZ = 65536;

#define FILEDEC(fn, fh, buf, typ, do_open) \
	char buf [BUFSZ]; \
	FILE * fh = NULL; \
	if(do_open) { \
		fh = fopen( fn .c_str(), "wb"); \
		if(fh == NULL) { \
			cerr << "Could not open output " << typ << " file \"" << fn << "\"" << endl; \
			return -1; \
		} \
		setvbuf(fh, buf, _IOFBF, BUFSZ); \
	}

/**
 * Read next prediction from the prediction file and popular line and mapq out
 * paramters appropriately.
 */
static bool next_prediction(FILE *fh, size_t& line, float& mapq) {
	char linebuf[BUFSZ];
	static size_t last_line = 0;
	if(fgets(linebuf, BUFSZ, fh) == NULL) {
		return false; // done
	}
	sscanf(linebuf, "%llu,%f", (unsigned long long*)&line, &mapq);
	assert(last_line == 0 || line > last_line);
	last_line = line;
	assert(mapq >= 0.0f);
	return true;
}

/**
 * Write a new line of SAM (buf) to output filehandle (osam_fh) replacing the
 * existing MAPQ with the predicted one (mapq).
 */
static void rewrite(FILE *fh, char *buf, float mapq) {
	char orig[10];
	char *orig_cur = orig;
	for(int i = 0; i < 4; i++) {
		while(*buf != '\t') fputc(*buf++, fh);
		fputc(*buf++, fh);
	}
	// Replace MAPQ with our new one
	int mapq_rounded = (int)(mapq + 0.5f);
	fprintf(fh, "%d", mapq_rounded);
	while(*buf != '\t') {
		// Copy orginal to buffer?
		*orig_cur++ = *buf++;
	}
	*orig_cur = '\0';
	while(*buf != '\n' && *buf != '\r') {
		if(!keep_ztz && *buf == '\t') {
			if(strncmp(buf+1, "ZT:Z:", 5) == 0) { // Remove the ZT:Z
				buf += 6;
				while(true) {
					if(*buf == '\t' || *buf == '\n' || *buf == '\r') {
						break;
					}
					buf++;
				}
			}
			if(*buf == '\n' || *buf == '\r') {
				break;
			}
		}
		fputc(*buf, fh);
		buf++;
	}
	if(write_orig_mapq) {
		fprintf(fh, "\t%s:%s", orig_mapq_flag, orig);
	}
	if(write_precise_mapq) {
		fprintf(fh, "\t%s:%0.3f", precise_mapq_flag, mapq);
	}
	fputc('\n', fh);
}

int main(int argc, char **argv) {

	if(argc == 1) {
		// print which arguments from ts.py should pass through to here
		cout << "orig-mapq-flag "
		     << "precise-mapq-flag "
		     << "write-orig-mapq "
		     << "write-precise-mapq "
		     << "keep-ztz" << endl;
		return 0;
	}

	string fn;
	string prefix;
	vector<string> sams;
	vector<string> preds;
	string osam_fn;
	
	// All arguments except last are SAM files to parse.  Final argument is
	// prefix for output files.
	{
		int section = 0, prefix_set = 0;
		for(int i = 1; i < argc; i++) {
			if(strcmp(argv[i], "--") == 0) {
				section++;
				continue;
			}
			if(section == 0) {
				if(i == argc-1) {
					cerr << "Error: odd number of arguments in options section" << endl;
					throw 1;
				}
				if(strcmp(argv[i], "orig-mapq-flag") == 0) {
					orig_mapq_flag = argv[++i];
				}
				if(strcmp(argv[i], "precise-mapq-flag") == 0) {
					precise_mapq_flag = argv[++i];
				}
				if(strcmp(argv[i], "write-orig-mapq") == 0) {
					write_orig_mapq = strcmp(argv[++i], "True") == 0;
				}
				if(strcmp(argv[i], "write-precise-mapq") == 0) {
					write_precise_mapq = strcmp(argv[++i], "True") == 0;
				}
				if(strcmp(argv[i], "keep-ztz") == 0) {
					keep_ztz = strcmp(argv[++i], "True") == 0;
				}
			} else if(section == 1) {
				sams.push_back(string(argv[i]));
			} else if(section == 2) {
				preds.push_back(string(argv[i]));
			} else {
				prefix = argv[i];
				prefix_set++;
				osam_fn = prefix + string(".sam");
			}
			if(prefix_set > 1) {
				cerr << "Warning: More than output prefix specified; using last one: \"" << prefix << "\"" << endl;
			}
		}
		if(sams.empty() || !prefix_set) {
			cerr << "Usage: qtip_rewrite" << endl;
		}
	}

	// Output SAM file
	FILEDEC(osam_fn, osam_fh, osam_buf, "SAM", true);

	char buf_input_sam[BUFSZ];
	char buf_input_pred[BUFSZ];

	for(size_t i = 0; i < sams.size(); i++) {
		// Input SAM file
		cerr << "Parsing SAM file \"" << sams[i] << "\"" << endl;
		FILE *fh_sam = fopen(sams[i].c_str(), "rb");
		if(fh_sam == NULL) {
			cerr << "Could not open input SAM file \"" << sams[i] << "\"" << endl;
			return -1;
		}
		setvbuf(fh_sam, buf_input_sam, _IOFBF, BUFSZ);
		
		// Input prediction file
		cerr << "Parsing prediction file \"" << preds[i] << "\"" << endl;
		FILE *fh_pred = fopen(preds[i].c_str(), "rb");
		if(fh_pred == NULL) {
			cerr << "Could not open input prediction file \"" << preds[i] << "\"" << endl;
			return -1;
		}
		setvbuf(fh_pred, buf_input_pred, _IOFBF, BUFSZ);
		
		bool done_with_predictions = false;
		bool done_with_sam = false;
		char linebuf[BUFSZ];
		size_t nline = 0, nhead = 0;
		while(!done_with_predictions || !done_with_sam) {
			size_t pred_line = 0;
			float mapq = 0.0f;
			if(!next_prediction(fh_pred, pred_line, mapq)) {
				done_with_predictions = true;
			}
			while(true) {
				// Handle line of sam
				if(fgets(linebuf, BUFSZ, fh_sam) == NULL) {
					assert(done_with_predictions);
					done_with_sam = true;
					break;
				}
				nline++;
				assert(done_with_predictions || nline <= pred_line);
				if(linebuf[0] == '@') {
					nhead++;
					fputs(linebuf, osam_fh);
					continue; // skip header
				}
				if(done_with_predictions || pred_line > nline) {
					fputs(linebuf, osam_fh); // no prediction for this line
					continue;
				}
				assert(nline == pred_line); // there is a prediction
				rewrite(osam_fh, linebuf, mapq);
				break; // get next prediction
			}
		}
		assert(done_with_predictions && done_with_sam);
		
		fclose(fh_sam);
		fclose(fh_pred);
	}
}
