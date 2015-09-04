#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

/**
 * Implementations of the various passes that Qsim makes over SAM files.
 *
 * NOTE: uses strtok so not multithreaded.
 */

/* 64K buffer for all input and output */
#define BUFSZ 65536

struct Alignment {
    char *rest_of_line;
    int valid;
    char *qname;
    int flag;
    char m1_flag;
    char *rname;
    int pos;
    int mapq;
    char *cigar;
    char *rnext;
    int pnext;
    int tlen;
    char *seq;
    int len;
    char *qual;
    char *mdz;
    int best_score;

    char *rf_aln_buf;
    size_t rf_aln_buf_len;
    char *rd_aln_buf;
    size_t rd_aln_buf_len;
};

static void clear_alignment(struct Alignment *al) {
    al->rest_of_line = NULL;
    al->valid = 0;
    al->qname = NULL;
    al->flag = 0;
    al->m1_flag = 0;
    al->rname = NULL;
    al->pos = 0;
    al->mapq = 0;
    al->cigar = NULL;
    al->rnext = NULL;
    al->pnext = 0;
    al->tlen = 0;
    al->seq = NULL;
    al->len = 0;
    al->qual = NULL;
    al->mdz = NULL;
    al->best_score = 0;
}

static void init_alignment(struct Alignment *al) {
    clear_alignment(al);
    al->rf_aln_buf = malloc(BUFSZ);
    al->rf_aln_buf_len = BUFSZ;
    al->rd_aln_buf = malloc(BUFSZ);
    al->rd_aln_buf_len = BUFSZ;
}

static void destroy_alignment(struct Alignment *al) {
    if(al->rf_aln_buf != NULL) {
        free(al->rf_aln_buf);
    }
    if(al->rd_aln_buf != NULL) {
        free(al->rd_aln_buf);
    }
}

#define IS_ALIGNED(al) ((al->flag & 4) == 0)

/**
 * strtok already used to parse up to rname.  Parse the rest and return char *
 * to the extra flags.
 *
 * This function uses strtok so it will interrrupt any strtoks already in
 * progress.  But when it returns, it's finished with its stateful strtok
 * session, so there's nothing keeping the caller from starting another
 * strtok.
 */
static char * parse_from_rname_on(struct Alignment *al) {
    assert(al->rest_of_line != NULL);
    al->rname = strtok(al->rest_of_line, "\t"); assert(al->rname != NULL);
    char *pos_str = strtok(NULL, "\t"); assert(pos_str != NULL);
    al->pos = atoi(pos_str);
    char *mapq_str = strtok(NULL, "\t"); assert(mapq_str != NULL);
    al->mapq = atoi(mapq_str);
    assert(al->mapq < 256);
    al->cigar = strtok(NULL, "\t"); assert(al->cigar != NULL);
    al->rnext = strtok(NULL, "\t"); assert(al->rnext != NULL);
    char *pnext_str = strtok(NULL, "\t"); assert(pnext_str != NULL);
    al->pnext = atoi(pnext_str);
    char *tlen_str = strtok(NULL, "\t"); assert(tlen_str != NULL);
    al->tlen = atoi(tlen_str);
    al->seq = strtok(NULL, "\t"); assert(al->seq != NULL);
    al->len = strlen(al->seq);
    al->qual = strtok(NULL, "\t"); assert(al->qual != NULL);
    al->rest_of_line = al->qual + strlen(al->qual) + 1;
    return al->rest_of_line;
}

/**
 * No guarantee about state of strtok upon return.
 */
static char * parse_extra_up_to_ztz(char *extra) {
    char *ztz = NULL;
    extra = strtok(extra, "\t");
    /* we're skipping over everything but ZT:Z */
    while(extra != NULL) {
        if(strncmp(extra, "ZT:Z:", 5) == 0) {
            ztz = extra + 5;
            break;
        }
        extra = strtok(NULL, "\t");
    }
    assert(ztz != NULL);
    char *ztz_tok = strtok(ztz, ",");
    assert(ztz_tok != NULL);
    return ztz_tok;
}

/**
 * No guarantee about state of strtok upon return.
 */
static void print_unpaired(struct Alignment *al, int ordlen, FILE *fh_model, FILE *fh_recs) {
    assert(IS_ALIGNED(al));
    char *extra = parse_from_rname_on(al);
    char *ztz_tok = parse_extra_up_to_ztz(extra);
    al->best_score = atoi(ztz_tok);
    char fw_flag = ((al->flag & 16) == 0) ? 'T' : 'F';

    /* TODO: add stacked alignment info */
    /* TODO: add correct/incorrect info */

    if(fh_model != NULL) {
        /* Output information relevant to input model */
        fprintf(fh_model, "%d\t%c\t%s\t%d\t%c\t%d\n",
                al->best_score,
                fw_flag,
                al->qual,
                al->len,
                al->m1_flag,
                ordlen);
    }

    if(fh_recs != NULL) {
        /* Output information relevant to MAPQ model */
        fprintf(fh_recs, "%d\t%d\t%d",
                al->len,
                al->mapq,
                al->tlen);

        /* ... including all the ZT:Z fields */
        while(ztz_tok != NULL) {
            fprintf(fh_recs, "\t%s", ztz_tok);
            ztz_tok = strtok(NULL, ",");
        }
    }
}

/**
 * No guarantee about state of strtok upon return.
 */
static void print_paired(struct Alignment *al1, struct Alignment *al2, FILE *fh_model, FILE *fh_recs) {
    assert(IS_ALIGNED(al1));
    assert(IS_ALIGNED(al2));

    char *extra1 = parse_from_rname_on(al1);
    char *extra2 = parse_from_rname_on(al2);

    char *ztz_tok1 = parse_extra_up_to_ztz(extra1);
    al1->best_score = atoi(ztz_tok1);
    char fw_flag1 = ((al1->flag & 16) == 0) ? 'T' : 'F';

    int fraglen = 0;
    int upstream1 = 0;

    /* TODO: add stacked alignment info */
    /* TODO: add correct/incorrect info */

    if(fh_recs != NULL) {

        /*
         * Mate 1
         */

        /* Output information relevant to input model */
        fprintf(fh_recs, "%d\t%d\t%d\t%d\t%d",
                al1->len,
                al1->mapq,
                al2->len,
                al2->mapq,
                fraglen);
        /* ... including all the ZT:Z fields */
        while(ztz_tok1 != NULL) {
            int toklen = strlen(ztz_tok1);
            /* remove trailing whitespace */
            while(ztz_tok1[toklen-1] == '\n' || ztz_tok1[toklen-1] == '\r') {
                ztz_tok1[toklen-1] = '\0';
                toklen--;
            }
            fprintf(fh_recs, "\t%s", ztz_tok1);
            ztz_tok1 = strtok(NULL, ",");
        }
    }

    char *ztz_tok2 = parse_extra_up_to_ztz(extra2);
    al2->best_score = atoi(ztz_tok2);
    char fw_flag2 = ((al2->flag & 16) == 0) ? 'T' : 'F';

    if(fh_recs != NULL) {

        /*
         * Mate 2
         */

        /* Output information relevant to input model */
        fprintf(fh_recs, "\t%d\t%d\t%d\t%d\t%d",
                al2->len,
                al2->mapq,
                al1->len,
                al1->mapq,
                fraglen);
        /* ... including all the ZT:Z fields */
        while(ztz_tok2 != NULL) {
            fprintf(fh_recs, "\t%s", ztz_tok2);
            ztz_tok2 = strtok(NULL, ",");
        }
    }

    if(fh_model != NULL) {
        /* Output information relevant to input model */
        fprintf(fh_model, "%d\t%c\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\n",
                al1->best_score + al2->best_score,
                fw_flag1,
                al1->qual,
                al1->best_score,
                al1->len,
                fw_flag2,
                al2->qual,
                al2->best_score,
                al2->len,
                fraglen,
                upstream1);
    }
}

/**
 *
 */
static void parse_cigar_mdz(const char *cigar, const char *mdz) {
    while(1) {
        if(cigar[0] != '\0') {
            return;
        }
        /* step 1: parse op and run from CIGAR string */
        int op = 0, run = 0;
        while(1) {
            int c = *cigar;
            if(isdigit(c)) {
                run *= 10;
                run += c;
            } else if(isalpha(c)) {
                op = c;
                break;
            }
            cigar++;
        }
    }
}

/*
    """ Takes parsed cigar and parsed MD:Z and generates a stacked alignment:
        a pair of strings with gap characters inserted (possibly) and where
        characters at at the same offsets are opposite each other in the
        alignment.  Returns tuple of 2 parallel strings: read string, ref
        string.

        TODO: Fix so it works in local mode. """
    mdp = mdz.mdz_list[:]
    rds, rfs = "", ""
    mdo, rdoff = 0, 0
    for op, run in cigar.cigar_list:
        skipping = (op == 4)
        assert skipping or mdo < len(mdp), "%d, %d" % (mdo, len(mdp))
        if op == 0:   # M
            # Look for block matches and mismatches in MD:Z string
            mdrun = 0
            runleft = run
            while runleft > 0 and mdo < len(mdp):
                op_m, run_m, st_m = mdp[mdo]
                run_comb = min(runleft, run_m)
                runleft -= run_comb
                assert op_m == 0 or op_m == 1
                rds += seq[rdoff:rdoff + run_comb]
                if op_m == 0:
                    # match?
                    rfs += seq[rdoff:rdoff + run_comb]
                else:
                    #mismatch?
                    assert len(st_m) == run_comb
                    rfs += st_m
                    for i in xrange(0, run_comb):
                        assert st_m[i] != seq[rdoff+i:rdoff+i+1] or st_m[i] == 'N'
                mdrun += run_comb
                rdoff += run_comb
                # A stretch of matches in MD:Z could span M and I sections of
                # CIGAR
                if run_comb < run_m:
                    assert op_m == 0
                    mdp[mdo][1] -= run_comb
                else:
                    mdo += 1
        elif op == 1:  # I
            rds += seq[rdoff:rdoff + run]
            rfs += "-" * run
            rdoff += run
        elif op == 2:  # D
            op_m, run_m, st_m = mdp[mdo]
            assert op_m == 2
            assert run == run_m
            assert len(st_m) == run
            mdo += 1
            rds += "-" * run
            rfs += st_m
        elif op == 3:  # N
            # If this is really long, this probably isn't the best way to do
            # this
            rds += "-" * run
            rfs += "-" * run
        elif op == 4:  # S
            rdoff += run
        elif op == 5:  # H
            pass
        elif op == 6:  # P
            assert False
        elif op == 7:  # =
            assert False
        elif op == 8:  # X
            assert False
        else:
            assert False
    assert mdo == len(mdp)
    return rds, rfs
*/

static void parse_wgsim_correct() {
/*
    _wgsimex_re = re.compile('(.+)_([^_]+)_([^_]+)_([^:]+):([^:]+):([^_]+)_([^:]+):([^:]+):([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^/]+).*')
    #                           1   2       3       4       5       6       7       8       9       10      11      12      13

    def isExtendedWgsim(nm):
        ret = _wgsimex_re.match(nm) is not None
        return ret

    def parseExtendedWgsim(al):
        ''' Note: this is my extended version of wgsim's output '''
        nm = al.name
        res = _wgsimex_re.match(nm)
        refid, fragst1, fragen1 = res.group(1), int(res.group(2))-1, int(res.group(3))-1
        len1, len2 = int(res.group(10)), int(res.group(11))
        flip = res.group(12) == '1'
        mate1 = al.mate1 or not al.paired
        ln = len1 if mate1 else len2
        if (not flip) == mate1:
            return fragst1, refid, True
        else:
            return fragen1 - (ln-1), refid, False
*/
}

/**
 * Read the input SAM file while simultaneously writing out records used to
 * train a MAPQ model as well as records used to build an input model.
 */
static int sam_test_pass1(char const *input_sam,
                          const char *orec_u_fn, const char *omod_u_fn,
                          const char *orec_b_fn, const char *omod_b_fn,
                          const char *orec_c_fn, const char *omod_c_fn,
                          const char *orec_d_fn, const char *omod_d_fn,
                          int quiet)
{
    char buf_input_sam[BUFSZ];
    FILE *fh = fopen(input_sam, "rb");
    if(fh == NULL) {
        fprintf(stderr, "Could not open SAM file '%s'\n", input_sam);
        return -1;
    }
    setvbuf(fh, buf_input_sam, _IOFBF, BUFSZ);

    /* unpaired */

    char orec_u_buf[BUFSZ];
    FILE *orec_u_fh = fopen(orec_u_fn, "wb");
    if(orec_u_fh == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", orec_u_fn);
        return -1;
    }
    setvbuf(orec_u_fh, orec_u_buf, _IOFBF, BUFSZ);

    char omod_u_buf[BUFSZ];
    FILE *omod_u_fh = fopen(omod_u_fn, "wb");
    if(omod_u_fh == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", omod_u_fn);
        return -1;
    }
    setvbuf(omod_u_fh, omod_u_buf, _IOFBF, BUFSZ);

    /* bad-end */

    char orec_b_buf[BUFSZ];
    FILE *orec_b_fh = fopen(orec_b_fn, "wb");
    if(orec_b_fh == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", orec_b_fn);
        return -1;
    }
    setvbuf(orec_b_fh, orec_b_buf, _IOFBF, BUFSZ);

    char omod_b_buf[BUFSZ];
    FILE *omod_b_fh = fopen(omod_b_fn, "wb");
    if(omod_b_fh == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", omod_b_fn);
        return -1;
    }
    setvbuf(omod_b_fh, omod_b_buf, _IOFBF, BUFSZ);

    /* concordant */

    char orec_c_buf[BUFSZ];
    FILE *orec_c_fh = fopen(orec_c_fn, "wb");
    if(orec_c_fh == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", orec_c_fn);
        return -1;
    }
    setvbuf(orec_c_fh, orec_c_buf, _IOFBF, BUFSZ);

    char omod_c_buf[BUFSZ];
    FILE *omod_c_fh = fopen(omod_c_fn, "wb");
    if(omod_c_fh == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", omod_c_fn);
        return -1;
    }
    setvbuf(omod_c_fh, omod_c_buf, _IOFBF, BUFSZ);

    /* discordant */

    char orec_d_buf[BUFSZ];
    FILE *orec_d_fh = fopen(orec_d_fn, "wb");
    if(orec_d_fh == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", orec_d_fn);
        return -1;
    }
    setvbuf(orec_d_fh, orec_d_buf, _IOFBF, BUFSZ);

    char omod_d_buf[BUFSZ];
    FILE *omod_d_fh = fopen(omod_d_fn, "wb");
    if(omod_d_fh == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", omod_d_fn);
        return -1;
    }
    setvbuf(omod_d_fh, omod_d_buf, _IOFBF, BUFSZ);

    /* Advise the kernel of our access pattern.  */
    /* posix_fadvise(fd, 0, 0, 1); */ /* FDADVICE_SEQUENTIAL */

    char linebuf1[BUFSZ], linebuf2[BUFSZ];
    int line1 = 1;

    struct Alignment al1, al2;
    init_alignment(&al1);
    init_alignment(&al2);

    int al_cur1 = 1;

    int nline = 0, nignored = 0, npair = 0, nunp = 0;
    int nunp_al = 0, nunp_unal = 0, npair_badend = 0, npair_conc = 0, npair_disc = 0, npair_unal = 0;

    while(1) {
        char *line = line1 ? linebuf1 : linebuf2;
        if(fgets(line, BUFSZ, fh) == NULL) {
            break; /* done */
        }
        if(line[0] == '@') {
            continue; // skip header
        }
        nline++;
        char *qname = strtok(line, "\t"); assert(qname != NULL);
        assert(qname == line);
        char *flag_str = strtok(NULL, "\t"); assert(flag_str != NULL);
        int flag = atoi(flag_str);
        if((flag & 2048) != 0) {
            nignored++;
            continue;
        }

        /* switch which buffer "line" points to */
        line1 = !line1;

        /* somehow switch between alignments? */
        struct Alignment *al_cur  = al_cur1 ? &al1 : &al2;
        assert(!al_cur->valid);
        clear_alignment(al_cur);
        struct Alignment *al_prev = al_cur1 ? &al2 : &al1;
        al_cur1 = !al_cur1;

        al_cur->rest_of_line = flag_str + strlen(flag_str) + 1; /* for re-parsing */
        al_cur->qname = qname;
        al_cur->flag = flag;
        al_cur->m1_flag = ((flag & 64) != 0) ? '1' : (((flag & 128) != 0) ? '2' : '0');

        /* If we're able to mate up ends at this time, do it */
        struct Alignment *mate1 = NULL, *mate2 = NULL;
        if(al_cur->m1_flag != '0' && al_prev->valid) {
            if(al_cur->m1_flag == '1') {
                assert(al_prev->m1_flag == '2');
                mate1 = al_cur;
                mate2 = al_prev;
            } else {
                assert(al_cur->m1_flag == '2');
                assert(al_prev->m1_flag == '1');
                mate1 = al_prev;
                mate2 = al_cur;
            }
            mate1->valid = mate2->valid = 0;
            npair++;
        }

        if(al_cur->m1_flag == '0') {
            nunp++;

            /* Case 1: Current read is unpaired and unlineigned, we can safely skip */
            if(!IS_ALIGNED(al_cur)) {
                nunp_unal++;
                continue;
            }

            /* Case 2: Current read is unpaired and aligned */
            else {
                nunp_al++;
                print_unpaired(al_cur, 0, omod_u_fh, orec_u_fh);
            }
        }

        else if(mate1 != NULL) {
            /* Case 3: Current read is paired and unlineigned, opposite mate is also unlineigned; nothing more to do! */
            assert(mate2 != NULL);
            if(!IS_ALIGNED(mate1) && !IS_ALIGNED(mate2)) {
                npair_unal++;
                continue;
            }

            /* Case 4: Current read is paired and aligned, opposite mate is unlineigned */
            /* Case 5: Current read is paired and unlineigned, opposite mate is aligned */
            /* we handle both here */
            else if(IS_ALIGNED(mate1) != IS_ALIGNED(mate2)) {
                npair_badend++;
                print_unpaired(
                    IS_ALIGNED(mate1) ? mate1 : mate2,
                    IS_ALIGNED(mate1) ? mate2->len : mate1->len,
                    omod_b_fh, orec_b_fh);
            }

            else {
                assert((mate1->flag & 2) == (mate1->flag & 2));
                int concordant = ((mate1->flag & 2) != 0);

                if(concordant) {
                    /* Case 6: Current read is paired and both mates aligned, concordantly */
                    npair_conc++;
                    print_paired(mate1, mate2, omod_c_fh, orec_c_fh);
                }

                else {
                    /* Case 7: Current read is paired and both mates aligned, not condordantly */
                    npair_disc++;
                    print_paired(mate1, mate2, omod_d_fh, orec_d_fh);
                }
            }
        }

        else {
            /* This read is paired but we haven't seen the mate yet */
            assert(al_cur->m1_flag != '0');
            al_cur->valid = 1;
        }
    }

    if(!quiet) {
        fprintf(stderr, "%d lines\n", nline);
        fprintf(stderr, "%d ignored b/c secondary\n", nignored);
        fprintf(stderr, "%d unpaired\n", nunp);
        fprintf(stderr, "    %d aligned\n", nunp_al);
        fprintf(stderr, "    %d unaligned\n", nunp_unal);
        fprintf(stderr, "%d paired-end\n", npair);
        fprintf(stderr, "    %d concordant\n", npair_conc);
        fprintf(stderr, "    %d discordant\n", npair_disc);
        fprintf(stderr, "    %d bad-end\n", npair_badend);
        fprintf(stderr, "    %d unaligned\n", npair_unal);
    }

    fclose(fh);
    fclose(omod_u_fh);
    fclose(orec_u_fh);
    fclose(omod_b_fh);
    fclose(orec_b_fh);
    fclose(omod_c_fh);
    fclose(orec_c_fh);
    fclose(omod_d_fh);
    fclose(orec_d_fh);

    destroy_alignment(&al1);
    destroy_alignment(&al2);

    return 0;
}

/**
 * Read the input SAM file in tandem with the MAPQ prediction file and replace
 * old MAPQs with new ones.
 */
int sam_test_pass2(char const *input_sam, const char *predictions) {
    return 0;
}

/**
 * Read the training SAM file (alignments for tandem reads) and
 */
int sam_training_pass(char const *input_sam, const char *predictions) {
    return 0;
}

int main(int argc, char **argv) {
    if(argc == 2) {
        return sam_test_pass1(
            argv[1],
            "orec_u.tsv", "omod_u.tsv",
            "orec_b.tsv", "omod_b.tsv",
            "orec_c.tsv", "omod_c.tsv",
            "orec_d.tsv", "omod_d.tsv",
            0);
    } else if(argc == 10) {
        return sam_test_pass1(
            argv[1],
            argv[2], argv[3],
            argv[4], argv[5],
            argv[6], argv[7],
            argv[8], argv[9],
            0);
    } else {
        fprintf(stderr, "Not enough arguments\n");
    }
}
