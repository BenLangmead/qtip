#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define BUFSZ 65536

struct Alignment {
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
    int best_score;
};

/**
 * strtok already used to parse up to rname.
 */
static char * parse_from_rname_on(Alignment *al) {
    al->rname = strtok(NULL, "\t"); assert(al->rname != NULL);
    char *pos_str = strtok(NULL, "\t"); assert(pos_str != NULL);
    al->pos = atoi(pos_str);
    char *mapq_str = strtok(NULL, "\t"); assert(mapq_str != NULL);
    al->mapq = atoi(mapq_str);
    al->cigar = strtok(NULL, "\t"); assert(al->cigar != NULL);
    al->rnext = strtok(NULL, "\t"); assert(al->rnext != NULL);
    char *pnext_str = strtok(NULL, "\t"); assert(pnext_str != NULL);
    al->pnext = atoi(pnext_str);
    char *tlen_str = strtok(NULL, "\t"); assert(tlen_str != NULL);
    al->tlen = atoi(tlen_str);
    al->seq = strtok(NULL, "\t"); assert(al->seq != NULL);
    al->len = strlen(al->seq);
    al->qual = strtok(NULL, "\t"); assert(al->qual != NULL);
    char *extra = strtok(NULL, "\t");
    char *should_be_null = strtok(NULL, "\t");
    assert(should_be_null == NULL);
    return extra;
}

/**
 *
 */
static int parse_ztz_up_to_best_score(char *extra) {
    char *ztz = NULL;
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

#define IS_ALIGNED(al) ((al->flag & 4) == 0)

/**
 * Read the input SAM file while simultaneously writing out records used to
 * train a MAPQ model as well as records used to build an input model.
 */
static int sam_test_pass1(char const *input_sam, const char *output_recs, const char *output_model) {

    char buf_input_sam[BUFSZ];
    FILE *fh = fopen(input_sam, "rb");
    if(fh == NULL) {
        fprintf(stderr, "Could not open SAM file '%s'\n", input_sam);
        return -1;
    }
    setvbuf(fh, buf_input_sam, _IOFBF, BUFSZ);

    char buf_output_recs[BUFSZ];
    FILE *fh_output_recs = fopen(output_recs, "wb");
    if(fh_output_recs == NULL) {
        fprintf(stderr, "Could not open output record file '%s'\n", output_recs);
        return -1;
    }
    setvbuf(fh_output_recs, buf_output_recs, _IOFBF, BUFSZ);

    char buf_output_model[BUFSZ];
    FILE *fh_output_model = fopen(output_model, "wb");
    if(fh_output_model == NULL) {
        fprintf(stderr, "Could not open output model file '%s'\n", output_model);
        return -1;
    }
    setvbuf(fh_output_model, buf_output_model, _IOFBF, BUFSZ);

    /* Advise the kernel of our access pattern.  */
    /* posix_fadvise(fd, 0, 0, 1); */ /* FDADVICE_SEQUENTIAL */

    char linebuf1[BUFSZ], linebuf2[BUFSZ];
    int line1 = 1;

    Alignment al1, al2;
    int al_cur1 = 1;

    int nal = 0, nignored = 0. npair = 0, nunp = 0;

    while(fgets(line, BUFSZ, fh) != NULL) {
        if(line[0] == '@') {
            continue; // skip header
        }
        nal++;
        char *qname = strtok(line, "\t"); assert(qname != NULL);
        char *flag_str = strtok(NULL, "\t"); assert(flag_str != NULL);
        int flag = atoi(flag_str);
        if((flag & 2048) != 0) {
            nignored++;
            continue;
        }

        /* switch which buffer "line" points to */
        char *line = line1 ? linebuf1 : linebuf2;
        line1 = !line1;

        /* somehow switch between alignments? */
        Alignment *al_cur  = al_cur1 ? &al1 : &al2;
        Alignment *al_prev = al_cur1 ? &al2 : &al1;
        al_cur1 = !al_cur1;

        memset(al_cur, 0, sizeof(Alignment));

        al_cur->qname = qname;
        al_cur->flag = flag;
        al_cur->m1_flag = ((flag & 64) != 0) ? '1' : '0';
        if(al_cur->m1_flag != '1') {
            al_cur->m1_flag = ((flag & 128) != 0) ? '2' : '0';
        }

        /* If we're able to mate up ends at this time, do it */
        Alignment *mate1 = NULL, *mate2 = NULL;
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
        }

        int aligned = IS_ALIGNED(al_cur);
        char fw_flag = ((flag & 16) == 0) ? 'T' : 'F';

        /* Case 1: Current read is unpaired and unaligned; nothing more to do! */
        if(al_cur->m1_flag == '0' && !aligned) {
            continue;
        }

        /* Case 2: Current read is unpaired and aligned */
        if(al_cur->m1_flag == '0' && aligned) {
            char *ztz_tok = parse_ztz_up_to_best_score(parse_from_rname_on(al_cur));
            al_cur->best_score = atoi(ztz_tok);

            /* Output information relevant to input model */
            fprintf(fh_output_model, "%d\t%d\t%c\t%s\t%c\t%d\t%d\n",
                    al_cur->best_score,
                    al_cur->len,
                    fw_flag,
                    al_cur->qual,
                    al_cur->m1_flag,
                    al_cur->tlen,
                    al_cur->pos);

            /* Output information relevant to MAPQ model */
            fprintf(fh_output_recs, "%d\t%d\t%d",
                    al_cur->len,
                    al_cur->mapq,
                    al_cur->tlen);

            /* ... including all the ZT:Z fields */
            while(ztz_tok != NULL) {
                fprintf(fh_output_recs, "\t%s", ztz_tok);
                ztz_tok = strtok(NULL, ",");
            }
        }

        /* Case 3: Current read is paired and unaligned, opposite mate is also unaligned; nothing more to do! */
        if(mate1 != NULL) {
            assert(mate2 != NULL);
            if(!IS_ALIGNED(mate1) && !IS_ALIGNED(mate2)) {
                continue;
            }
        }

        /* Case 4: Current read is paired and aligned, opposite mate is unaligned */
        if(mate1 != NULL) {
            assert(mate2 != NULL);
            if(IS_ALIGNED(mate1) && !IS_ALIGNED(mate2)) {
                /* output info to bad-end file */
            }
        }

        /* Case 5: Current read is paired and unaligned, opposite mate is aligned */

        /* Case 6: Current read is paired and both mates aligned, concordantly */

        /* Case 7: Current read is paired and both mates aligned, not condordantly */


        if((flag & 4) == 0) {
        }

        /* Set valid = 0 iff al_cur might be mated up in future */
        al_cur->valid = 0;
        if(mate1 != NULL) {
            /* these two were already mated up */
            al_prev->valid = 0;
        } else if(al_cur->m1_flag != '0') {
            al_cur->valid = 1;
        }
    }

    fclose(fh);
    fclose(fh_output_model);
    fclose(fh_output_recs);

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
    if(argc > 3) {
        return sam_test_pass1(argv[1], argv[2], argv[3]);
    } else {
        fprintf(stderr, "Not enough arguments\n");
    }
}
