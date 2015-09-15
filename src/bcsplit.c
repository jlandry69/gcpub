#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <libgen.h>
#include <time.h>
#include "array.h"
#include "klib/kseq.h"
#include "klib/khash.h"

#define DNA_ALPHA "ACGT"
#define DEFAULT_NUM_MISMATCHES 0
#define DEFAULT_NUM_SPACER_BASES 1 // e.g. for T spacer
#define DEFAULT_OUTPUT_PREFIX "bcsplit"

/*
 * barcode record. holds:
 *   - barcode ID
 *   - barcode sequence
 *   - counter of found tags
 *   - file pointers to output files
 */
typedef struct {
    char *id, *seq;
    unsigned num_found;
    gzFile *fp1, *fp2;
} BcRec;

ARRAY_DECLARE(BcRec);
KSEQ_INIT(gzFile, gzread);
KHASH_MAP_INIT_STR(str, int);

/* {{{ print_main_usage() */
void print_main_usage(char *progname)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s <mode> [arguments]\n\n", progname);
    fprintf(stderr, "Mode:  pe  paired-end\n");
    fprintf(stderr, "       se  single-end\n");
    fprintf(stderr, "\n");
}
/* }}} */

/* {{{ print_pe_usage() */
void print_pe_usage(char *progname)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   %s pe [options] <barcodes.txt> <fq1.gz> <fq2.gz>\n", progname);
    fprintf(stderr, "\n");
    fprintf(stderr, "Options: -m INT   number mismatches 0/1 [%u]\n", DEFAULT_NUM_MISMATCHES);
    fprintf(stderr, "         -s INT   number spacer bases [%u]\n", DEFAULT_NUM_SPACER_BASES);
    fprintf(stderr, "         -o STR   output prefix [\"%s\"]\n", DEFAULT_OUTPUT_PREFIX);
    fprintf(stderr, "         -c       only count\n");
    fprintf(stderr, "\n");
}
/* }}} */

/* {{{ print_se_usage() */
void print_se_usage(char *progname)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   %s se [options] <barcodes.txt> <fq.gz>\n", progname);
    fprintf(stderr, "\n");
    fprintf(stderr, "Options: -m INT   number mismatches 0/1 [%u]\n", DEFAULT_NUM_MISMATCHES);
    fprintf(stderr, "         -s INT   number spacer bases [%u]\n", DEFAULT_NUM_SPACER_BASES);
    fprintf(stderr, "         -o STR   output prefix [\"%s\"]\n", DEFAULT_OUTPUT_PREFIX);
    fprintf(stderr, "         -c       only count\n");
    fprintf(stderr, "\n");
}
/* }}} */

/* {{{ get_num_mismatch */
int get_num_mismatch(char *s, unsigned *nm)
{
    if (sscanf(s, "%u", nm) != 1 || (*nm  != 0 && *nm != 1)) {
        fprintf(stderr, "argument <num-mismatch> has to be 0 or 1\n");
        return -1;
    }
    return 0;
}       
/* }}} */

int split_pe(int argc, char *argv[], char *progname)
{
    int c, bc_len = -1, ret, i, j, bc_idx, only_count = 0;
    unsigned num_mismatches = DEFAULT_NUM_MISMATCHES, 
             num_spacer_bases = DEFAULT_NUM_SPACER_BASES, 
             dna_alpha_len = strlen(DNA_ALPHA),
             num_undetermined = 0;
    char *out_prefix = NULL, *fn, **sptr,
         bc_id[1024], bc_seq[1024], bc_seq_cpy[1024]; // hello, buffer overflow
    clock_t t = clock();
    BcRec bc;
    ArrayBcRec bcs;
    FILE *fp;
    kseq_t *seq1, *seq2;
    khash_t(str) *h = kh_init(str);
    khint_t k, k2;
    gzFile *fp1, *fp2;

    ARRAY_INIT(&bcs, BcRec, 1000);

    while ((c = getopt(argc, argv, "m:s:o:c")) >= 0) {
        switch (c) {
            case 'm': if (sscanf(optarg, "%u", &num_mismatches) != 1) {
                          fprintf(stderr, "Error: option -m expects unsigned int\n");
                          return -1;
                      }
                      break;
            case 's': if (sscanf(optarg, "%u", &num_spacer_bases) != 1) {
                          fprintf(stderr, "Error: option -s expects unsigned int\n");
                          return -1;
                      }
                      break;
            case 'o': out_prefix = strdup(optarg);
                      break;
            case 'c': only_count = 1;
                      break;
        }
    }

    if (optind + 3 != argc) {
        print_pe_usage(progname);
        return -1;
    }

    if (num_mismatches != 0 && num_mismatches != 1) {
        fprintf(stderr, "Error: argument -m has to be 0 or 1\n");
        return -1;
    }

    if (out_prefix == NULL) {
        out_prefix = strdup(DEFAULT_OUTPUT_PREFIX);
    }

    for (sptr = argv+optind; sptr-argv<argc; sptr++) {
        if (access(*sptr, F_OK) == -1) {
            fprintf(stderr, "Error: file %s does not exist\n", *sptr);
            return -1;
        }
    }

    fprintf(stderr, "[barcode file: %s]\n", argv[optind]);
    fprintf(stderr, "[fastq file1: %s]\n", argv[optind+1]);
    fprintf(stderr, "[fastq file2: %s]\n", argv[optind+2]);
    fprintf(stderr, "[number of mismatches allowed: %u]\n", num_mismatches);
    fprintf(stderr, "[number of spacer bases: %u]\n", num_spacer_bases);
    fprintf(stderr, "[output prefix: %s]\n", out_prefix);
    fprintf(stderr, "[only count: %s]\n", only_count ? "true" : "false");

    /* read barcode file */
    if ((fp = fopen(argv[optind], "r")) == NULL) {
        fprintf(stderr, "Error: cannot open barcode file %s\n", argv[optind]);
        return -1;
    }

    while (fscanf(fp, "%s %s", bc_id, bc_seq) == 2) {
        bc_len = strlen(bc_seq);
        bc.id = strdup(bc_id);
        bc.seq = strdup(bc_seq);
        bc.num_found = 0;
        if (!only_count) {
            fn = (char*)calloc(strlen(out_prefix) + 3 + strlen(bc_id) + 6 + 1, sizeof(char));
            strcpy(fn, out_prefix);
            strcat(fn, "_1_");
            strcat(fn, bc_id);
            strcat(fn, ".fq.gz");
            bc.fp1 = gzopen(fn, "w");
            fn[strlen(out_prefix)+1] = '2';
            bc.fp2 = gzopen(fn, "w");
            free(fn);
        } else {
            bc.fp1 = NULL;
            bc.fp2 = NULL;
        }
        ARRAY_PUSH(&bcs, BcRec, bc);
        k = kh_put(str, h, strdup(bc_seq), &ret);
        if (num_mismatches == 0) {
            kh_val(h, k) = bcs.nextfree - 1;
            //printf("setting %s to %lu (%s %s)\n", bc_seq, bcs.nextfree - 1, bcs.elems[bcs.nextfree - 1].seq, bcs.elems[bcs.nextfree - 1].id);
        } else {
            for (i=0; i<strlen(bc_seq); i++) {
                strcpy(bc_seq_cpy, bc_seq);
                for (j=0; j<dna_alpha_len; j++) {
                    bc_seq_cpy[i] = DNA_ALPHA[j];
                    k = kh_put(str, h, strdup(bc_seq_cpy), &ret);
                    kh_val(h, k) = bcs.nextfree - 1;
                    //printf("setting %s to %lu (%s %s)\n", bc_seq_cpy, bcs.nextfree - 1, bcs.elems[bcs.nextfree - 1].seq, bcs.elems[bcs.nextfree - 1].id);
                }
            }
        }
    }

    fclose(fp);

    if (bc_len == -1) {
        fprintf(stderr, "Error: could not find any barcodes in file %s\n", argv[optind]);
        return -1;
    }

    fp1 = gzopen(argv[optind+1], "r");
    seq1 = kseq_init(fp1);
    fp2 = gzopen(argv[optind+2], "r");
    seq2 = kseq_init(fp2);

    while (kseq_read(seq1) >= 0) {
        strncpy(bc_seq, seq1->seq.s, bc_len);
        k = kh_get(str, h, bc_seq);
        kseq_read(seq2);
        strncpy(bc_seq, seq2->seq.s, bc_len);
        k2 = kh_get(str, h, bc_seq);
        if (k != kh_end(h) || k2 != kh_end(h)) {
            bc_idx = k2 != kh_end(h) ? kh_val(h, k2) : kh_val(h, k);
            if (!only_count) {
                gzprintf(bcs.elems[bc_idx].fp1, "@%s %s\n%s\n+\n%s\n"
                         , seq1->name.s
                         , seq1->comment.s
                         , seq1->seq.s+bc_len+num_spacer_bases
                         , seq1->qual.s+bc_len+num_spacer_bases);
                gzprintf(bcs.elems[bc_idx].fp2, "@%s %s\n%s\n+\n%s\n"
                         , seq2->name.s
                         , seq2->comment.s
                         , seq2->seq.s+bc_len+num_spacer_bases
                         , seq2->qual.s+bc_len+num_spacer_bases);
            }
            bcs.elems[bc_idx].num_found += 2;
        } else {    
            num_undetermined += 2;
        }
    }
    
    gzclose(fp1);
    gzclose(fp2);
    kseq_destroy(seq1);
    kseq_destroy(seq2);

    for (i=0; i<bcs.nextfree; i++) {
        printf("%s\t%s\t%u\n", bcs.elems[i].id, bcs.elems[i].seq, bcs.elems[i].num_found); 
        if (!only_count) {
            gzclose(bcs.elems[i].fp1);
            gzclose(bcs.elems[i].fp2);
        }
    }

    printf("UNDETERMINED\tNONE\t%u\n", num_undetermined);

    ARRAY_FREE(&bcs);
    kh_destroy(str, h);

    fprintf(stderr, "[CPU time: %.2f sec]\n", (float)(clock() - t) / CLOCKS_PER_SEC);

    return 0;
}

int split_se(int argc, char *argv[], char *progname)
{
    int c, bc_len = -1, ret, i, j, bc_idx, only_count = 0;
    unsigned num_mismatches = DEFAULT_NUM_MISMATCHES, 
             num_spacer_bases = DEFAULT_NUM_SPACER_BASES, 
             dna_alpha_len = strlen(DNA_ALPHA),
             num_undetermined = 0;
    char *out_prefix = NULL, *fn = NULL, **sptr,
         bc_id[1024], bc_seq[1024], bc_seq_cpy[1024]; // hello, buffer overflow
    clock_t t = clock();
    BcRec bc;
    ArrayBcRec bcs;
    FILE *fp;
    kseq_t *seq;
    khash_t(str) *h = kh_init(str);
    khint_t k;
    gzFile *gzFp;

    ARRAY_INIT(&bcs, BcRec, 1000);

    while ((c = getopt(argc, argv, "m:s:o:c")) >= 0) {
        switch (c) {
            case 'm': if (sscanf(optarg, "%u", &num_mismatches) != 1) {
                          fprintf(stderr, "Error: option -m expects unsigned int\n");
                          return -1;
                      }
                      break;
            case 's': if (sscanf(optarg, "%u", &num_spacer_bases) != 1) {
                          fprintf(stderr, "Error: option -s expects unsigned int\n");
                          return -1;
                      }
                      break;
            case 'o': out_prefix = strdup(optarg);
                      break;
            case 'c': only_count = 1;
                      break;
        }
    }

    if (optind + 2 != argc) {
        print_se_usage(progname);
        return -1;
    }

    if (num_mismatches != 0 && num_mismatches != 1) {
        fprintf(stderr, "Error: argument -m has to be 0 or 1\n");
        return -1;
    }

    if (out_prefix == NULL) {
        out_prefix = strdup(DEFAULT_OUTPUT_PREFIX);
    }

    for (sptr = argv+optind; sptr-argv<argc; sptr++) {
        if (access(*sptr, F_OK) == -1) {
            fprintf(stderr, "Error: file %s does not exist\n", *sptr);
            return -1;
        }
    }

    fprintf(stderr, "[barcode file: %s]\n", argv[optind]);
    fprintf(stderr, "[fastq file: %s]\n", argv[optind+1]);
    fprintf(stderr, "[number of mismatches allowed: %u]\n", num_mismatches);
    fprintf(stderr, "[number of spacer bases: %u]\n", num_spacer_bases);
    fprintf(stderr, "[output prefix: %s]\n", out_prefix);
    fprintf(stderr, "[only count: %s]\n", only_count ? "true" : "false");

    /* read barcode file */
    if ((fp = fopen(argv[optind], "r")) == NULL) {
        fprintf(stderr, "Error: cannot open barcode file %s\n", argv[optind]);
        return -1;
    } 

    while (fscanf(fp, "%s %s", bc_id, bc_seq) == 2) {
        bc_len = strlen(bc_seq);
        bc.id = strdup(bc_id);
        bc.seq = strdup(bc_seq);
        bc.num_found = 0;
        if (!only_count) {
            fn = (char*)calloc(strlen(out_prefix) + 1 + strlen(bc_id) + 6 + 1, sizeof(char));
            strcpy(fn, out_prefix);
            strcat(fn, "_");
            strcat(fn, bc_id);
            strcat(fn, ".fq.gz");
            bc.fp1 = gzopen(fn, "w");
        } else {
            bc.fp1 = NULL;
        }
        bc.fp2 = NULL;
        free(fn);
        ARRAY_PUSH(&bcs, BcRec, bc);
        k = kh_put(str, h, strdup(bc_seq), &ret);
        if (num_mismatches == 0) {
            kh_val(h, k) = bcs.nextfree - 1;
        } else {
            for (i=0; i<strlen(bc_seq); i++) {
                strcpy(bc_seq_cpy, bc_seq);
                for (j=0; j<dna_alpha_len; j++) {
                    bc_seq_cpy[i] = DNA_ALPHA[j];
                    k = kh_put(str, h, strdup(bc_seq_cpy), &ret);
                    kh_val(h, k) = bcs.nextfree - 1;
                }
            }
        }
    }

    fclose(fp);

    if (bc_len == -1) {
        fprintf(stderr, "Error: could not find any barcodes in file %s\n", argv[optind]);
        return -1;
    }

    gzFp = gzopen(argv[optind+1], "r");
    seq = kseq_init(gzFp);

    while (kseq_read(seq) >= 0) {
        strncpy(bc_seq, seq->seq.s, bc_len);
        k = kh_get(str, h, bc_seq);
        if (k != kh_end(h)) {
            bc_idx = kh_val(h, k);
            if (!only_count) {
                gzprintf(bcs.elems[bc_idx].fp1, "@%s %s\n%s\n+\n%s\n"
                         , seq->name.s
                         , seq->comment.s
                         , seq->seq.s+bc_len+num_spacer_bases
                         , seq->qual.s+bc_len+num_spacer_bases);
            }
            bcs.elems[bc_idx].num_found += 1;
        } else {
            num_undetermined += 1;
        }
    }

    gzclose(gzFp);
    kseq_destroy(seq);

    for (i=0; i<bcs.nextfree; i++) {
        printf("%s\t%s\t%u\n", bcs.elems[i].id, bcs.elems[i].seq, bcs.elems[i].num_found);
        if (!only_count) {
            gzclose(bcs.elems[i].fp1);
        }
    }

    printf("UNDETERMINED\tNONE\t%u\n", num_undetermined);

    ARRAY_FREE(&bcs);
    kh_destroy(str, h);

    fprintf(stderr, "[CPU time: %.2f sec]\n", (float)(clock() - t) / CLOCKS_PER_SEC);

    return 0;
}

int main(int argc, char *argv[])
{
    int ret;

    if (argc == 1) {
        print_main_usage(basename(argv[0]));
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "pe") == 0) {
        ret = split_pe(argc-1, argv+1, basename(argv[0])) != 0;
    } else if (strcmp(argv[1], "se") == 0) {
        ret = split_se(argc-1, argv+1, basename(argv[0])) != 0;
    } else {
        fprintf(stderr, "\nError: unknown mode \"%s\"\n", argv[1]);
        print_main_usage(argv[0]);
        return EXIT_FAILURE;
    }

    return ret == 0 ? EXIT_SUCCESS: EXIT_FAILURE;
}
