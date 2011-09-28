/**
   Copyright (c) 2010, Pacific Biosciences of California, Inc.

   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

   * Neither the name of Pacific Biosciences nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED

   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS
   BE LIABLE FOR ANY

   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
   GOODS OR SERVICES;

   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
   LIABILITY, OR TORT

   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/

#include <Rinternals.h>    
#include <R.h>

/**
   Returns the insertion point in a vector of values which preserves
   sorted order. If val exists in vec then the return value is to the
   left of the position of val in vec.
**/
 int _bisect(int val, int* vec,  int l,  int r) {
    if ((r - l) <= 3) {
	// a cheat :)
	while (((r - 1) >= 0) && (vec[r-1] >= val)) {
	    r--;
	}
	return r;
    }
    else {
	int i = (l + r)/2;

	if (vec[i] > val) {
	    return _bisect(val, vec, l, i);
	} else {
	    return _bisect(val, vec, i, r);
	}
    }
}

int left_bin_search(int val, int* vec,  int l) {
    int s = _bisect(val, vec, 0, l);
    int v = -1;
    
    if (s == 0) {
	return s;
    }
    else if (s == l) {
	v = vec[--s];
    } 
    else {
	v = vec[s];
    }

    if (v > val) {
	s--;
    }
    while (s > 0 && vec[s-1] == vec[s])	s--;

    return s;
}

int right_bin_search(int val, int* vec,  int l) {
    int s = _bisect(val, vec, 0, l);
    
    if (s == l)
	return(s);
    while (s + 1 < l && vec[s + 1] == val)
	s++;
    return s;
}

SEXP R_bisect(SEXP val, SEXP vec) {
    return ScalarInteger( _bisect(INTEGER(val)[0], INTEGER(vec), 0,  length(vec)));
}

int _overlaps(int tstart, int tend, int start, int end) {
    if (end >= tstart && end <= tend)
	return 1;
    if (start >= tstart && start <= tend)
	return 1;
    if (start >= tstart && end <= tend)
	return 1;
    if (start <= tstart && end >= tend)
	return 1;
    return 0;
}

SEXP PBR_indices_within_range(SEXP _t_start, SEXP _t_end, SEXP _n_back, 
			      SEXP _n_over, SEXP _start, SEXP _end) {
    int* t_start = INTEGER(_t_start);
    int* t_end   = INTEGER(_t_end);
    int* n_back  = INTEGER(_n_back);
    int* n_over  = INTEGER(_n_over);
    int start    = INTEGER(_start)[0];
    int end      = INTEGER(_end)[0];
    int len      = length(_t_start);

    int lm, rm, i, s = 0;

    lm = left_bin_search(start, t_start, len);

    //    Rprintf("leftmost bin: %d\n", lm);

    lm = lm - (n_back[lm] > lm ? lm : n_back[lm]);

    //    Rprintf("leftmost bin: %d\n", lm);

    rm = right_bin_search(end, t_start, len);
    
    //    Rprintf("rightmost bin: %d\n", rm);
    
    /** First, I determine how much to allocate. **/
    for (i = lm; i <= rm; i++) {
	// This is unclean - my searches don't exactly want
	// to fulfill their contracts.
	if (i >= 0 && i < len) 
	    s += _overlaps(t_start[i], t_end[i], start, end);
    }

    SEXP res = R_NilValue;
    PROTECT(res = allocVector(INTSXP, s));
    int* rptr = INTEGER(res);

    /* Rprintf("Allocated vectors of size: %d, writing in range (%d, %d)\n",  */
    /* 	    s, 0, rm - lm - 1); */

    int k = 0;
    for (i = lm; i <= rm; i++) {
	// Ibid, unclean.
	if (i >= 0 && i < len) {
	    if (_overlaps(t_start[i], t_end[i], start, end)) {
		rptr[k++] = i;
	    }
	}
    }
    
    UNPROTECT(1);
    
    return res;
}

#define _MIN_(a, b) (((a < b) ? a : b))
#define _MAX_(a, b) (((a < b) ? b : a))

SEXP PBR_coverage_within_range(SEXP _t_start, SEXP _t_end, SEXP _n_back, 
			       SEXP _n_over, SEXP _start, SEXP _end) {
    int* t_start  = INTEGER(_t_start);
    int* t_end    = INTEGER(_t_end);
    int start     = INTEGER(_start)[0];
    int end       = INTEGER(_end)[0];
    int len       = length(_t_start);
    SEXP _res     = R_NilValue; 
    SEXP _indices = PBR_indices_within_range(_t_start, _t_end, _n_back,
					     _n_over, _start, _end);
    int* indices = INTEGER(_indices);

    // Rprintf("Got indices, length: %d\n", length(_indices));

    int i, j, min_s, max_e;

    PROTECT(_res = allocVector(INTSXP, end - start + 1)); 
    int* res = INTEGER(_res);
    for (i = 0; i < length(_res); i++)
	res[i] = 0;

    // Rprintf("Allocated: %d\n", end - start + 1);

    for (i = 0; i < length(_indices); i++) {
	min_s = _MAX_(t_start[indices[i]], start) - start;
	max_e = _MIN_(t_end[indices[i]], end) - start;
	
	// Rprintf("(%d, %d)\n", min_s, max_e);

	if (max_e >= min_s) {
	    for (j = min_s; j <= max_e; j++) {
		res[j] = res[j] + 1;
	    }
	}
    }
    UNPROTECT(1);

    return(_res);
}

typedef struct _PBR_tnode_s {
    struct _PBR_tnode_s* parent;
    struct _PBR_tnode_s** children;
    char letter; 
    int count;
    int nchildren;
    int index;
} PBR_node;


PBR_node* _PBR_make_node(PBR_node* parent, char character) {
    PBR_node* node = (PBR_node*) R_alloc(1, sizeof(PBR_node));
    node->letter = character;
    node->parent = parent;
    node->count = 0;
    node->children = NULL;
    node->nchildren = 0;
    node->index = 0;
    return node;
}

void _PBR_init_tree(PBR_node* tree_node, const char* alphabet, int a_length, int k) {
    int i;

    if (k > 0) {
	tree_node->children = (PBR_node**) R_alloc(a_length, sizeof(PBR_node*));
	tree_node->nchildren = a_length;

	for (i = 0; i < a_length; i++) {
	    tree_node->children[i] = _PBR_make_node(tree_node, alphabet[i]);
	    _PBR_init_tree(tree_node->children[i], alphabet, a_length, k - 1);
	}
    }
}

void _PBR_init_labels(PBR_node* tree, int* count) {
    int i;

    if (tree->children == NULL) {
	tree->index = (*count)++;
    }
    else {
	for (i = 0; i < tree->nchildren; i++) {
	    _PBR_init_labels(tree->children[i], count);
	}
    }
}

void _PBR_write_tree(PBR_node* tree, int* vector, int* level) {
    int i;

    if (tree->children == NULL) {
	vector[(*level)++] = tree->count;
    }
    else {
	for (i = 0; i < tree->nchildren; i++) {
	    _PBR_write_tree(tree->children[i], vector, level);
	}
    }
}

void _PBR_zero_tree(PBR_node* tree) {
    int i;

    if (tree->children == NULL) {
	tree->count = 0;
    }
    else {
	for (i = 0; i < tree->nchildren; i++) {
	    _PBR_zero_tree(tree->children[i]);
	}
    }
}

void _PBR_print_tree(PBR_node* tree) {
    int i;

    if (tree->children == NULL) {
	Rprintf("node: %d\n", tree->index);
    }
    else {
	for (i = 0; i < tree->nchildren; i++) {
	    _PBR_print_tree(tree->children[i]);
	}
    }
}

SEXP PBR_k_mer_tabulator(SEXP lst, SEXP _alphabet, SEXP _k) {
    int i,j,l,m;
    int location = 0;

    const char* alphabet = CHAR(STRING_ELT(_alphabet, 0));
    int k = INTEGER(_k)[0];
    int alen = strlen(alphabet);
    int tbl_length = (int) pow(alen, k);
    int llen = LENGTH(lst);

    /* Rprintf("Alphabet: %s\n", alphabet); */
    /* Rprintf("Alphabet Length: %d\n", alen); */
    /* Rprintf("K: %d\n", k); */
    /* Rprintf("List Length v:%d\n", llen); */

    SEXP r_lst;
    PROTECT(r_lst = allocVector(INTSXP, tbl_length*llen));

    PBR_node* tree = _PBR_make_node(NULL, '|');
    _PBR_init_tree(tree, alphabet, strlen(alphabet), k);
    
    for (m = 0; m < llen; m++) {
    	SEXP v   = VECTOR_ELT(lst, m);
    	int vlen = LENGTH(v);
    
    	PBR_node* current = tree;
    	for (i = 0; i < vlen - k + 1; i++) {
    	    for (j = 0; j < k; j++) {
    		for (l = 0; l < alen; l++) {
    		    if (current->children[l]->letter == CHAR(STRING_ELT(v, i + j))[0]) {
    			current = current->children[l];
    			break;
    		    }
    		}
    	    }
    	    current->count++;
    	    current = tree;
    	}
    	location = 0;
    	_PBR_write_tree(tree, &(INTEGER(r_lst)[m*tbl_length]), &location);
	_PBR_zero_tree(tree);
    }
    UNPROTECT(1);

    return r_lst;
}

SEXP PBR_k_mer_labeler(SEXP lst, SEXP _alphabet, SEXP _k) {
    int i,j,l,m;
    int location = 0;

    const char* alphabet = CHAR(STRING_ELT(_alphabet, 0));
    int k = INTEGER(_k)[0];
    int alen = strlen(alphabet);
    int tbl_length = (int) pow(alen, k);
    int llen = LENGTH(lst);

    /* Rprintf("Alphabet: %s\n", alphabet); */
    /* Rprintf("Alphabet Length: %d\n", alen); */
    /* Rprintf("K: %d\n", k); */
    /* Rprintf("List Length v:%d\n", llen); */

    SEXP r_lst;
    PROTECT(r_lst = allocVector(VECSXP, llen));

    PBR_node* tree = _PBR_make_node(NULL, '|');
    _PBR_init_tree(tree, alphabet, strlen(alphabet), k);

    int labels = 0;
    _PBR_init_labels(tree, &labels);
    //    _PBR_print_tree(tree);

    for (m = 0; m < llen; m++) {
    	SEXP v      = VECTOR_ELT(lst, m);
    	int vlen    = LENGTH(v);
	SEXP r_vec  = allocVector(INTSXP, vlen);
	int* _r_vec = INTEGER(r_vec);
    
    	PBR_node* current = tree;
    	for (i = 0; i < vlen - k + 1; i++) {
    	    for (j = 0; j < k; j++) {
    		for (l = 0; l < alen; l++) {
    		    if (current->children[l]->letter == CHAR(STRING_ELT(v, i + j))[0]) {
    			current = current->children[l];
    			break;
    		    }
    		}
    	    }
	    _r_vec[i] = current->index;
    	    current   = tree;
    	}
	for (i = vlen - 1; i > vlen - k; i--)
	    _r_vec[i] = NA_INTEGER;

	SET_VECTOR_ELT(r_lst, m, r_vec);
    }
    UNPROTECT(1);

    return r_lst;
}



#define IS_DELETION(I) (((I) >> 4) == 0)
#define IS_INSERTION(I) (((I) & 0x0f) == 0)

SEXP PBR_advance_time(SEXP start_lst, SEXP aln_lst) {
    int m, llen, slen, j, del_count, k, i;
    SEXP res, start_vec, advance_vec, aln_vec;
    double r_time, l_time;
    double* _advance_vec;
    double* _start_vec;
    int* _aln_vec;

    llen = LENGTH(start_lst);
    PROTECT(res = allocVector(VECSXP, llen));
    
    for (m = 0; m < llen; m++) {
    	start_vec    = VECTOR_ELT(start_lst, m);
	aln_vec      = VECTOR_ELT(aln_lst, m);
    	slen         = LENGTH(start_vec);
	advance_vec  = allocVector(REALSXP, slen);

	_advance_vec = REAL(advance_vec);
	_start_vec   = REAL(start_vec);
	_aln_vec     = INTEGER(aln_vec);
	
	_advance_vec[0] = NA_REAL;

	j = slen - 1;

	while (j > 0) {
	    r_time    = _start_vec[j];
	    del_count = 0;
	    i         = j - 1;

	    while (IS_INSERTION(_aln_vec[i]) && i >= 0) {
		_advance_vec[i--] = NA_REAL;
	    }
	    while (IS_DELETION(_aln_vec[i]) && i >= 0) {
		del_count++;
		i--;
	    }
	    l_time = _start_vec[i];

	    if (del_count > 0) {
		for (k = i + 1; k <= j; k++)
		    _advance_vec[k] = (r_time - l_time)/((double) del_count);
	    } else {
		_advance_vec[j] = (r_time - l_time);
	    }
	    j = i;
	}

	/** Total Hack because I have seen gaps starting alignments. **/
	for (j = 0; j < slen; j++) 
	    _advance_vec[j] = (_advance_vec[j] < 0) ? NA_REAL : _advance_vec[j];

	SET_VECTOR_ELT(res, m, advance_vec);
    }
    UNPROTECT(1);

    return res;
}


SEXP PBR_get_template_position(SEXP lst, SEXP strands, SEXP starts, SEXP ends) {
    int i,m,tpos,llen = LENGTH(lst);
    SEXP r_lst;
    PROTECT(r_lst = allocVector(VECSXP, llen));

    for (m = 0; m < llen; m++) {
    	SEXP aln    = VECTOR_ELT(lst, m);
	int strand  = INTEGER(strands)[m];
    	int vlen    = LENGTH(aln);
	SEXP r_vec  = allocVector(INTSXP, vlen);
	int* _r_vec = INTEGER(r_vec);

	int start = INTEGER(starts)[m] - 1;
	int end   = INTEGER(ends)[m] + 1;

	if (strand == 1) {
	    tpos = end;
	    for (i = 0; i < vlen; i++) {
		if (! IS_INSERTION(INTEGER(aln)[i])) {
		    tpos--;
		}
		_r_vec[i] = tpos;
	    }
	}
	else {
	    tpos = start;
	    for (i = 0; i < vlen; i++) {
		int insert = IS_INSERTION(INTEGER(aln)[i]);
		if (! insert) {
		    tpos++;
		}
		_r_vec[i] = (insert) ? tpos + 1 : tpos;
	    }
	}
	SET_VECTOR_ELT(r_lst, m, r_vec);
    }

    UNPROTECT(1);
    return r_lst;
}

/**
 * All this is necessary for the consensus code and to process raw alignments.
 * 
 */
char PB_READ_MAP[] = { '-','-','-','*','-','*','*','*','-','*','*','*','*','*','*','-','A','A','A','*','A','*','*','*','A','*','*','*',
		       '*','*','*','A','C','C','C','*','C','*','*','*','C','*','*','*','*','*','*','C','*','*','*','*','*','*','*','*','*','*',
		       '*','*','*','*','*','*','G','G','G','*','G','*','*','*','G','*','*','*','*','*','*','G','*','*','*','*','*','*','*','*',
		       '*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*',
		       '*','*','*','*','*','*','*','*','T','T','T','*','T','*','*','*','T','*','*','*','*','*','*','T','*','*','*','*','*','*',
		       '*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*',
		       '*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*',
		       '*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','N','N','N','*','N',
		       '*','*','*','N','*','*','*','*','*','*','N' };
int RC_READ_MAP[] = {
    0,8,4,3,2,5,6,7,1,9,10,11,12,13,14,15,128,136,132,19,130,21,22,23,129,25,26,27,28,29,30,143,64,72,68,35,66,37,38,
    39,65,41,42,43,44,45,46,79,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,32,40,36,67,34,69,70,71,33,73,74,75,76,
    77,78,47,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,
    112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,16,24,20,131,18,133,134,135,17,137,138,139,140,141,142,
    31,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,
    173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,
    203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,
    233,234,235,236,237,238,239,240,248,244,243,242,245,246,247,241,249,250,251,252,253,254,255
};
#define BASECALL(code) (PB_READ_MAP[(code)])

void _PBR_reverse_complement(int* original, int* results, int length) {
    int k;
    for (k = 0; k < length; k++) {
	results[k] = RC_READ_MAP[original[length - k - 1]];
    }
}

typedef struct _PBR_base_node {
    struct _PBR_base_node* next;
    char base;
} PBR_base_node;

typedef struct _PBR_evt_node {
    struct _PBR_evt_node* next;
    struct _PBR_base_node* basecalls;
    int nbases;
} PBR_evt_node;

typedef struct _PBR_evt_list {
    struct _PBR_evt_node* list;
    int nevents;
} PBR_evt_list;

PBR_evt_list* _PBR_make_evt_list() {
    PBR_evt_list* list = (PBR_evt_list*) R_alloc(1, sizeof(PBR_evt_list));
    list->nevents = 0;
    list->list = NULL;
    return list;
}

void _print_events(PBR_evt_list* evt) {
    if (evt == NULL) {
	Rprintf("Event list is null\n");
	return;
    }

    PBR_evt_node* node = evt->list;
    while (node != NULL) {
	PBR_base_node* basecalls = node->basecalls;
	while (basecalls != NULL) {
	    Rprintf("%c", basecalls->base);
	    basecalls = basecalls->next;
	}
	Rprintf("[%d] ", node->nbases);
	node = node->next;
    }
}

SEXP PBR_compute_consensus(SEXP alignments, SEXP _refStart, SEXP _refEnd, SEXP _strands, SEXP _starts, 
			   SEXP _ends) {
    int i,j,k;
    
    int ref_start = INTEGER(_refStart)[0];
    int ref_end   = INTEGER(_refEnd)[0];
    int n_bases   = ref_end - ref_start + 1;
    
    PBR_evt_list** all_events = (PBR_evt_list**) R_alloc(n_bases, sizeof(PBR_evt_list*));
    for (i = 0; i < n_bases; i++) all_events[i] = _PBR_make_evt_list();
  
    for (i = 0; i < LENGTH(alignments); i++) {
    	SEXP _aln   = VECTOR_ELT(alignments, i);
	int vlen    = LENGTH(_aln);
	int strand  = INTEGER(_strands)[i];
	int start   = INTEGER(_starts)[i] - ref_start;
	int* aln    = INTEGER(_aln);
	
	if (strand == 1) {
	    int* rcaln = (int*) R_alloc(vlen, sizeof(int));
	    _PBR_reverse_complement(aln, rcaln, vlen);
	    aln = rcaln;
	}

	/** j indexes positions in the alignment. 
	    k indexes shifted positions in the reference. **/
	j = 0; k = start;
	while (j < vlen) {
	    if (k >= n_bases) {
		break;
	    }
	    if (k < 0) {
		if (! IS_INSERTION(aln[j])) 
		    k++;
		j++;
		continue;
	    }

	    /** get the event list corresponding to this reference position. **/
	    PBR_evt_list* events = all_events[k];
	    PBR_evt_node* current_node = (PBR_evt_node*) R_alloc(1, sizeof(PBR_evt_node));

	    if (events->list == NULL) {
		events->list = current_node;
		events->list->next = NULL;
		events->nevents = 1;
	    } else { 
		current_node->next = events->list;
		events->list = current_node;
		events->nevents++;
	    }
	    current_node = events->list;
	    
	    /** Initialize basecall for this position.**/
	    PBR_base_node* current_base = (PBR_base_node*) R_alloc(1, sizeof(PBR_base_node));
	    current_base->base = BASECALL(aln[j]);
	    current_base->next = NULL;
	    
	    current_node->basecalls = current_base;
	    current_node->nbases    = 1;

	    j++;
	    while (j < vlen && IS_INSERTION(aln[j])) { 
		PBR_base_node* next_base = (PBR_base_node*) R_alloc(1, sizeof(PBR_base_node));
		next_base->next = NULL;
		next_base->base = BASECALL(aln[j]);
		current_base->next = next_base;
		current_base = next_base;

		current_node->nbases++;
	    	j++;
	    }
	    k++;
	}
    }
    SEXP r_lst;
    SEXP bevents;
    int p_stack = 1;
    
    PROTECT(r_lst = allocVector(VECSXP, n_bases));
    
    for (i = 0; i < n_bases; i++) {
	PBR_evt_list* events = all_events[i];
	
	if (events != NULL && events->nevents > 0) {
	    bevents = allocVector(STRSXP, events->nevents);
	    SET_VECTOR_ELT(r_lst, i, bevents);
	 
	    j = 0;
	    PBR_evt_node* event = events->list;

	    while (event != NULL) {
		if (event->basecalls == NULL || event->nbases == 0) {
		    break;
		}
		PBR_base_node* basecall = event->basecalls;
		char* newstr = (char*) R_alloc(event->nbases + 1, sizeof(char));
		newstr[event->nbases] = '\0';
		
		k = 0;
		while(basecall != NULL) {
		    newstr[k] = basecall->base;
		    basecall  = basecall->next;
		    k++;
		}

		SEXP s = mkChar(newstr);
		SET_STRING_ELT(bevents, j, s);
		event = event->next;
		j++;
	    }
	    
	}
    }
    UNPROTECT(p_stack);
    return r_lst;
}


int _PBR_cstring_cmp(const void *a, const void *b) { 
    const char **ia = (const char **)a;
    const char **ib = (const char **)b;
    return strcmp(*ia, *ib);
} 

char* _choose_event(char** bevents, int n) {
    if (n == 1) {
	return bevents[0];
    } else {
	char* current_evt = bevents[0];
	char* max_evt = current_evt;
	int current_count = 1;
	int max_count = 1;
	qsort(bevents, n, sizeof(char*), _PBR_cstring_cmp);
	for (int k = 1; k < n; k++) {
	    if (strcmp(current_evt, bevents[k]) == 0) {
		if (++current_count > max_count) {
		    max_count++;
		    max_evt = current_evt;
		}
	    } else {
		current_evt = bevents[k];
		current_count = 1;
	    }
	}
	return max_evt;
    }
}


struct _char_ll {
    char* dta;
    struct _char_ll* nxt;
} typedef char_ll;

char_ll* make_char_ll(char* d) {
    char_ll* newd = (char_ll*) R_alloc(1, sizeof(char_ll));
    newd->dta = d;
    newd->nxt = NULL;
    return newd;
}

char* _choose_event_2(char** bevents, int n) {
    if (n == 1) {
	return bevents[0];
    } else if (n == 2) {
	if (unif_rand() > .5) {
	    return bevents[0];
	} else {
	    return bevents[1];
	}
    } else {
	char* current_evt = bevents[0];
	char_ll* max_lst  = NULL;
	int current_count = 1;
	int max_count = 1;
	int nmaxes = 0;
	qsort(bevents, n, sizeof(char*), _PBR_cstring_cmp);
	for (int k = 1; k < n; k++) {
	    if (strcmp(current_evt, bevents[k]) == 0) {
		++current_count;
		if (current_count == max_count) {
		    nmaxes++;
		    if (max_lst == NULL) {
			max_lst = make_char_ll(current_evt);
		    } else {
			max_lst->nxt = make_char_ll(current_evt);
		    }
		}
		if (current_count > max_count) {
		    nmaxes = 1;
		    max_count = current_count;
		    max_lst = make_char_ll(current_evt);
		}
	    } else {
		current_evt = bevents[k];
		current_count = 1;
	    }
	}
	if (nmaxes == 0) {
	    double u = unif_rand();
    	    double z = 1.0 / n;
    	    for (int k = 0; k < n; k++) {
    		if (u > k*z && u <= (k+1)*z) {
    		    return bevents[k];
    		}
	    }
	    // This shouldn't occur. 
	    return bevents[0];
	} else if (nmaxes == 1) {
	    return max_lst->dta;
	} else {
	    double u = unif_rand();
    	    double z = 1.0 / nmaxes;
    	    for (int k = 0; k < nmaxes; k++) {
    		if (u > k*z && u <= (k+1)*z) {
    		    return max_lst->dta;
    		} else {
    		    if (max_lst->nxt == NULL) {
    			return max_lst->dta;
    		    } else {
    			max_lst = max_lst->nxt;
    		    }
    		}
	    }
	}
    }
}

SEXP PBR_compute_consensus2(SEXP alignments, SEXP _refStart, SEXP _refEnd, SEXP _strands, SEXP _starts, 
			    SEXP _ends) {
    int i,j,k;
    int ref_start = INTEGER(_refStart)[0];
    int ref_end   = INTEGER(_refEnd)[0];
    int n_bases   = ref_end - ref_start + 1;

    // each base gets an event-list - a linked list. 
    PBR_evt_list** all_events = (PBR_evt_list**) R_alloc(n_bases, sizeof(PBR_evt_list*));
    for (i = 0; i < n_bases; i++) all_events[i] = _PBR_make_evt_list();
  
    // go through each alignment and add events to each event list. 
    for (i = 0; i < LENGTH(alignments); i++) {
    	SEXP _aln   = VECTOR_ELT(alignments, i);
	int vlen    = LENGTH(_aln);
	int strand  = INTEGER(_strands)[i];
	int start   = INTEGER(_starts)[i] - ref_start;
	int* aln    = INTEGER(_aln);
	
	// reverse complement if need be - this might be a point of optimization.
	if (strand == 1) {
	    int* rcaln = (int*) R_alloc(vlen, sizeof(int));
	    _PBR_reverse_complement(aln, rcaln, vlen);
	    aln = rcaln;
	}

	/** j indexes positions in the alignment. 
	    k indexes shifted positions in the reference. **/
	j = 0; k = start;
	while (j < vlen) {
	    if (k >= n_bases) {
		break;
	    }
	    if (k < 0) {
		if (! IS_INSERTION(aln[j])) 
		    k++;
		j++;
		continue;
	    }

	    /** get the event list corresponding to this reference position. **/
	    PBR_evt_list* events = all_events[k];
	    PBR_evt_node* current_node = (PBR_evt_node*) R_alloc(1, sizeof(PBR_evt_node));

	    if (events->list == NULL) {
		events->list = current_node;
		events->list->next = NULL;
		events->nevents = 1;
	    } else { 
		current_node->next = events->list;
		events->list = current_node;
		events->nevents++;
	    }
	    current_node = events->list;
	    
	    /** Initialize basecall for this position.**/
	    PBR_base_node* current_base = (PBR_base_node*) R_alloc(1, sizeof(PBR_base_node));
	    current_base->base = BASECALL(aln[j]);
	    current_base->next = NULL;
	    current_node->basecalls = current_base;
	    current_node->nbases    = 1;

	    j++;
	    while (j < vlen && IS_INSERTION(aln[j])) { 
		PBR_base_node* next_base = (PBR_base_node*) R_alloc(1, sizeof(PBR_base_node));
		next_base->next = NULL;
		next_base->base = BASECALL(aln[j]);
		current_base->next = next_base;
		current_base = next_base;
		current_node->nbases++;
	    	j++;
	    }
	    k++;
	}
    }
    
    // Now we are setting up the return value of the function. 
    SEXP r_lst;
    SEXP bevents;
    PROTECT(r_lst = allocVector(STRSXP, n_bases));
    
    for (i = 0; i < n_bases; i++) {
	PBR_evt_list* events = all_events[i];
	if (events != NULL && events->nevents > 0) {
	    char** bevents = (char**) R_alloc(events->nevents, sizeof(char*));
	    PBR_evt_node* event = events->list;
	    j = 0;
	    while (event != NULL) {
		if (event->basecalls == NULL || event->nbases == 0) {
		    break;
		}
		PBR_base_node* basecall = event->basecalls;
		char* newstr = (char*) R_alloc(event->nbases + 1, sizeof(char));
		newstr[event->nbases] = '\0';
		k = 0; while(basecall != NULL) {
		    newstr[k] = basecall->base;
		    basecall  = basecall->next;
		    k++;
		}
		bevents[j] = newstr;
		event = event->next;
		j++;
	    }
	    char* max_evt = _choose_event_2(bevents, events->nevents);
	    SET_STRING_ELT(r_lst, i, mkCharLen(max_evt, strlen(max_evt)));
	} else {
	    SET_STRING_ELT(r_lst, i, NA_STRING);
	}
    }

    UNPROTECT(1);
    return r_lst;
}
