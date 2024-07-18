#include <iostream>
#include <random>


/* #include <flint/flint.h> */
/* #include <flint/fmpz.h>		/\* large integers *\/ */
/* #include <flint/fq.h>		/\* finite fields *\/ */

/* #include "flint/fq_poly.h"	/\* pol. in finite fields *\/ */
/* #include "flint/fmpz_poly.h"	/\* pol. in integers *\/ */
/* #include "flint/fmpz_vec.h"	/\* vectors integers *\/ */
/* #include "flint/fq_vec.h"	/\* vectors finite fields *\/ */
/* #include "flint/perm.h"		/\* permutations *\/ */
/* #include "flint/fq_mat.h"	/\* matrix / finite fields *\/ */


#include "tools.hpp"



using namespace std;



void
fisher_yates(int* perm, const int n) {
    int j;
  
    std::random_device rd;
    std::mt19937 rand_gen(rd());
  
  
    for (int i = 0; i < n; ++i) {
	std::uniform_int_distribution<std::mt19937::result_type> dis(0, i);
	j = dis(rand_gen);
	if (i != j) {
	    perm[i] = perm[j];
	}
	perm[j] = i;
    }
}




slong*
random_indices(const slong len, flint_rand_t state) {
    slong* perm = _perm_init(len);
    _perm_randtest(perm, len, state);
    return(perm);
}

/* check repetitions in a naive way */
int
int_check_repeat(const int* a, const int len) {
    for (int i = 0; i < len; i++){
	for (int j = i+1; j < len; j++) {
	    if (a[i] == a[j]) return 1;
	}
    }
    return 0;
}




/*  computes t random and pairwise distinct indexes in [1,n] as in Classic
    McEliece
    draw random τ > t elts in [1, N[ and take the t first ones that are in [1,n]
    return the number of elts in [1,n] effectively computed
    Typically in CM, N = 2^m
*/
int
cm_random_indices(int* res, const int n, const int t, const int N, const int tau) {
    int ind = 0;
    int count_loop = 0;
    int a;

    // rand with mersenne
    std::random_device rd;
    std::mt19937 rand_gen(rd());
    std::uniform_int_distribution<std::mt19937::result_type> dis(0, N-1);
    
    while (ind < t && count_loop < tau) {
	a = dis(rand_gen);
	if (a < n) {
	    res[ind] = a;
	    ind++;
	}
	count_loop++;
    }
    return ind;
}


/* **************************************************************************** */
/*                              ERROR GENERATION                                */
/* **************************************************************************** */

/* computes e as in Classic McEliece, more or less */
void
cm_gen_e(int* e, const int n, const int t, const int N, const int tau) {
    int len, repet = 1;
    int inds[t];
    while (repet) {
	len = cm_random_indices(inds, n, t, N, tau); /* compute random set of indices */
	if (len == t) {				 /* check if we obtained enough indices */
	    repet = int_check_repeat(inds, t);
	}
    }
    for (int i = 0; i <t; i++) {
	e[inds[i]] = 1;
    }
}




/* **************************************************************************** */
/*                              FINITE FIELDS                                   */
/* **************************************************************************** */


/*
  Go from finite field vec to integral vector
  Assume a is 0 or 1 */
void
_fq_vec_2_int(int* res, const fq_struct* a, const int len, const fq_ctx_t ctx) {
    fmpz_t tmp;
    fmpz_init(tmp);
    int b;
    for (int i = 0 ; i < len; ++i) {
	b = fq_get_fmpz(tmp, &a[i], ctx);
	res[i] = fmpz_get_ui(tmp);
    }
}


/*
  Go from integral vector  to finite field vec 
  Assume a is 0 or 1 */
void
_int_vec_2_fq(fq_struct* res, const int* a, const int len, const fq_ctx_t ctx) {
    for (int i = 0 ; i < len; ++i) {
	fq_set_ui(&res[i], a[i], ctx);
    }
}


/* check repetitions in finite fields in a naive way */
int
fq_check_repeat(const fq_struct* a, const int len, const fq_ctx_t ctx) {
    for (int i = 0; i < len; i++){
	for (int j = i+1; j < len; j++) {
	    if (fq_equal(&a[i], &a[j], ctx))  return 1; 
	}
    }
    return 0;
}



/* generate a random vector whose coeffs are in a finite field and all distincts */
void
fq_vec_rand_distinct(fq_struct* res, const int len, const fq_ctx_t ctx, flint_rand_t state) {
    _fq_vec_randtest(res, state, len, ctx); /* NOT GOOD AT ALL → change it */
    while (fq_check_repeat(res, len, ctx)) {
	_fq_vec_randtest(res, state, len, ctx);
    }
}

/* generate a random vector whose coeffs are in a finite field and all distincts
   WORK ONLY for *BINARY* fields
   NOT THE MOST ELEGANT WAY -- MAYBE CHANGE */
void
fq_vec_rand_distinct_2(fq_struct* res, const int len, const fq_ctx_t ctx, flint_rand_t state) {
    slong m = fq_ctx_degree(ctx);
    slong d = 1 << m;
    slong* inds = random_indices(d, state); /* here is the binary stuff */
    fmpz_t tmp;
    fmpz_init(tmp);
    for (int i = 0; i < len; i++) {
	fmpz_set_ui(tmp, inds[i]);
	fq_bit_unpack(&res[i], tmp, 1, ctx); 
    }
    fmpz_clear(tmp);
}



/* **************************************************************************** */
/*                           POLYNOMIAL MANIPUTATION                            */
/* **************************************************************************** */

/* check if f has a root in alpha */
int
fq_poly_eval_zero(const fq_poly_t f, const fq_struct *alpha, const int len,
		  const fq_ctx_t ctx) {
    int b = 0;
    fq_t a;
    fq_init(a, ctx);
    for (int i = 0; i < len; i++) {
	fq_poly_evaluate_fq(a, f, &alpha[i], ctx);
	b = fq_is_zero(a, ctx);
	if (b) break;
    }
    fq_clear(a, ctx);
    return (b);
}


/* computes the unique polynomial */
void
fq_poly_interpolate(fq_poly_t res, const fq_struct* alpha, const fq_struct*  beta, const int len,
		    const fq_ctx_t ctx) {
    /* declare and initialise local variables */
    fq_poly_t A, B, tmp;
    fq_poly_init(A, ctx);
    fq_poly_init(B, ctx);
    fq_poly_init(tmp, ctx);
    fq_poly_gen(tmp, ctx);
    fq_poly_one(A, ctx);

    fq_t q; fq_init(q, ctx);
    fq_t a; fq_init(a, ctx);

    for (int i = 0; i < len; ++i) {
	fq_neg(a, &alpha[i], ctx);
	fq_poly_set_coeff(tmp, 0, a, ctx);
	fq_poly_mul(A, A, tmp, ctx);
    }

    fq_poly_zero(res, ctx);	/* put res to zero */
    for (int i = 0; i < len; i++) {
	fq_neg(a, &alpha[i], ctx);
	fq_poly_set_coeff(tmp, 0, a, ctx);
	fq_poly_div(B, A, tmp, ctx);
	fq_poly_evaluate_fq(q, B, &alpha[i], ctx);
	fq_div(q, &beta[i], q, ctx);
	fq_poly_scalar_mul_fq(B, B, q, ctx);
	fq_poly_add(res, res, B, ctx);
    }

    /* clear memory */
    fq_poly_clear(tmp, ctx);
    fq_poly_clear(A, ctx);
    fq_poly_clear(B, ctx);
    fq_clear(q, ctx);
    fq_clear(a, ctx);
}


/* compute an irreducible polynomial of degree "deg" as in CM ie so that it generates a Goppa code
   together with the roots alpha Γ(alpha, res)  
*/
void
cm_fq_poly_irr_pol(fq_poly_t& res, const int deg, const fq_struct* alpha, const int len,
		   const fq_ctx_t ctx, flint_rand_t state) {
    fq_poly_randtest_irreducible(res, state, deg+1, ctx);
    while (fq_poly_eval_zero(res, alpha, len, ctx)) {
	fq_poly_randtest_irreducible(res, state, deg+1, ctx);
    }
}


/* computes early abort extended gcd of a and b in finite field defined by context ctx */
void
xgcd_abort(fq_poly_t u, fq_poly_t v, fq_poly_t d, const fq_poly_t a, const fq_poly_t b,
	   const slong k, const fq_ctx_t ctx) {
    /* initialisation */
    fq_poly_t u0, u1, v0, v1, d0, d1, q, r, tmp;
    fq_poly_init(u0, ctx);
    fq_poly_init(u1, ctx);
    fq_poly_init(v0, ctx);
    fq_poly_init(v1, ctx);
    fq_poly_init(d0, ctx);
    fq_poly_init(d1, ctx);
    fq_poly_init(q, ctx);
    fq_poly_init(r, ctx);
    fq_poly_init(tmp, ctx);

    /* set the values */
    fq_poly_zero(u1, ctx);
    fq_poly_zero(v0, ctx); 
    fq_poly_one(u0, ctx);
    fq_poly_one(v1, ctx);
    fq_poly_set(d0, a, ctx);
    fq_poly_set(d1, b, ctx);
    
    /* start the loop */
    while (fq_poly_degree(d0, ctx) > k && !fq_poly_is_zero(d1, ctx)) {
	fq_poly_divrem(q, r, d0, d1, ctx); /* quotient */

	fq_poly_set(d0, d1, ctx);
	fq_poly_set(d1, r, ctx);
    
	fq_poly_mul(tmp, u1, q, ctx); /* update of u factors */
	fq_poly_sub(tmp, u0, tmp, ctx);
	fq_poly_set(u0, u1, ctx);
	fq_poly_set(u1, tmp, ctx);
    
	fq_poly_mul(tmp, v1, q, ctx); /* update of v factors */
	fq_poly_sub(tmp, v0, tmp, ctx);
	fq_poly_set(v0, v1, ctx);
	fq_poly_set(v1, tmp, ctx);
    }

    /* put everything in the result variables */
    fq_poly_set(u, u0, ctx);
    fq_poly_set(v, v0, ctx);
    fq_poly_set(d, d0, ctx);

    /* clear everything */
    fq_poly_clear(tmp, ctx);
    fq_poly_clear(u0, ctx);
    fq_poly_clear(u1, ctx);
    fq_poly_clear(v0, ctx);
    fq_poly_clear(v1, ctx);
    fq_poly_clear(d0, ctx);
    fq_poly_clear(d1, ctx);
}



/* **************************************************************************** */
/*                                 MATRICES                                     */
/* **************************************************************************** */

/* expansion of a matrix over F_q^m to a matrix over FF_q
 * NEED TO OBTAIN OUTPUT AS A MATRIX IN FF_q
 */
void
fq_matrix_expand(fq_mat_t res, const fq_mat_t H, const fq_ctx_t ctx, const fq_ctx_t ctx_q) {
    int m = fq_ctx_degree(ctx);
    int nr, nc, tmp;
    fq_t tmp_q;
    fq_init(tmp_q, ctx_q);
    fmpz_poly_t pol;
    fmpz_poly_init(pol);
    
    /* nb of rows / columns of H */
    nr = fq_mat_nrows(H, ctx);
    nc = fq_mat_ncols(H, ctx);

    /* is there another way that this double loop ? */
    for (int j = 0; j < nc; j++ ) {
	for (int i = 0; i < nr; i++) {
	    fq_get_fmpz_poly(pol, fq_mat_entry(H,  i,  j), ctx);
	    for (int k = 0; k < m; k++) {
		tmp = fmpz_poly_get_coeff_ui(pol, k);
		fq_set_ui(tmp_q, tmp, ctx_q);
		fq_mat_entry_set(res, i*m + k, j, tmp_q, ctx_q);
	    }
	}
    }

    /* clearing memory */
    fq_clear(tmp_q, ctx_q);
    fmpz_poly_clear(pol);
}
