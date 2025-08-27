#include <iostream>
#include <random>

#include "Tools/tools_fq.hpp"
#include "Tools/tools.hpp"

using namespace std;


/* **************************************************************************** */
/*                              MISCELLANEOUS                                   */
/* **************************************************************************** */

/**
 * Hamming weight of a vector.
*/
int
hamming_weight(const fq_struct* v, const int len, const fq_ctx_t ctx) {
    int count = 0;
    for (int i = 0; i < len; ++i) {
	if (!fq_is_zero(&v[i], ctx)) {
	    count++;
	}
    }
    return count;
}


/**
 * distance between 2 vectors.
*/
int
hamming_distance(const fq_struct* v, const fq_struct* w, const int len, const fq_ctx_t ctx) {
    int count = 0;
    for (int i = 0; i < len; ++i) {
      if (!fq_equal(&v[i], &w[i], ctx)) {
	    count++;
	}
    }
    return count;
}


/* **************************************************************************** */
/*                              FINITE FIELDS                                   */
/* **************************************************************************** */

/**
 * Print vector
 */
void
_fq_vec_print_pretty(const fq_struct* v, const int len, const fq_ctx_t ctx) {
    for (int i = 0; i < len; ++i) {
	fq_print_pretty(&v[i], ctx);
	printf(" ");
    }
    printf("\n");
}


/**
 * Go from finite field vec to integral vector
 */
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


/**
 * Go from integral vector to finite field vec 
*/
void
_int_vec_2_fq(fq_struct* res, const int* a, const int len, const fq_ctx_t ctx) {
    _fq_vec_zero(res, len, ctx);
    for (int i = 0 ; i < len; ++i) {
	fq_set_ui(&res[i], a[i], ctx);
    }
}

/**
 * Go from integral vector (fmpz) to finite field vec 
 */
void
_fmpz_vec_2_fq(fq_struct* res, const fmpz* a, const int len, const fq_ctx_t ctx) {
    for (int i = 0 ; i < len; ++i) {
	fq_set_fmpz(&res[i], &a[i], ctx);
    }
}

/**
 * Go from finite field vector to  integral vector (fmpz) 
 */
void
_fq_vec_2_fmpz(fmpz* res, const fq_struct* a, const int len, const fq_ctx_t ctx) {
    int b;
    for (int i = 0 ; i < len; ++i) {
	b = fq_get_fmpz(&res[i], &a[i], ctx);
    }
}


/**
 * Computes the coeffs in F_q of an elt 'a' of F_q^m / F_q
 */
void fq_get_coeffs(fq_struct *res, const fq_struct a, const int m,
                   const fq_ctx_t ctx, const fq_ctx_t ctx_q) {

    fmpz_poly_t tmp_poly; fmpz_poly_init(tmp_poly);
    fmpz* tmp_vec = _fmpz_vec_init(m);

    fq_get_fmpz_poly(tmp_poly, &a, ctx);
    fmpz_poly_get_coeffs(tmp_vec, tmp_poly, m);
    _fmpz_vec_2_fq(res, tmp_vec, m, ctx_q);

    _fmpz_vec_clear(tmp_vec, m);
    fmpz_poly_clear(tmp_poly);
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
fq_vec_rand(fq_struct* res, const int len, const fq_ctx_t ctx,
		     flint_rand_t state) {
    for (int i =  0; i < len; ++i) {
	fq_randtest(&res[i], state, ctx);
    }
}

/* generate a random vector whose coeffs are in a finite field and all distincts */
void
fq_vec_rand_distinct(fq_struct* res, const int len, const fq_ctx_t ctx,
		     flint_rand_t state) {
    _fq_vec_randtest(res, state, len, ctx); /* NOT GOOD AT ALL → change it */
    while (fq_check_repeat(res, len, ctx)) {
	_fq_vec_randtest(res, state, len, ctx);
    }
}


/**
 * Generates a random vector whose coeffs are in a finite field and all distincts
 * WORK ONLY for *BINARY* fields
 * NOT THE MOST ELEGANT WAY -- MAYBE CHANGE
 */
void
fq_vec_rand_distinct_2(fq_struct* res, const int len, const fq_ctx_t ctx,
		       flint_rand_t state) {
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

void
fq_vec_shorten(fq_struct *res, const fq_struct *v, const int len,
                    const fq_ctx_t ctx) {
    _fq_vec_set(res, v, len, ctx);
}

/**
 * Expand vector from length len1 to length len2 with zero coefficients.
 */
void
fq_vec_expand(fq_struct *res, const fq_struct *v, const int len1, const int len2,
	      const fq_ctx_t ctx) {
    _fq_vec_set(res, v, len1, ctx);
    for (int i = len1; i < len2; ++i) {
	fq_zero(&res[i], ctx);
    }
}

/* **************************************************************************** */
/*                           POLYNOMIAL MANIPUTATION                            */
/* **************************************************************************** */

/**
 * Set the coeffs of polynomial
 */
void
fq_poly_set_coeffs(fq_poly_t f, const fq_struct* alpha, const int len,
		   const fq_ctx_t ctx) {
    fq_poly_zero(f, ctx);
    for (int i = 0; i < len; i++) {
	fq_poly_set_coeff(f, i, &alpha[i], ctx);
    }
}



/**
 * Get the coefficients of a polynomial over finite field and store it in a
 * vector.
 */
void
fq_poly_get_coeffs(fq_struct* res, const fq_poly_t f, const int len, const fq_ctx_t ctx) {
    for (int i = 0; i < len; i++) {
	fq_poly_get_coeff(&res[i], f, i, ctx);
    }
}

/**
 * Set the coefficients of a polynomial over integers (fmpz).
 */
void
fmpz_poly_set_coeffs(fmpz_poly_t f, const fmpz* a, const int len) {
    for (int i = 0; i < len; i++) {
      fmpz_poly_set_coeff_fmpz(f, i, &a[i]);
    }
}


/**
 * Get the coefficients of a polynomial over integers (fmpz) and store it in a
 * vector.
 */
void
fmpz_poly_get_coeffs(fmpz* res, const fmpz_poly_t f, const int len) {
    for (int i = 0; i < len; i++) {
	fmpz_poly_get_coeff_fmpz(&res[i], f, i);
    }
}


/**
 * Compute the positions of non-zero coefficients.
 * Assume pos has the good length.
 */
void
fq_poly_nonzero_coeffs(int* pos, const fq_poly_t f, const int r,
		       const fq_ctx_t ctx) {
    int j = 0;
    fq_t tmp; fq_init(tmp, ctx);
    for (int i = 0; i < r; ++i) {
	fq_poly_get_coeff(tmp, f, i, ctx);
	if (fq_is_one(tmp, ctx)) {
            pos[j] = i;
	    j++;
        }
    }
}


/**
 * Sets res to "cyclic polynomial generator", i.e. X^d - 1 (over FF_2)
 */
void
fq_poly_set_cyclic(fq_poly_t res, const int d, const fq_ctx_t ctx) {
    fq_poly_zero(res, ctx);
    fq_t tmp; fq_init(tmp, ctx);
    fq_one(tmp, ctx);
    fq_poly_set_coeff(res, d, tmp, ctx);
    fq_poly_set_coeff(res, 0, tmp, ctx);
    fq_clear(tmp, ctx);
}

/**
 * Sets res to "linear polynomial", i.e. X + alpha 
 */
void
fq_poly_set_linear(fq_poly_t res, const fq_t alpha, const fq_ctx_t ctx) {
    fq_poly_zero(res, ctx);
    fq_t tmp; fq_init(tmp, ctx);
    fq_one(tmp, ctx);
    fq_poly_set_coeff(res, 1, tmp, ctx);
    fq_poly_set_coeff(res, 0, alpha, ctx);
    fq_clear(tmp, ctx);
}

/**
 * Sets res to the product of "linear polynomials", i.e. prod_alpha X + alpha 
 */
void
fq_poly_set_linear_product(fq_poly_t res, const fq_struct* alpha, const int len,
			   const fq_ctx_t ctx) {
    fq_poly_one(res, ctx);

    fq_poly_t tmp_pol;
    fq_poly_zero(tmp_pol, ctx);


    fq_t tmp; fq_init(tmp, ctx);    
    fq_one(tmp, ctx);
    fq_poly_set_coeff(tmp_pol, 1, tmp, ctx);
    
    for (int i = 0; i < len; ++i) {
	fq_poly_set_coeff(tmp_pol, 0, &alpha[i], ctx);
	fq_poly_mul_classical(res, res, tmp_pol, ctx);
    }

    fq_poly_clear(tmp_pol, ctx);
    fq_clear(tmp, ctx);
}



/**
 * Check if f has a root in alpha.
 */
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


/**
 * Computes the unique polynomial P defined by points α and β,
 * i.e. s.t P(αᵢ) = βᵢ
 */
void
fq_poly_interpolate(fq_poly_t res, const fq_struct* alpha, const fq_struct*  beta,
		    const int len, const fq_ctx_t ctx) {
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


/**
 * Compute an irreducible polynomial of degree "deg" as in CM ie so that it
 * generates a Goppa code together with the roots alpha Γ(alpha, res)  
*/
void
cm_fq_poly_irr_pol(fq_poly_t& res, const int deg, const fq_struct* alpha, const int len,
		   const fq_ctx_t ctx, flint_rand_t state) {
    fq_poly_randtest_irreducible(res, state, deg+1, ctx);
    while (fq_poly_eval_zero(res, alpha, len, ctx)) {
	fq_poly_randtest_irreducible(res, state, deg+1, ctx);
    }
}



/**
 * Computes early abort extended gcd of a and b in finite field
 * defined by context ctx.
 */
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

/**
 * Expansion of a matrix over F_q^m to a matrix over FF_q
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




/**
 * Computes multiplication matrix of h in FF_q[X] / (P(X))
 */
void
fq_mult_matrix(fq_mat_t res, const fq_poly_t h, const fq_poly_t P, const fq_ctx_t ctx) {
    int i, j, d = fq_poly_degree(P, ctx);
    fq_t tmp; fq_init(tmp, ctx);
    fq_poly_t tmp_pol, gen;
    fq_poly_init(tmp_pol, ctx);     fq_poly_init(gen, ctx);
    fq_poly_set(tmp_pol, h, ctx);   fq_poly_gen(gen, ctx);
    
    for (i = 0; i < d; ++i) {
	// printf("it is the %dth index in mult mat comp. \t ", i);
	for (j = 0; j < d; ++j) {
	    fq_poly_get_coeff(tmp, tmp_pol, j, ctx);
	    fq_mat_entry_set(res, i, j, tmp, ctx);
	}
	fq_poly_mulmod(tmp_pol, tmp_pol, gen, P, ctx);
    }
    fq_clear(tmp, ctx);
    fq_poly_clear(tmp_pol, ctx);
    fq_poly_clear(gen, ctx);
}



/* **************************************************************************** */
/*                                   BIKE                                       */
/* **************************************************************************** */



/** 
 * Counter function as in Bike documentation.
 * Computes the number of agreeing bits in vector v and jth column vector of H.
 */
int
ctr(const fq_struct* v, const fq_mat_t& H,  const int j, const fq_ctx_t ctx) {
    int count = 0;
    
    int r = fq_mat_nrows(H, ctx);

    /* fq_equal(&v[i], fq_mat_entry(H, i, j), ctx)  */
    for (int i = 0; i < r; ++i) {
	if (fq_is_one(&v[i], ctx) && fq_is_one(fq_mat_entry(H, i, j), ctx)) {
	    count++;
	}
    }
    printf("count is: %d and r is %d \n", count, r);
    return count;
}


/** 
 * Counter function as in Bike documentation.
 * Computes the number of agreeing bits in vector v and jth column vector of H.
 */
int
ctrv2(const fq_poly_t sp,  const int* pos, const int j, const int weight, const int r,
      const fq_ctx_t ctx) {

    int count = 0;
    int k = 0;
    fq_t tmp; fq_init(tmp, ctx);
    
    for (int i = 0; i < weight; ++i) {
	k = (j + pos[i]) % r ;
	/* printf("%d ", k); */
	fq_poly_get_coeff(tmp, sp, k, ctx);
	if (fq_is_one(tmp, ctx)) {
	    count++;
	}
    }
    /* printf("\n"); */
    fq_clear(tmp, ctx);
    return count;
}


/**
 * Computes black and gray positions as in Bike specification doc.
 */
void
BFIter(fq_struct* e, int* black, int* gray, const fq_struct* s,
       const fq_mat_t& H, const int T, const int tau, const fq_ctx_t ctx) {
    int i, j, count;

    int n = fq_mat_ncols(H, ctx);

    fq_t one;  fq_init(one, ctx); fq_set_ui(one, 1, ctx);
    
    /* initialisation black and gray vectors */
    for (i = 0; i < n; ++i) {
	black[i] = 0;
	gray[i] = 0;
    }

    for (j = 0; j < n; ++j) {
	count = ctr(s, H, j, ctx);
	if (count >= T) {
	    fq_add(&e[j], &e[j], one, ctx);
	    black[j] = 1;
	} else if (count >= T - tau) {
	    gray[j] = 1;
	}
    }
    fq_clear(one, ctx);
}


/**
 * Computes black and gray positions as in Bike specification doc.
 */
void
BFIterv2(fq_poly_t e0, fq_poly_t e1, int* black, int* gray, const fq_poly_t sp,
	 const int* pos0, const int* pos1, const int weight, const int r, const int T,
	 const int tau, const fq_ctx_t ctx) {
    
    int i, j, count;

    int n = 2*r;
    
    fq_t one;  fq_init(one, ctx); fq_set_ui(one, 1, ctx);
    fq_t tmp; fq_init(tmp, ctx);
    
    /* initialisation black and gray vectors */
    for (i = 0; i < n; ++i) {
	black[i] = 0;
	gray[i] = 0;
    }
    
    for (j = 0; j < n; ++j) {
        if (j < r) {
            count = ctrv2(sp, pos0, j, weight, r, ctx);
	} else {
	    count = ctrv2(sp, pos1, j - r, weight, r, ctx);
	}
	/* printf("count is: %d and threshold is %d \n", count, T); */
        if (count >= T) {
	    // printf("count is %d and threshold is %d \n", count, T);
	    if (j < r) {
		fq_poly_get_coeff(tmp, e0, j, ctx);
		fq_add(tmp, tmp, one, ctx);
		fq_poly_set_coeff(e0, j, tmp, ctx);
	    } else {
		fq_poly_get_coeff(tmp, e1, j - r, ctx);
		fq_add(tmp, tmp, one, ctx);
		fq_poly_set_coeff(e1, j-r, tmp, ctx);
	    }
	    black[j] = 1;
	} else if (count >= T - tau) {
	    gray[j] = 1;
	}
    }
    fq_clear(one, ctx);
    fq_clear(tmp, ctx);
}


/**
 * Modify e using black or gray positions as in Bike specification doc.
 */
void
BFMaskedIter(fq_struct* e, const fq_struct* s, const fq_mat_t& H, const int T,
	     const int* mask, const fq_ctx_t ctx) {
    int j;

    int n = fq_mat_ncols(H, ctx);    

    fq_t tmp;
    fq_init(tmp, ctx);
    
    for (j = 0; j < n; ++j) {
	if (ctr(s, H, j, ctx) >= T) {
	    fq_set_ui(tmp, mask[j], ctx);
	    fq_add(&e[j], &e[j], tmp, ctx);
	} 
    }
    fq_clear(tmp, ctx);
}


/**
 * Modify e using black or gray positions as in Bike specification doc.
 */
void
BFMaskedIterv2(fq_poly_t e0, fq_poly_t e1, const fq_poly_t sp, const int* pos0,
	       const int* pos1, const int weight, const int r, const int T,
	       const int* mask, const fq_ctx_t ctx) {
    int j, count;

    int n = 2*r;    

    fq_t tmp, tmp1;
    fq_init(tmp, ctx); fq_init(tmp1, ctx);
    
    for (j = 0; j < n; ++j) {
	if (j < r) {
            count = ctrv2(sp, pos0, j, weight, r, ctx);
        } else {
	    count = ctrv2(sp, pos1, j - r, weight, r, ctx);
        }
	
	if (count >= T) {
	    fq_set_ui(tmp, mask[j], ctx);

	    if (j < r) {
		fq_poly_get_coeff(tmp1, e0, j, ctx);
		fq_add(tmp1, tmp1, tmp, ctx);
		fq_poly_set_coeff(e0, j, tmp1, ctx);
	    } else {
		fq_poly_get_coeff(tmp1, e1, j - r, ctx);
		fq_add(tmp1, tmp1, tmp, ctx);
		fq_poly_set_coeff(e1, j-r, tmp1, ctx);
	    }

	} 
    }
    fq_clear(tmp, ctx);
    fq_clear(tmp1, ctx);
}


/* ***************************************************************************** */
/*              Flint <-> GF2X Conversion functions (for testing)                */
/* ***************************************************************************** */

// Converts Flint polynom to uint64 vector for GF2X use
void fq_poly_to_uint64(const fq_poly_t poly, const fq_ctx_t ctx, uint64_t* result, slong* num_blocks_out) {
    slong degree = fq_poly_degree(poly, ctx); // Polynom degree
    slong num_blocks = (degree / 64) + 1;     // Number of 64 bits blocks needed

    // Initialisation
    for (slong i = 0; i < num_blocks; ++i) {
        result[i] = 0;
    }

    // Conversion
    for (slong i = 0; i <= degree; ++i) {
        fq_t x; fq_init(x, ctx);
        fq_poly_get_coeff(x, poly, i, ctx);
        if (!fq_is_zero(x, ctx)) {
            result[i / 64] |= (1 << (i % 64)); 
        }
    }
    *num_blocks_out = num_blocks; // Useless
}

// FLINT F_2 polynomial to gf2x
GF2X fq_poly_to_gf2x(const fq_poly_t in_poly, const fq_ctx_t ctx){
    slong degree = fq_poly_degree(in_poly, ctx);      // Polynom degree
    GF2X out_poly;
    out_poly.SetMaxLength(degree+1);
    out_poly.SetLength(degree+1);
    // Conversion
    fq_t x; fq_init(x, ctx);        // Temporary variable 
    //printf("Degree : %i\n", degree);
    for (slong i = 0; i <= degree; ++i) {
        fq_poly_get_coeff(x, in_poly, i, ctx);
        if (!fq_is_zero(x, ctx)) {
          SetCoeff(out_poly, i); 
        }
    }
    out_poly.normalize(); // Maybe not necessary
    fq_clear(x,ctx);
    return out_poly;
}

// FLINT F_2^m polynomial to gf2ex 
GF2EX fq_poly_to_gf2ex(const fq_poly_t in_poly, const fq_ctx_t ctx){ //, GF2EX* out_poly){
    slong degree = fq_poly_degree(in_poly, ctx);      // Polynom degree
    slong m = fq_ctx_degree(ctx);                     // Extension degree
    
    GF2EX out_poly;
    out_poly.SetMaxLength(degree+1); // Memory allocation (Mandatory)
    out_poly.SetLength(degree+1);    // Poolynom degree (current)

    // Conversion
    fq_t x; fq_init(x, ctx);         // Temporary variable 

    for (slong i = 0; i <= degree; ++i) {
        GF2X poly_coeff;
        fmpz_poly_t tmp_fmpz;
        fmpz_poly_init(tmp_fmpz);

        fq_poly_get_coeff(x, in_poly, i, ctx);
        fq_get_fmpz_poly(tmp_fmpz,x,ctx);

        for(int j = 0; j < m ; j++){
          if (fmpz_poly_get_coeff_ui(tmp_fmpz,j)) {
            SetCoeff(poly_coeff, j); 
          }
        }
        GF2E alpha = to_GF2E(poly_coeff); 

        out_poly[i] = alpha;
        fmpz_poly_clear(tmp_fmpz);
    }

    fq_clear(x,ctx);
    return out_poly;
}

// Covnerts each coeff of the GF2EX polynomial into an array
void GF2EX_to_fq_poly(const GF2EX& in_pol, fq_poly_t f, int len,  const fq_ctx_t ctx){
    vec_GF2E alpha = VectorCopy(in_pol,len);
    fq_struct* tmp = _fq_vec_init(len, ctx);
    vec_GF2E_to_fq_vec(tmp,len,alpha,ctx);
    fq_poly_set_coeffs(f,tmp,len,ctx);
}


// FLINT F_2 Matrix -> GF2 Matrix (GF2! not GF2X, 0-1 only)
mat_GF2 fq_mat_to_gf2_mat(const fq_mat_t& T, const fq_ctx_t ctx_q){
    int n = fq_mat_nrows(T, ctx_q);
    int m = fq_mat_ncols(T, ctx_q);

    mat_GF2 M;
    M.SetDims(n, m); 

    for(int row = 0; row < n ; row ++){
        for(int col = 0; col < m ; col++){
        M[row][col] = fq_is_one(fq_mat_entry(T, row, col), ctx_q);
        }
    }

    return M;
}

// FLINT F_2^m Matrix -> GF2E Matrix
mat_GF2E fq_mat_to_gf2e_mat(const fq_mat_t& T, const fq_ctx_t& ctx){
    int row = fq_mat_nrows(T, ctx);
    int col = fq_mat_ncols(T, ctx);
    long m = fq_ctx_degree(ctx); // Extension degree

    mat_GF2E M;
    M.SetDims(row, col); 

    fmpz_poly_t poly_repr;
    fmpz_poly_init(poly_repr);

    for (slong i = 0; i < row; i++) {
        for(slong j = 0; j<col; j++){
            GF2X poly_coeff;
            fq_get_fmpz_poly(poly_repr, fq_mat_entry(T,i,j), ctx);  

            for (long k = 0; k < m; k++) {
                if (fmpz_poly_get_coeff_ui(poly_repr, k)) {
                    SetCoeff(poly_coeff,k); 
                }
            }
            GF2E tmp_GF2E = to_GF2E(poly_coeff);
            M[i][j]= tmp_GF2E;

        }
    }
    fmpz_poly_clear(poly_repr);

    return M;
}

// FLINT F_2 vector to GF2 vector
vec_GF2 fq_vec_to_vec_GF2(const fq_struct* fq_vec, slong len, const fq_ctx_t ctx) {
    long m = fq_ctx_degree(ctx); // Extension degree
    
    vec_GF2 result;
    result.SetLength(len * m); 

    fmpz_poly_t poly_repr;
    fmpz_poly_init(poly_repr);

    for (slong i = 0; i < len; i++) {
        fq_get_fmpz_poly(poly_repr, fq_vec + i, ctx);

        for (long j = 0; j < m; j++) {
            if (fmpz_poly_get_coeff_ui(poly_repr, j)) {
                result[i * m + j] = 1; 
            }
        }
    }

    fmpz_poly_clear(poly_repr);
    return result;
}


// Flint F_2^m vector to GF2E vector
vec_GF2E fq_vec_to_vec_GF2E(const fq_struct* fq_vec, slong len, const fq_ctx_t ctx){
    long m = fq_ctx_degree(ctx); // Extension degree
    vec_GF2E result;
    result.SetLength(len);  

    fmpz_poly_t poly_repr;
    fmpz_poly_init(poly_repr);

    for (slong i = 0; i < len; i++) {
    
        GF2X poly_coeff;
        fq_get_fmpz_poly(poly_repr, fq_vec + i, ctx); 

        // Store each coefficient
        for (long j = 0; j < m; j++) {
            if (fmpz_poly_get_coeff_ui(poly_repr, j)) {
                SetCoeff(poly_coeff,j); 
            }
        }
        poly_coeff.normalize();
        GF2E tmp_GF2E = to_GF2E(poly_coeff);
        result[i]= tmp_GF2E;
    }

    fmpz_poly_clear(poly_repr);
    return result;

}

// vec_GF2E vector to Flint F_2^m vector 
void vec_GF2E_to_fq_vec(fq_struct* fq_vec, slong len, const vec_GF2E& input_vec, const fq_ctx_t ctx) {
    long m = fq_ctx_degree(ctx); // Extension degree

    fmpz_poly_t poly_repr;
    fmpz_poly_init(poly_repr);

    for (slong i = 0; i < len; i++) {
        GF2X poly_coeff = rep(input_vec[i]); // GF2X Représentation of the GF2E pol.

        fmpz_poly_zero(poly_repr); 

        for (long j = 0; j < m; j++) {
            if (IsOne(coeff(poly_coeff, j))) { 
                fmpz_poly_set_coeff_ui(poly_repr, j, 1); 
            }
        }

        fq_set_fmpz_poly(&fq_vec[i], poly_repr, ctx);
 
    }
    fmpz_poly_clear(poly_repr);

}
