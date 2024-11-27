#include "codes.hpp"
#include <iostream>
#include <flint/fq_mat.h>


/* **************************************************************************** */
/*                                   DECODERS                                   */
/* **************************************************************************** */



/**
 * Decoding algorithm for Reed-Solomon codes
 * Decode received message c wrt to RS_k(alpha)
*/
int
RS_decoder(fq_poly_t res, const fq_struct* c, const fq_struct* alpha, const int n, const int k,
	   const fq_ctx_t ctx) {

    /* declaration and init of local variables */
    fq_poly_t A, P, tmp, u, v, D;
    fq_poly_init(A, ctx); fq_poly_one(A, ctx);
    fq_poly_init(tmp, ctx); fq_poly_gen(tmp, ctx);
    fq_poly_init(P, ctx);
    fq_poly_init(u, ctx);		/* u, v, D for xgcd */
    fq_poly_init(v, ctx);
    fq_poly_init(D, ctx);

    fq_t a; fq_init(a, ctx);

    for (int i = 0; i < n; ++i) {
	fq_neg(a, &alpha[i], ctx);
	fq_poly_set_coeff(tmp, 0, a, ctx);
	fq_poly_mul(A, A, tmp, ctx);
    }
    fq_poly_interpolate(P, alpha, c, n, ctx);

    int d = k + (n-k +1 ) / 2 - 1; /* computes k + ceil((n-k)/2) -1 */
    xgcd_abort(u, v, D, A, P, d, ctx);
    int b = fq_poly_divides(res, D, v, ctx);
    b = b && (fq_poly_degree(res, ctx) <= k);
    
    /* clear memory */
    fq_clear(a, ctx);
    fq_poly_clear(A, ctx);
    fq_poly_clear(P, ctx);
    fq_poly_clear(tmp, ctx);
    fq_poly_clear(u, ctx);
    fq_poly_clear(v, ctx);
    fq_poly_clear(D, ctx);
  
    return b;
}


/**
 * Decoding algorithm for Generalised Reed-Solomon codes
 */
int
GRS_decoder(fq_struct* res, const fq_struct* c, const fq_struct* alpha, const fq_struct* beta,
	    const int n, const int k, const fq_ctx_t ctx) { 
    /* declaration and init of local variables */
    fq_t tmp; fq_init(tmp, ctx);
    fq_struct* c1 = _fq_vec_init(n, ctx);
    fq_poly_t p;   fq_poly_init(p, ctx);

    /* goes from GRS to RS rep. */
    for (int i = 0; i < n; i++) {
	fq_inv(tmp, &beta[i], ctx);
	fq_mul(&c1[i], &c[i], tmp, ctx);
    }

    /* use the RS decoder */
    int b = RS_decoder(p, c1, alpha, n, k, ctx);
    if (b != 0) {
	for(int ii =0; ii < n; ii++) {
	    fq_poly_evaluate_fq(&res[ii], p, &alpha[ii], ctx);
	    fq_mul(&res[ii], &res[ii], &beta[ii], ctx);
	}
    }
    
    /* clear memory */
    fq_clear(tmp, ctx);
    fq_poly_clear(p, ctx);
    _fq_vec_clear(c1, n, ctx);
  
    return b;
}


/* Decoding algorithm for Goppa codes */
int
Goppa_decoder(fq_struct* res, const fq_struct* c, const fq_struct* alpha, const fq_poly_t g,
	      const int n, const int k, const fq_ctx_t ctx) {

    /* declaration and init of local variables */
    fq_poly_t A, tmp_poly, Aprime;
    fq_poly_init(A, ctx); fq_poly_one(A, ctx);
    fq_poly_init(Aprime, ctx);
    fq_poly_init(tmp_poly, ctx); fq_poly_gen(tmp_poly, ctx);
    fq_struct* beta = _fq_vec_init(n, ctx);
    fq_t tmp; fq_init(tmp, ctx);
  

    for (int i = 0; i < n; ++i) {
	fq_neg(tmp, &alpha[i], ctx);
	fq_poly_set_coeff(tmp_poly, 0, tmp, ctx);
	fq_poly_mul(A, A, tmp_poly, ctx);
    }
    fq_poly_derivative(Aprime, A, ctx);
  
    for (int i = 0; i < n; i++) {
	fq_poly_evaluate_fq(&beta[i], g, &alpha[i], ctx);
	fq_poly_evaluate_fq(tmp, Aprime, &alpha[i], ctx);
	fq_div(&beta[i], &beta[i], tmp, ctx);
    }

    int b = GRS_decoder(res, c, alpha, beta, n, k, ctx);

    /* clearing memory */
    fq_poly_clear(A, ctx);
    fq_poly_clear(Aprime, ctx);
    fq_poly_clear(tmp_poly, ctx);
    _fq_vec_clear(beta, n, ctx);
    fq_clear(tmp, ctx);
  
    return b;
}


/* Decoding algorithm for *BINARY* Goppa codes */
int
Goppa_decoder_bin(fq_struct* res, const fq_struct* c, const fq_struct* alpha,
		  const fq_poly_t g, const int n, const int k, const fq_ctx_t ctx) {

    /* declaration and init of local variables */
    fq_poly_t A, tmp_poly, Aprime;
    fq_poly_init(A, ctx); fq_poly_one(A, ctx);
    fq_poly_init(Aprime, ctx);
    fq_poly_init(tmp_poly, ctx); fq_poly_gen(tmp_poly, ctx);
    fq_struct* beta = _fq_vec_init(n, ctx);
    fq_t tmp; fq_init(tmp, ctx);
 

    for (int i = 0; i < n; ++i) {
	fq_neg(tmp, &alpha[i], ctx);
	fq_poly_set_coeff(tmp_poly, 0, tmp, ctx);
	fq_poly_mul(A, A, tmp_poly, ctx);
    }
    fq_poly_derivative(Aprime, A, ctx);
  
    for (int i = 0; i < n; i++) {
	fq_poly_evaluate_fq(&beta[i], g, &alpha[i], ctx);
	fq_sqr(&beta[i], &beta[i], ctx);
	fq_poly_evaluate_fq(tmp, Aprime, &alpha[i], ctx);
	fq_div(&beta[i], &beta[i], tmp, ctx);
    }
    int b = GRS_decoder(res, c, alpha, beta, n, k, ctx);
    return b;
}




/* **************************************************************************** */
/*                                   MISC                                       */
/* **************************************************************************** */

/* random codeword of the Goppa code Gamma(alpha, g) */
void
Goppa_codeword_random(fq_struct* res, const fq_struct* alpha, const fq_poly_t g, const int len,
     		      const fq_ctx_t ctx, flint_rand_t state) {
    fq_poly_t A, Aprime, poltmp;
    fq_struct* beta = _fq_vec_init(len, ctx);
    int i, pr;
    fq_t tmp;

    fq_init(tmp, ctx);
    fq_one(tmp, ctx);
  
    fq_poly_init(A, ctx);
    fq_poly_init(Aprime, ctx);
    fq_poly_init(poltmp, ctx);

    _fq_poly_set_length(poltmp, 1, ctx);
    fq_poly_set_coeff(poltmp, 1, tmp, ctx);
    fq_poly_one(A, ctx);

    fq_t a; fq_init(a, ctx);
    for (i = 0; i < len; i++) {
	fq_neg(a, &alpha[i], ctx);
	fq_poly_set_coeff(poltmp, 0, a, ctx);
	fq_poly_mul(A, A, poltmp, ctx);
    }
  
    fq_poly_derivative(Aprime, A, ctx);  
  
    for (i = 0; i < len; i++) {
	fq_poly_evaluate_fq(&beta[i], g, &alpha[i], ctx);
	fq_poly_evaluate_fq(tmp, Aprime, &alpha[i], ctx);
	fq_div(&beta[i], &beta[i], tmp, ctx);
    }
  
  
    fq_poly_randtest(poltmp, state, len-fq_poly_degree(g, ctx)-1, ctx);

    for (i = 0; i < len; i++) {
	fq_poly_evaluate_fq(tmp, poltmp, &alpha[i], ctx);
	fq_mul(tmp, tmp, &beta[i], ctx);
	fq_set(&res[i], tmp, ctx);
    }
  

    /* clearing memory */
    fq_poly_clear(A, ctx);
    fq_poly_clear(Aprime, ctx);
    fq_poly_clear(poltmp, ctx);
    _fq_vec_clear(beta, len, ctx);
    fq_clear(a, ctx);
    fq_clear(tmp, ctx);
}


void
Goppa_codeword_random_bin(fq_struct* res, const fq_struct* alpha, const fq_poly_t g,
			  const int len, const fq_ctx_t ctx, flint_rand_t state) {
    fq_poly_t g2;
    fq_poly_init(g2, ctx);
    fq_poly_sqr(g2, g, ctx);
  
    Goppa_codeword_random(res, alpha, g2, len, ctx, state);
  
    fq_poly_clear(g2, ctx);
}


/* Computes the parity check matrix of a Goppa code */
void
Goppa_parity_check(fq_mat_t res, const fq_struct* alpha, const fq_poly_t g, const fq_ctx_t ctx) {
    int i, j, n, m;
    n = fq_mat_nrows(res, ctx);
    m = fq_mat_ncols(res, ctx);
    fq_t tmp;   fq_init(tmp, ctx);
    fq_struct* beta = _fq_vec_init(m, ctx);

    /* put 1/g(alpha_i) in a vector */
    for (j = 0; j < m; j++) {
	fq_poly_evaluate_fq(tmp, g, &alpha[j], ctx);
	fq_inv(&beta[j], tmp, ctx);
    }
  
    /* computes the first line */
    for (j = 0; j < m; j++) {
	fq_mat_entry_set(res, 0, j, &beta[j], ctx);
    }
    
    /* computes the following lines by multiplying by alpha_i */
    for (i = 1; i < n; i++) {
	for (j = 0; j < m; j++) {
	    fq_mul(tmp, fq_mat_entry(res, i-1, j), &alpha[j], ctx);
	    fq_mat_entry_set(res, i, j, tmp, ctx);
	}
    }

    /* clearing memory */
    fq_clear(tmp, ctx);
    _fq_vec_clear(beta, m, ctx);
}


/**
 * Computes the parity check matrix of a Goppa code knowing it is *BINARY*
 */
void
Goppa_parity_check_bin(fq_mat_t res, const fq_struct* alpha, const fq_poly_t g, const fq_ctx_t ctx) {
    int i, j, n, m;
    fq_poly_t g2;
    fq_poly_init(g2, ctx);
    fq_poly_sqr(g2, g, ctx);
  
    Goppa_parity_check(res, alpha, g2, ctx);

    fq_poly_clear(g2, ctx);
}


/* **************************************************************************** */
/*                                   ENCODERS                                   */
/* **************************************************************************** */

/**
 * Encoding for RS codes.
 */
void
RS_encoding(fq_struct* res, const fq_poly_t f, const fq_struct* alpha, const int n,
	    const fq_ctx_t ctx) {
    int i = 0;
    for (i = 0; i < n; ++i) {
	fq_poly_evaluate_fq(&res[i], f, &alpha[i], ctx);
    }
}


/**
 * Encoding for Reed-Muller codes of parameters (1, m)
 * IN: alpha is a vector over FF_2
 * OUT: vector of 1 + sum_{0 < i < m} alpha_i *  X_i evaluated in all elements of FF_2^n
 */
void
RM_encoding(fq_struct* res, const fq_struct* alpha, const int m, const fq_ctx_t ctx) {
    fq_t tmp; fq_init(tmp, ctx);
    fq_t tmp1; fq_init(tmp1, ctx);
    
    for (int i=0; i < (1 << m); ++i) {
	 fq_set(tmp, &alpha[0], ctx);
	 fq_zero(tmp1, ctx);
	 for (int j = 0; j < m; ++j) {
	     fq_set_ui(tmp1, (i>>j)&1, ctx);
	     fq_mul(tmp1, tmp1, &alpha[j+1], ctx);
	     fq_add(tmp, tmp, tmp1, ctx);
	 }
	 fq_set(&res[i], tmp, ctx);
    }  
}

/**
 * Encoding for duplicated Reed-Muller codes of parameters (1, m)
 * IN: alpha is a vector over FF_2
 * OUT: vector of 1 + sum_{0 < i < m} alpha_i *  X_i evaluated in all elements of FF_2^n
 * duplicated r times
 */
void
RM_encoding_duplicated(fq_struct* res, const fq_struct* alpha, const int m, const int r,
		       const fq_ctx_t ctx) {
    fq_struct* tmp_vec = _fq_vec_init(1<<m, ctx);
    RM_encoding(tmp_vec, alpha, m, ctx);
    for (int i = 0; i < r; ++i) {
	for (int j = 0; j < (1<<m); ++j) {
	    fq_set(&res[i*(1<<m) + j], &tmp_vec[j], ctx);
	}
    }
}


/**
 * Encodes with concatenated Reed-Solomon and Reed-Muller codes.
 * Here `f` is a polynomial in F_2^m[X], `alpha` is a vector of elements in F_2^m,
 * `n` is the length of the RS code, we are using RM(1, m-1) duplicated `r` times
 * and `ctx` / `ctx_q` is F_2^m / F_2.
 */
void
RS_RM_concatenated_encoding(fq_struct* res, const fq_poly_t f,
			    const fq_struct* alpha, const int n, const fq_ctx_t ctx,
			    const int r, const fq_ctx_t ctx_q) {

    int i, j;
    int m = fq_ctx_degree(ctx);
    int len = (1 << (m-1));	/* length of RM(1, m-1) code */
    
    /* first use the RS code to encode the message f */
    fq_struct* tmp_m1 = _fq_vec_init(n, ctx);
    RS_encoding(tmp_m1, f, alpha, n, ctx);

    

    for (i = 0; i < n; ++i) {
	/* transform each coordinate of tmp_m1 / F_q^k into a vector in F_q^k */
	fq_struct* tmp_m2 = _fq_vec_init(m, ctx_q);
	fq_get_coeffs(tmp_m2, tmp_m1[i], m, ctx, ctx_q);

	/* encode said vector using a duplicated RM code */
	fq_struct* tmp_m3 = _fq_vec_init(r * len, ctx_q);
	RM_encoding_duplicated(tmp_m3, tmp_m2, m-1, r, ctx_q);
	
	for (j = 0; j < r * len; j++) {
	    fq_set(&res[i * (r*len) + j], &tmp_m3[j], ctx_q);
	}
    }

}



/**
 * Encoding as in CM.
 * Compute the syndrome associated s = H * e associated to e.
*/
void
CM_encoding(fq_struct* res, const fq_struct* e, const fq_mat_t& T, const int len, 
	    const fq_ctx_t& ctx_q)  {

    /* nb of rows and columns in T */
    int n = fq_mat_nrows(T, ctx_q);
    int m = fq_mat_ncols(T, ctx_q);

    /* create parity check matrix H = [Id | T] */   
    fq_mat_t I, H;
    fq_mat_init(I, n, n, ctx_q);
    fq_mat_one(I, ctx_q);
    fq_mat_init(H, n, len, ctx_q);
    fq_mat_concat_horizontal(H, I, T, ctx_q);

    fq_mat_clear(I, ctx_q); /* clearing memory */

    /* computes H * e */
    fq_mat_mul_vec(res, H, e, len, ctx_q);
    
    fq_mat_clear(H, ctx_q); /* clearing memory */
}



/**
 * Encoding as in Bike
 * Compute s = e0 + e1 * h modulo X^r - 1
*/
void 
Bike_encoding(fq_poly_t res, const int* e, const fq_poly_t h, const int r, const fq_ctx_t& ctx_q) {
    /* init P to X^r -1 */
    fq_poly_t P;
    fq_poly_init(P, ctx_q); 
    fq_poly_zero(P, ctx_q);
    fq_t tmp; fq_init(tmp, ctx_q);
    fq_one(tmp, ctx_q);
    fq_poly_set_coeff(P, r, tmp, ctx_q);
    fq_poly_set_coeff(P, 0, tmp, ctx_q);
    
    fq_poly_t e0, e1;
    fq_poly_init(e0, ctx_q);
    fq_poly_init(e1, ctx_q);
    
    fq_struct* tmp_vec = _fq_vec_init(r, ctx_q);

    int E[r];
    
    /* Set e0 to e[0..r-1] */
    for (int i = 0 ; i < r; ++i)
	E[i] = e[i];
    _int_vec_2_fq(tmp_vec, E, r, ctx_q);
    fq_poly_set_coeffs(e0, tmp_vec, r, ctx_q);

    /* set e1 to e[r..2r-1] */
    for (int i = 0 ; i < r; ++i)
	E[i] = e[i+r];
    _int_vec_2_fq(tmp_vec, E, r, ctx_q);
    fq_poly_set_coeffs(e1, tmp_vec, r, ctx_q);


    fq_poly_mulmod(res, e1, h, P, ctx_q);
    fq_poly_add(res, res, e0, ctx_q);


    fq_clear(tmp, ctx_q);
    _fq_vec_clear(tmp_vec, r, ctx_q);
    fq_poly_clear(e0, ctx_q);
    fq_poly_clear(e1, ctx_q);
    fq_poly_clear(P, ctx_q);
}


/**
 * Encoding as in HQC
 */
void
HQC_encoding( fq_poly_t res1, fq_poly_t res2, const fq_poly_t m, const fq_poly_t h,
	     const fq_poly_t s, const fq_struct* alpha, const int n, const int k,
	     const int n1, const int r, const int we, const int wr,
	      const fq_ctx_t ctx, const fq_ctx_t ctx_q, flint_rand_t state) {

    int _m = fq_ctx_degree(ctx);
    int len = (1 << (_m - 1));	/* length of RM(1, m-1) code */
    
    fq_struct* tmp_vec = _fq_vec_init(n, ctx_q);
    fq_struct* tmp_vec_2 = _fq_vec_init(r * n1 * len, ctx_q);
    fq_poly_t tmp_pol; fq_poly_init(tmp_pol, ctx_q);
    
    fq_poly_t r1; fq_poly_init(r1, ctx_q);
    fq_poly_t r2; fq_poly_init(r2, ctx_q);
    fq_poly_t P; fq_poly_init(P, ctx_q); 
    fq_poly_set_cyclic(P, n, ctx_q);
    
    int v[n];
    
    /* computes r1, r2 with fixed weight `wr`  */
    for (int i = 0 ; i < n; ++i) {
	v[i] = 0;
    }
    hqc_gen_e(v, n, wr);	/* r1 */
    
    _int_vec_2_fq(tmp_vec, v, n, ctx_q);
    fq_poly_set_coeffs(r1, tmp_vec, n, ctx_q);
    
    for (int i = 0 ; i < n; ++i) {
	v[i] = 0;
    }
    _fq_vec_zero(tmp_vec, n, ctx_q);
    hqc_gen_e(v, n, wr);	/* r2 */
    _int_vec_2_fq(tmp_vec, v, n, ctx_q);
    fq_poly_set_coeffs(r2, tmp_vec, n, ctx_q);
    
    /* put the result of r1 + h * r2 into res1 */
    fq_poly_mulmod(res1, h, r2, P, ctx_q);
    fq_poly_add(res1, res1, r1, ctx_q);

    fq_poly_mulmod(res2, s, r2, P, ctx_q);
    
    _fq_vec_zero(tmp_vec_2, r*len*n1,  ctx_q); fq_poly_zero(tmp_pol, ctx_q);
    RS_RM_concatenated_encoding(tmp_vec_2, m, alpha, n1, ctx, r, ctx_q);
    fq_poly_set_coeffs(tmp_pol, tmp_vec_2, r*len*n1, ctx_q);
   
    fq_poly_add(res2, res2, tmp_pol, ctx_q);

    /* computes e with fixed weight `we`  */
    for (int i = 0 ; i < n; ++i) {
	v[i] = 0;
    }
    _fq_vec_zero(tmp_vec, n, ctx_q); fq_poly_zero(tmp_pol, ctx_q);
    hqc_gen_e(v, n, we);
    _int_vec_2_fq(tmp_vec, v, n, ctx_q);
    fq_poly_set_coeffs(tmp_pol, tmp_vec, n, ctx_q);
    fq_poly_add(res2, res2, tmp_pol, ctx_q);
    
    fq_poly_truncate(res2, r*len*n1, ctx_q);
     
    
    /* clear memory */
    fq_poly_clear(tmp_pol, ctx_q);
    fq_poly_clear(P, ctx_q);
    fq_poly_clear(r1, ctx_q);
    fq_poly_clear(r2, ctx_q);
    _fq_vec_clear(tmp_vec, n, ctx_q);
}

/* **************************************************************************** */
/*                                   DECODERS                                   */
/* **************************************************************************** */



/* decode syndrome as in CM */
int
CM_syndrome_decoding(fq_struct* res, const fq_struct* s, const fq_struct* alpha,
		     const fq_poly_t g, const int len, const int t, const fq_ctx_t ctx) {
    int d = fq_poly_degree(g, ctx);
    int m = fq_ctx_degree(ctx);

    /* first needs to expand syndrome s */
    fq_struct* ss = _fq_vec_init(len, ctx);
    _fq_vec_zero(ss, len, ctx);
    for (int i = 0; i < d * m; i++) {
	if (!fq_is_zero(&s[i], ctx)) fq_set_ui(&ss[i], 1, ctx);
    }

    /* then we use the Goppa decoder */
    int b = Goppa_decoder(res, ss, alpha, g, len, len-2*t,  ctx);

    /* clearing memory */
    _fq_vec_clear(ss, len, ctx);
  
    return b;
}


/* decode syndrome as in CM */
int
CM_syndrome_decoding_bin(fq_struct* res, const fq_struct* s, const fq_struct* alpha,
			 const fq_poly_t g, const int len, const int t, const fq_ctx_t ctx) {
    int d = fq_poly_degree(g, ctx);
    int m = fq_ctx_degree(ctx);

    /* first needs to expand syndrome s */
    fq_struct* ss = _fq_vec_init(len, ctx);
    _fq_vec_zero(ss, len, ctx);
    for (int i = 0; i < d * m; i++) {
	if (!fq_is_zero(&s[i], ctx)) fq_set_ui(&ss[i], 1, ctx);
    }
 
    /* then we use the Goppa decoder */
    int b = Goppa_decoder_bin(res, ss, alpha, g, len, len-2*t, ctx);

    /* clearing memory */
    _fq_vec_clear(ss, len, ctx);
  
    return b;
}



/**
 * Decoding as in Bike : use the Black-Gray-Flip (BGF) algorithm
*/
int
Bike_decoding(fq_struct* res, const fq_struct* s, const fq_mat_t& H,
	      const int weight, const int NbIter, const int tau,
	      const fq_ctx_t ctx_q) {
  
    int i, T, w;
    int r = fq_mat_nrows(H, ctx_q);
    int n = 2*r;
    int d = (weight/2 + 1)/2 + 1;
    int b = 0;
    fq_struct* e = _fq_vec_init(n, ctx_q);
    fq_struct* tmp = _fq_vec_init(r, ctx_q);
    fq_struct* synd = _fq_vec_init(r, ctx_q);
    int black[n]; int gray[n];

    
    for (i = 0; i < NbIter; ++i) {
	printf("This is the %dth iteration !! \n", i);
	fq_mat_mul_vec(tmp, H, e, n, ctx_q); /* CHANGE THIS : need polynomial
						manipulation (?) */
	_fq_vec_add(tmp, tmp, s, r, ctx_q);
	w = hamming_weight(tmp, r, ctx_q);
	T = compute_threshold(w, i, r);
	printf("The threshold is %d\n", T);
	
	BFIter(e, black, gray, tmp , H, T, tau, ctx_q);
	
	if (i==0) {
	    BFMaskedIter(e, tmp, H, d, black, ctx_q); /* CHANGE THIS : change
							 dependency on H */
	    BFMaskedIter(e, tmp, H, d, gray, ctx_q);
	}
	

    }
    
    fq_mat_mul_vec(synd, H, e, n, ctx_q); /* computation of syndrome : remove H */

    if (_fq_vec_equal(s, synd, r, ctx_q)) {
	_fq_vec_set(res, e, n, ctx_q);
	b = 1;
    }
    _fq_vec_clear(synd, r, ctx_q);
    _fq_vec_clear(tmp, r, ctx_q);
    _fq_vec_clear(e, n, ctx_q);

    return b;
}


/**
 * Decoding as in Bike : use the Black-Gray-Flip (BGF) algorithm.
 * Comment : we use a polynomial representation as we should.
*/
int
Bike_decoding_v2(fq_struct* res, const fq_struct* s, const fq_poly_t h0,
		 const fq_poly_t h1, const int r, const int weight,
		 const int NbIter, const int tau,
		 const fq_ctx_t ctx_q) {
    int i, T, w;
    int n = 2*r;
    int d = (weight + 1)/2 + 1;
    int b = 0;

    /* P = X^r - 1 */
    fq_poly_t P; fq_poly_init(P, ctx_q); fq_poly_set_cyclic(P, r, ctx_q);
    
    /* position of non-zero coefficients of h0, h1 (secret key) */
    int pos0[weight], pos1[weight];
    fq_poly_nonzero_coeffs(pos0, h0, r, ctx_q);
    fq_poly_nonzero_coeffs(pos1, h1, r, ctx_q);

    /* polynomials used in computation */
    fq_poly_t e0; fq_poly_init(e0, ctx_q); fq_poly_zero(e0, ctx_q);
    fq_poly_t e1; fq_poly_init(e1, ctx_q); fq_poly_zero(e1, ctx_q);
    fq_poly_t tmp; fq_poly_init(tmp, ctx_q); fq_poly_zero(tmp, ctx_q);
    fq_poly_t tmp2; fq_poly_init(tmp2, ctx_q); fq_poly_zero(tmp2, ctx_q);
    fq_poly_t sp;  fq_poly_init(sp, ctx_q); fq_poly_zero(sp, ctx_q);
    fq_poly_set_coeffs(sp, s, r, ctx_q);

    
    fq_poly_t synd; fq_poly_init(synd, ctx_q); fq_poly_zero(synd, ctx_q);
    
    int black[n]; int gray[n];
    
    for (i = 0; i < NbIter; ++i) {
       	
	/* Computes s' = s + e*h */
	fq_poly_mulmod(tmp, e0, h0, P, ctx_q);
	fq_poly_mulmod(tmp2, e1, h1, P, ctx_q);
	fq_poly_add(tmp, tmp, tmp2, ctx_q);
	fq_poly_add(tmp, tmp, sp, ctx_q);
	
	/* Hamming weight of s' + threshold */
	w = fq_poly_hamming_weight(tmp, ctx_q);
	T = compute_threshold(w, i, r);
	
	/* BFIter with polynomials */
	BFIterv2(e0, e1, black, gray, tmp, pos0, pos1, weight, r, T, tau, ctx_q);
	
	
	if (i==0) {
	    BFMaskedIterv2(e0, e1, tmp, pos0, pos1, weight, r, d, black, ctx_q); 
	    BFMaskedIterv2(e0, e1, tmp, pos0, pos1, weight, r, d, gray, ctx_q);
	}
	
    }

    /* Computes the syndrome */
    fq_poly_mulmod(tmp, e0, h0, P, ctx_q);
    fq_poly_mulmod(tmp2, e1, h1, P, ctx_q);
    fq_poly_add(synd, tmp, tmp2, ctx_q);

    fq_struct* tmp0 = _fq_vec_init(r, ctx_q);
    fq_struct* tmp1 = _fq_vec_init(r, ctx_q);

    if (fq_poly_equal(sp, synd, ctx_q)) {
	fq_poly_get_coeffs(tmp0, e0, r, ctx_q);
	fq_poly_get_coeffs(tmp1, e1, r, ctx_q);
	for (i = 0; i < n; ++i) {
	    if (i < r) {
		fq_set(&res[i], &tmp0[i], ctx_q);
            } else {
		fq_set(&res[i], &tmp1[i-r], ctx_q);
            }
        }
	b = 1;
    }

    
    _fq_vec_clear(tmp0, r, ctx_q);
    _fq_vec_clear(tmp1, r, ctx_q);
    fq_poly_clear(e0, ctx_q);
    fq_poly_clear(e1, ctx_q);
    fq_poly_clear(tmp, ctx_q);
    fq_poly_clear(tmp2, ctx_q);
    fq_poly_clear(sp, ctx_q);
    fq_poly_clear(P, ctx_q);
    
    return b;
}




/**
 * Decoding as in HQC.
 * We are doing it with closest vector : dimensions allow it.
 */
void
RM_decoding_duplicated(fq_struct* res, const fq_struct* c, const int m, const int r,
		       const fq_ctx_t ctx) {

    int i, j, k, dist, min_dist, len;
    len = r * (1<<m);
    fq_struct* tmp_m = _fq_vec_init(m + 1, ctx);
    fq_struct* tmp_c = _fq_vec_init(len, ctx);

    /* initialisation at 0 */
    _fq_vec_zero(tmp_m, m+1, ctx);
    RM_encoding_duplicated(tmp_c, tmp_m, m, r, ctx);
    min_dist = hamming_distance(c, tmp_c, len, ctx);
    _fq_vec_set(res, tmp_m, m+1, ctx);

    for (k=1; k < (1<<(m+1)); ++k) {
	for (i = 0; i < m+1; i++) {
	    fq_set_ui(&tmp_m[i], (k>>i)&1, ctx);
	}

	RM_encoding_duplicated(tmp_c, tmp_m, m, r, ctx);

	dist = hamming_distance(c, tmp_c, len, ctx);
	    
	if (dist < min_dist) {
	    min_dist = dist;
	    _fq_vec_set(res, tmp_m, m+1, ctx);
	}

	if (min_dist == 0) {
	    break;
	}
    }
}


void
RS_RM_concatenated_decoding(fq_poly_t res, const fq_struct* c,
				 const fq_struct* alpha, const int n, const int k, 
				 const fq_ctx_t ctx, const int r, const fq_ctx_t ctx_q) {

    int i, j;
    int m = fq_ctx_degree(ctx);
    int len = (1 << (m-1));
    
    fq_struct* tmp_m2 = _fq_vec_init(m, ctx_q);
    fq_struct* tmp_c2 = _fq_vec_init(r * len, ctx_q);
    fmpz* tmp_m3 = _fmpz_vec_init(m);
    fmpz_poly_t tmp_pol; fmpz_poly_init(tmp_pol);

    fq_struct* tmp_c1 = _fq_vec_init(n, ctx_q);
    
    /* will decode each block of r * 2^len bits */
    for (i = 0; i < n; i++) {
	/* Decode a block of r*len bits using RM(1, m-1) code */
	for (j = 0; j < r * len; j++) {
	    fq_set(&tmp_c2[j], &c[i * r * len + j], ctx_q);
	}
	/* printf("After one block \n"); */
	
	RM_decoding_duplicated(tmp_m2, tmp_c2, m-1, r, ctx_q);
	/* printf("After RM decoding of one block \n"); */
	/* _fq_vec_print_pretty(tmp_m2, m, ctx_q); */
	
	/* Then transform the retrieved block in an elt of F_q^m  */
	_fq_vec_2_fmpz(tmp_m3, tmp_m2, m, ctx_q);
	fmpz_poly_set_coeffs(tmp_pol, tmp_m3, m);
	fq_set_fmpz_poly(&tmp_c1[i], tmp_pol, ctx);
    }
    
    RS_decoder(res, tmp_c1, alpha, n, k, ctx);


    _fq_vec_clear(tmp_c2, r * len, ctx_q);
    _fq_vec_clear(tmp_m2, m, ctx_q);
    _fmpz_vec_clear(tmp_m3, m);
    fmpz_poly_clear(tmp_pol);
}


void
HQC_decoding(fq_poly_t res, const fq_poly_t u, const fq_poly_t v, const fq_poly_t y,
	     const fq_struct* alpha, const int n, const int n1, const int k, const int r,
	     const fq_ctx_t ctx, const fq_ctx_t ctx_q) {

    int _m = fq_ctx_degree(ctx);
    int len = (1 << (_m - 1));	/* length of RM(1, m-1) code */
    
    fq_poly_t P; fq_poly_init(P, ctx_q); 
    fq_poly_set_cyclic(P, n, ctx_q);

    fq_struct* tmp_vec = _fq_vec_init(r*len*n1, ctx_q); _fq_vec_zero(tmp_vec, r*len*n1, ctx_q);
    fq_poly_t tmp_poly; fq_poly_init(tmp_poly, ctx_q); fq_poly_zero(tmp_poly, ctx_q);
    fq_poly_mulmod(tmp_poly, u, y, P, ctx_q);
    fq_poly_add(tmp_poly, tmp_poly, v, ctx_q);
    fq_poly_truncate(tmp_poly, r*len*n1, ctx_q);
    fq_poly_get_coeffs(tmp_vec, tmp_poly, r*len*n1, ctx_q);    
    

    RS_RM_concatenated_decoding(res, tmp_vec, alpha, n1, k, ctx, r, ctx_q);

    fq_poly_clear(P, ctx_q);
    fq_poly_clear(tmp_poly, ctx_q);
    _fq_vec_clear(tmp_vec, r*len*n1, ctx_q);
}
