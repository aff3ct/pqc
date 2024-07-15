#include "codes.hpp"


// /* **************************************************************************** */
// /*                            structures of keys                                */
// /* **************************************************************************** */


// /* structure of a secret key
//    alpha : Vector of elements in a finite field FF_{q^m}
//    g : irreducible polynomial over FF_{q^m} */
// struct CM_sk {
//     fq_struct* alpha;
//     fq_poly_t g;
// };


// /* initialisation of secret key */
// void
// init_CM_sk(struct CM_sk* sk, const int n, const fq_ctx_t ctx) {
//     sk->alpha = _fq_vec_init(n, ctx);
//     fq_poly_init(sk->g, ctx);
// }

// /* memory clear of secret key */
// void
// clear_CM_sk(struct CM_sk* sk, const int n, const fq_ctx_t ctx) {
//     _fq_vec_clear(sk->alpha, n, ctx);
//     fq_poly_clear(sk->g, ctx);
// }


// /* structure of public key
//    matrix of elements in a finite field
// */
// struct CM_pk  {
//     fq_mat_t T;			/* this is the true PK */
// };


// /* initialisation of public key */
// void
// init_CM_pk(struct CM_pk* pk, const int n,  const int t, const fq_ctx_t ctx,
// 	const fq_ctx_t ctx_q) {
//     fq_mat_init(pk->T, t * fq_ctx_degree(ctx), n - t * fq_ctx_degree(ctx), ctx_q);
// }

// /* memory clear of public key */
// void
// clear_CM_pk(struct CM_pk* pk, const int n,  const int t, const fq_ctx_t ctx,
// 	 const fq_ctx_t ctx_q) {
//     fq_mat_clear(pk->T, ctx_q);
// }



/* **************************************************************************** */
/*                                   DECODERS                                   */
/* **************************************************************************** */



/* Decoding algorithm for Reed-Solomon codes
   Decode received message c wrt to RS_k(alpha)
*/
int
RS_decoder(fq_poly_t res, const fq_struct* c, const fq_struct* alpha, const int n, const int k,
	   const fq_ctx_t ctx, const fq_ctx_t ctx_q) {

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
	fq_neg(a, &
	       alpha[i], ctx);
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


/* Decoding algorithm for Generalised Reed-Solomon codes */
int
GRS_decoder(fq_struct* res, const fq_struct* c, const fq_struct* alpha, const fq_struct* beta,
	    const int n, const int k, const fq_ctx_t ctx,
	    const fq_ctx_t ctx_q) {
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
    int b = RS_decoder(p, c1, alpha, n, k, ctx, ctx_q);

    if (b) {
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
	      const int n, const int k, const fq_ctx_t ctx, const fq_ctx_t ctx_q) {

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

    int b = GRS_decoder(res, c, alpha, beta, n, k, ctx, ctx_q);

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
Goppa_decoder_bin(fq_struct* res, const fq_struct* c, const fq_struct* alpha, const fq_poly_t g,
		  const int n, const int k, const fq_ctx_t ctx, const fq_ctx_t ctx_q) {

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
    int b = GRS_decoder(res, c, alpha, beta, n, k, ctx, ctx_q);
    return b;
}


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
Goppa_codeword_random_bin(fq_struct* res, const fq_struct* alpha, const fq_poly_t g, const int len,
			  const fq_ctx_t ctx, flint_rand_t state) {
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



/* computes the parity check matrix of a Goppa code knowing it is *BINARY* */
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


/* encoding as in CM
   Compute the syndrome associated s = H * e associated to e
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
