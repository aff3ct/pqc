#ifndef TOOLS_H
#define TOOLS_H


#include <flint/flint.h>
#include <flint/fmpz.h>	/* large integers */
#include <flint/ulong_extras.h>	/* word-size integers */
#include <flint/fq.h>		/* finite fields */
#include <flint/fq_nmod.h>		/* finite fields with word-size modulus */
#include <flint/fq_poly.h>	/* pol. in finite fields */
#include <flint/fmpz_poly.h>	/* pol. in integers */
#include <flint/fmpz_vec.h>	/* vectors integers */
#include <flint/fq_vec.h>	/* vectors finite fields */
#include <flint/fq_mat.h>	/* matrix / finite fields */
#include <flint/perm.h>	/* permutations */




/* MISCELLANEOUS */
int
random_suitable_integer(const int len);

int
compute_threshold(const int w, const int ind, const int r);

int
hamming_weight(const fq_struct* v, const int len, const fq_ctx_t ctx);

int
hamming_distance(const fq_struct* v, const fq_struct* w, const int len, const fq_ctx_t ctx);


/* ERROR GENERATION */
// Bike
void
fisher_yates(int* perm, const int n);

  
// Classic McEliece
slong*
random_indices(const slong len, flint_rand_t state);

int
int_check_repeat(const int* a, const int len);

int
cm_random_indices(int* res, const int n, const int t, const int N, const int tau);

void
random_bits(int *res, const int len);

void
cm_gen_e(int* e, const int n, const int t, const int N, const int tau);

void
bike_gen_e(int* e, const int len, const int t);

// HQC
void
hqc_gen_e(int* e, const int len, const int t);


/* FINITE FIELDS */
void
_fq_vec_print_pretty(const fq_struct* v, const int len, const fq_ctx_t ctx);

void
_fq_vec_2_int(int* res, const fq_struct* a, const int len, const fq_ctx_t ctx);

void
_int_vec_2_fq(fq_struct* res, const int* a, const int len, const fq_ctx_t ctx);

void
_fmpz_vec_2_fq(fq_struct* res, const fmpz* a, const int len, const fq_ctx_t ctx);

void
_fq_vec_2_fmpz(fmpz* res, const fq_struct* a, const int len, const fq_ctx_t ctx);


void fq_get_coeffs(fq_struct *res, const fq_struct a, const int m,
                   const fq_ctx_t ctx, const fq_ctx_t ctx_q);

int
fq_check_repeat(const fq_struct* a, const int len, const fq_ctx_t ctx);

void
fq_vec_rand(fq_struct* res, const int len, const fq_ctx_t ctx,
	    flint_rand_t state);
void
fq_vec_rand_distinct(fq_struct* res, const int len, const fq_ctx_t ctx,
		     flint_rand_t state);

void
fq_vec_rand_distinct_2(fq_struct* res, const int len, const fq_ctx_t ctx,
		       flint_rand_t state);

void
fq_vec_shorten(fq_struct *res, const fq_struct *v, const int len,
	       const fq_ctx_t ctx);

void
fq_vec_expand(fq_struct *res, const fq_struct *v, const int len1, const int len2,
	      const fq_ctx_t ctx);


/* POLYNOMIALS */
void
fq_poly_set_coeffs(fq_poly_t f, const fq_struct* alpha, const int len,
		   const fq_ctx_t ctx);

void
fq_poly_get_coeffs(fq_struct* res, const fq_poly_t f, const int len, const fq_ctx_t ctx);

void
fmpz_poly_set_coeffs(fmpz_poly_t f, const fmpz* a, const int len);
  
void
fmpz_poly_get_coeffs(fmpz* res, const fmpz_poly_t f, const int len);

  
void
fq_poly_nonzero_coeffs(int* pos, const fq_poly_t f, const int r,
			    const fq_ctx_t ctx);

void
fq_poly_set_cyclic(fq_poly_t res, const int d, const fq_ctx_t ctx);


void
fq_poly_set_linear(fq_poly_t res, const fq_t alpha, const fq_ctx_t ctx);


void
fq_poly_set_linear_product(fq_poly_t res, const fq_struct* alpha, const int len,
			   const fq_ctx_t ctx);

int
fq_poly_eval_zero(const fq_poly_t f, const fq_struct *alpha, const int len,
		  const fq_ctx_t ctx);

void
fq_poly_interpolate(fq_poly_t res, const fq_struct* alpha, const fq_struct*  beta,
		    const int len, const fq_ctx_t ctx);

void
cm_fq_poly_irr_pol(fq_poly_t& res, const int deg, const fq_struct* alpha, const int len,
			const fq_ctx_t ctx, flint_rand_t state);

void
xgcd_abort(fq_poly_t u, fq_poly_t v, fq_poly_t d, const fq_poly_t a, const fq_poly_t b,
		const slong k, const fq_ctx_t ctx);




 /* MATRICES */
void
fq_matrix_expand(fq_mat_t res, const fq_mat_t H, const fq_ctx_t ctx,
		      const fq_ctx_t ctx_q);


void
fq_mult_matrix(fq_mat_t res, const fq_poly_t h, const fq_poly_t P, const fq_ctx_t ctx);



/* BIKE */
int
ctr(const fq_struct *v, const fq_mat_t &H,  const int j, const fq_ctx_t ctx);


int
ctrv2(const fq_poly_t sp,  const int* pos, const int j, const int weight, const int r,
      const fq_ctx_t ctx);


void
BFIter(fq_struct* e, int* black, int* gray, const fq_struct* s,
       const fq_mat_t& H, const int T, const int tau, const fq_ctx_t ctx);



void
BFIterv2(fq_poly_t e0, fq_poly_t e1, int* black, int* gray, const fq_poly_t sp,
	 const int* pos0, const int* pos1, const int weight, const int r, const int T,
	 const int tau, const fq_ctx_t ctx);


void
BFMaskedIter(fq_struct* e, const fq_struct* s, const fq_mat_t& H, const int T,
	     const int* mask, const fq_ctx_t ctx);


void
BFMaskedIterv2(fq_poly_t e0, fq_poly_t e1, const fq_poly_t sp, const int* pos0,
	       const int* pos1, const int weight, const int r, const int T,
	       const int* mask, const fq_ctx_t ctx);


void
CM_params(int &m, int &n, int &t, const int level);

void
Bike_params(int& r, int& weight, int& error_weight, const int level);


void
BGF_params(int& NbIter, int& tau, const int level);


void
HQC_params(int& k, int &n1, int &n2, int &r, int &n, int &w, int &we, int &wr,
	   const int level);

#endif // TOOLS_H
