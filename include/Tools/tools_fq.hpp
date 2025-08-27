#ifndef TOOLS_FQ_H
#define TOOLS_FQ_H
/**
 * Obsolete tools functions using Flint Library
 * New GF2X implementations in tools.hpp
 */

// Flint Library
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

// GF2X
#include <NTL/GF2X.h>           // F_2 Finite field
#include <NTL/GF2E.h>           // F_2 mod(P[X]) Finite field
#include <NTL/GF2EX.h>          // F_2^m Extension Field
#include <NTL/mat_GF2.h>        // Matrix over GF2
#include <NTL/vec_GF2.h>        // GF2 Vectors
#include <NTL/vector.h>         // General vectors (necessary for GF2E vectors)
#include <NTL/mat_GF2E.h>       // Matrix over GF2X(modulo)
#include <NTL/GF2EXFactoring.h> // Factorization and irreductibility related functions

#include "Tools/tools.hpp"
#include <inttypes.h>

using namespace NTL;

/* MISCELLANEOUS */

int
hamming_weight(const fq_struct* v, const int len, const fq_ctx_t ctx);

int
hamming_distance(const fq_struct* v, const fq_struct* w, const int len, const fq_ctx_t ctx);

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



// Flint - GF2X conversion functions (for testing)
void fq_poly_to_uint64(const fq_poly_t poly, const fq_ctx_t ctx, uint64_t *result, slong *num_blocks_out);

GF2X fq_poly_to_gf2x(const fq_poly_t in_poly, const fq_ctx_t ctx);

GF2EX fq_poly_to_gf2ex(const fq_poly_t in_poly, const fq_ctx_t ctx);

void GF2EX_to_fq_poly(const GF2EX &in_vec, fq_poly_t result, int len, const fq_ctx_t ctx);

mat_GF2 fq_mat_to_gf2_mat(const fq_mat_t &M, const fq_ctx_t ctx_q);

mat_GF2E fq_mat_to_gf2e_mat(const fq_mat_t &T, const fq_ctx_t &ctx);

vec_GF2 fq_vec_to_vec_GF2(const fq_struct *fq_vec, slong len, const fq_ctx_t ctx);

void vec_GF2E_to_fq_vec(fq_struct *fq_vec, slong len, const vec_GF2E &input_vec, const fq_ctx_t ctx);

vec_GF2E fq_vec_to_vec_GF2E(const fq_struct *fq_vec, slong len, const fq_ctx_t ctx);

#endif // TOOLS_FQ_H
