#ifndef TOOLS_H
#define TOOLS_H


#include <flint/flint.h>
#include <flint/fmpz.h>	/* large integers */
#include <flint/fq.h>		/* finite fields */
#include <flint/fq_poly.h>	/* pol. in finite fields */
#include <flint/fmpz_poly.h>	/* pol. in integers */
#include <flint/fmpz_vec.h>	/* vectors integers */
#include <flint/fq_vec.h>	/* vectors finite fields */
#include <flint/fq_mat.h>	/* matrix / finite fields */
#include <flint/perm.h>	/* permutations */


void fisher_yates(int* perm, const int n);


// ERROR GENERATION
slong* random_indices(const slong len, flint_rand_t state);

int int_check_repeat(const int* a, const int len);

int cm_random_indices(int* res, const int n, const int t, const int N, const int tau);

void cm_gen_e(int* e, const int n, const int t, const int N, const int tau);


// FINITE FIELDS
int fq_check_repeat(const fq_struct* a, const int len, const fq_ctx_t ctx);

void fq_vec_rand_distinct(fq_struct* res, const int len, const fq_ctx_t ctx, flint_rand_t state);

void fq_vec_rand_distinct_2(fq_struct* res, const int len, const fq_ctx_t ctx, flint_rand_t state);


// POLYNOMIALS
int fq_poly_eval_zero(const fq_poly_t f, const fq_struct *alpha, const int len, const fq_ctx_t ctx);

void fq_poly_interpolate(fq_poly_t res, const fq_struct* alpha, const fq_struct*  beta, const int len,
			 const fq_ctx_t ctx);

void cm_fq_poly_irr_pol(fq_poly_t res, const int deg, const fq_struct* alpha, const int len,
			const fq_ctx_t ctx, flint_rand_t state);

void xgcd_abort(fq_poly_t u, fq_poly_t v, fq_poly_t d, const fq_poly_t a, const fq_poly_t b,
		const slong k, const fq_ctx_t ctx);




// MATRICES
void fq_matrix_expand(fq_mat_t res, const fq_mat_t H, const fq_ctx_t ctx,
		      const fq_ctx_t ctx_q);

#endif // TOOLS_H
