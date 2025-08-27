#ifndef TOOLS_H
#define TOOLS_H

#include <flint/flint.h>
#include <flint/fmpz.h>			/* large integers */
#include <flint/ulong_extras.h> /* word-size integers */
#include <flint/fq.h>			/* finite fields */
#include <flint/fq_nmod.h>		/* finite fields with word-size modulus */
#include <flint/fq_poly.h>		/* pol. in finite fields */
#include <flint/fmpz_poly.h>	/* pol. in integers */
#include <flint/fmpz_vec.h>		/* vectors integers */
#include <flint/fq_vec.h>		/* vectors finite fields */
#include <flint/fq_mat.h>		/* matrix / finite fields */
#include <flint/perm.h>			/* permutations */

#include <NTL/GF2X.h>           // F_2 Finite field
#include <NTL/GF2E.h>           // F_2 mod(P[X]) Finite field
#include <NTL/GF2EX.h>          // F_2^m Extension Field
#include <NTL/mat_GF2.h>        // Matrix over GF2
#include <NTL/vec_GF2.h>        // GF2 Vectors
#include <NTL/vector.h>         // General vectors (necessary for GF2E vectors)
#include <NTL/mat_GF2E.h>       // Matrix over GF2X(modulo)
#include <NTL/GF2EXFactoring.h> // Factorization and irreductibility related functions

using namespace NTL;

/* MISCELLANEOUS */
void extract_position(uint64_t* poly_in, int* pos ,int length);

void bit_to_uint(uint64_t *output, int *input, int length);

int random_suitable_integer(const int len);// Flint

int compute_threshold(const int w, const int ind, const int r);

/* ERROR GENERATION */
// Bike
void fisher_yates(int *perm, const int n);

// Classic McEliece
slong* random_indices(const slong len, flint_rand_t state);

int int_check_repeat(const int *a, const int len);

int cm_random_indices(int *res, const int n, const int t, const int N, const int tau);

void random_bits(int *res, const int len);

void random_bytes(int *res, const int len);

void cm_gen_e(int *e, const int n, const int t, const int N, const int tau);

void bike_gen_e(int *e, const int len, const int t);


void hqc_gen_e(int *e, const int len, const int t);

/* FINITE FIELDS */
void vec_rand_distinct_2(vec_GF2E &res, const int len,
							  flint_rand_t state);

int gf2ex_eval_zero(const GF2EX &f, const vec_GF2E &alpha, const int len);

void cm_poly_irr_pol(GF2EX &res, const int deg, const vec_GF2E &alpha, const int len);

void XGCD_abort(GF2EX &d, GF2EX &s, GF2EX &t, const GF2EX &a, const GF2EX &b, long k);

void rref(mat_GF2 &A);

/* MATRICES */

/* BIKE */
int ctrv2(const GF2X &sp, const int *pos, const int j, const int weight, const int r);

int ctrv2_v2(unsigned char *spBytes, const int *pos, const int j, const int weight, const int r);


void BFIterv2(GF2X &e0, GF2X &e1, int *black, int *gray, const GF2X &sp,
				   const int *pos0, const int *pos1, const int weight, const int r, const int T,
				   const int tau);

void BFMaskedIterv2(GF2X &e0, GF2X &e1, const GF2X &sp, const int *pos0,
						 const int *pos1, const int weight, const int r, const int T,
						 const int *mask);

void CM_params(int &m, int &n, int &t, GF2X &P, const int level);

void Bike_params(int &r, int &weight, int &error_weight, const int level);

void BGF_params(int &NbIter, int &tau, const int level);

void HQC_params(int &k, int &n1, int &n2, int &r, int &n, int &w, int &we, int &wr,
				const int level);

/* OTHER */

mat_GF2 concat_hor_mat_GF2(const mat_GF2 &M1, const mat_GF2 &M2);

mat_GF2 matrix_expand(mat_GF2E &M);

void vec_gf2_to_int(const Vec<GF2> &in_vec, int *result);

void vec_gf2e_to_int(const Vec<GF2E> &in_vec, int *result); // TO CHANGE

void GF2EX_to_int(const GF2EX &in_vec, int *result, int len);

void HQC_GF2EX_to_int(const GF2EX &in_pol, int *result, int len);

void GF2X_to_int(const GF2X &in_vec, int *result, int len);

GF2EX uint8_vec_to_GF2EX(int *in_vec, long len);

GF2EX uint8_vec_to_GF2EX_V2(unsigned char *in_vec, long len);

GF2EX int_vec_to_GF2EX(int *in_vec, long len);

GF2X int_vec_to_GF2X(int *in_vec, long len);

GF2X gf2x_set_cyclic(const int d);

#endif // TOOLS_H
