#ifndef CODES_H
#define CODES_H

#include "tools.hpp"

int
RS_decoder(fq_poly_t res, const fq_struct* c, const fq_struct* alpha, const int n,
	   const int k, const fq_ctx_t ctx);

int
GRS_decoder(fq_struct* res, const fq_struct* c, const fq_struct* alpha,
	    const fq_struct* beta, const int n, const int k, const fq_ctx_t ctx);

int
Goppa_decoder(fq_struct* res, const fq_struct* c, const fq_struct* alpha,
	      const fq_poly_t g, const int n, const int k, const fq_ctx_t ctx);



int
Goppa_decoder_bin(fq_struct* res, const fq_struct* c, const fq_struct* alpha,
		  const fq_poly_t g, const int n, const int k, const fq_ctx_t ctx);

void
Goppa_codeword_random(fq_struct* res, const fq_struct* alpha, const fq_poly_t g,
		      const int len, const fq_ctx_t ctx, flint_rand_t state);

void
Goppa_codeword_random_bin(fq_struct* res, const fq_struct* alpha, const fq_poly_t g,
			  const int len, const fq_ctx_t ctx, flint_rand_t state);

void
Goppa_parity_check(fq_mat_t res, const fq_struct* alpha, const fq_poly_t g,
		   const fq_ctx_t ctx);


void
Goppa_parity_check_bin(fq_mat_t res, const fq_struct* alpha, const fq_poly_t g,
		       const fq_ctx_t ctx);


// ENCODER
void
RS_encoding(fq_struct* res, const fq_poly_t f, const fq_struct* alpha, const int n,
	    const fq_ctx_t ctx);

void
RM_encoding(fq_struct* res, const fq_struct* alpha, const int m, const fq_ctx_t ctx);

void
RM_encoding_duplicated(fq_struct* res, const fq_struct* alpha, const int m, const int r, const fq_ctx_t ctx);

void
RS_RM_concatenated_encoding(fq_struct* res, const fq_poly_t f,
			    const fq_struct* alpha, const int n, const fq_ctx_t ctx,
			    const int r, const fq_ctx_t ctx_q);

void
CM_encoding(fq_struct* res, const fq_struct* e, const fq_mat_t& T, const int len, 
	    const fq_ctx_t& ctx_q);

void 
Bike_encoding(fq_poly_t res, const int* e, const fq_poly_t h, const int r, const fq_ctx_t& ctx_q);



void
HQC_encoding(fq_poly_t res1, fq_poly_t res2, const fq_poly_t m, const fq_poly_t h,
	     const fq_poly_t s, const fq_struct* alpha, const int n, const int k,
	     const int n1, const int r, const int we, const int wr,
	     const fq_ctx_t ctx, const fq_ctx_t ctx_q, flint_rand_t state);

// DECODER
int
CM_syndrome_decoding(fq_struct* res, const fq_struct* s, const fq_struct* alpha, const fq_poly_t g,
		     const int len, const int t, const fq_ctx_t ctx);

int
CM_syndrome_decoding_bin(fq_struct* res, const fq_struct* s, const fq_struct* alpha, const fq_poly_t g,
			 const int len, const int t, const fq_ctx_t ctx);


int
Bike_decoding(fq_struct* res, const fq_struct* s, const fq_mat_t& H,
	      const int weight, const int NbIter, const int tau,
	      const fq_ctx_t ctx_q);



int
Bike_decoding_v2(fq_struct* res, const fq_struct* s, const fq_poly_t h0,
		 const fq_poly_t h1, const int r, const int weight,
		 const int NbIter, const int tau,
		 const fq_ctx_t ctx_q);


void
RM_decoding_duplicated(fq_struct* res, const fq_struct* c, const int m, const int r, const fq_ctx_t ctx);




void
RS_RM_concatenated_decoding(fq_poly_t res, const fq_struct* c,
				 const fq_struct* alpha, const int n, const int k, 
				 const fq_ctx_t ctx, const int r, const fq_ctx_t ctx_q);



void
HQC_decoding(fq_poly_t res, const fq_poly_t u, const fq_poly_t v, const fq_poly_t y,
	     const fq_struct* alpha, const int n, const int n1, const int k, const int r,
	     const fq_ctx_t ctx, const fq_ctx_t ctx_q);

#endif // CODES_H
