#ifndef CODES_H
#define CODES_H

#include "tools.hpp"

int
RS_decoder(fq_poly_t res, const fq_struct* c, const fq_struct* alpha, const int n, const int k,
	       const fq_ctx_t ctx);

int
GRS_decoder(fq_struct* res, const fq_struct* c, const fq_struct* alpha, const fq_struct* beta,
	    const int n, const int k, const fq_ctx_t ctx);

int
Goppa_decoder(fq_struct* res, const fq_struct* c, const fq_struct* alpha, const fq_poly_t g,
	      const int n, const int k, const fq_ctx_t ctx);



int
Goppa_decoder_bin(fq_struct* res, const fq_struct* c, const fq_struct* alpha, const fq_poly_t g,
		  const int n, const int k, const fq_ctx_t ctx);

void
Goppa_codeword_random(fq_struct* res, const fq_struct* alpha, const fq_poly_t g, const int len,
     		      const fq_ctx_t ctx, flint_rand_t state);

void
Goppa_codeword_random_bin(fq_struct* res, const fq_struct* alpha, const fq_poly_t g, const int len,
			  const fq_ctx_t ctx, flint_rand_t state);

void
Goppa_parity_check(fq_mat_t res, const fq_struct* alpha, const fq_poly_t g, const fq_ctx_t ctx);


void
Goppa_parity_check_bin(fq_mat_t res, const fq_struct* alpha, const fq_poly_t g, const fq_ctx_t ctx);


// ENCODER
void
CM_encoding(fq_struct* res, const fq_struct* e, const fq_mat_t& T, const int len, 
	    const fq_ctx_t& ctx_q);

void 
Bike_encoding(fq_poly_t res, const int* e, const fq_poly_t h, const int r, const fq_ctx_t& ctx_q);


// DECODER
int
CM_syndrome_decoding(fq_struct* res, const fq_struct* s, const fq_struct* alpha, const fq_poly_t g,
		     const int len, const int t, const fq_ctx_t ctx);

int
CM_syndrome_decoding_bin(fq_struct* res, const fq_struct* s, const fq_struct* alpha, const fq_poly_t g,
			 const int len, const int t, const fq_ctx_t ctx);

#endif // CODES_H

int
Bike_decoding(fq_struct* res, const fq_struct* s, const fq_mat_t& H,
	      const int weight, const int NbIter, const int tau,
	      const fq_ctx_t ctx_q);
