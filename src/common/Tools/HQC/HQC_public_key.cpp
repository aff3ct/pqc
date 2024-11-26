#include "../tools.hpp"
#include "HQC_public_key.hpp"

#include <iostream>
#include <stdlib.h>
using namespace std;


HQC_public_key::HQC_public_key(int n, int n1, fq_ctx_t* ctx, fq_ctx_t* ctx_q) {
    this->n = n;
    this->n1 = n1;
    this->ctx = ctx;
    this->ctx_q = ctx_q;
    
    fq_poly_init((this->h), *(this->ctx_q)); // init h
    fq_poly_zero((this->h), *(this->ctx_q)); // put it to 0
    fq_poly_init((this->s), *(this->ctx_q)); // init s
    fq_poly_zero((this->s), *(this->ctx_q)); // put it to 0
    
    this->alpha = _fq_vec_init(n1, *ctx); // init alpha
    _fq_vec_zero(alpha, n1, *ctx);	 // put it to 0
}

HQC_public_key::~HQC_public_key() {
    fq_poly_clear(this->h, *ctx_q);
    fq_poly_clear(this->s, *ctx_q);
    _fq_vec_clear(this->alpha, n1, *ctx);
}


int
HQC_public_key:: get_n() const {
    return this->n;
}

int
HQC_public_key:: get_n1() const {
    return this->n1;
}


fq_ctx_t*
HQC_public_key:: get_ctx_q() const {
    return this->ctx_q;
}

fq_ctx_t*
HQC_public_key:: get_ctx() const {
    return this->ctx;
}


void
HQC_public_key::keygen(const HQC_secret_key& sk, flint_rand_t state) {
    /* init P and set it o X^n-1 */
    fq_poly_t P; fq_poly_init(P, *(this->ctx_q)); 
    fq_poly_set_cyclic(P, this->n, *(this->ctx_q));
    
    int e[this->n];
    for (int i = 0; i < n; ++i) {
      e[i] = 0;
    }

    fq_struct* tmp = _fq_vec_init(this->n, *(this->ctx_q));
    
    /* Computes h as a random element of F_2[X] / (X^n-1) */
    random_bits(e, this->n);
    _int_vec_2_fq(tmp, e, this->n, *(this->ctx_q));
    fq_poly_set_coeffs(this->h, tmp, n, *(this->ctx_q));

    /* Computes s = x + h*y */
    fq_poly_mulmod(s, this->h, sk.y, P, *(this->ctx_q));
    fq_poly_add(this->s, this->s, sk.x, *(this->ctx_q));

    /* Now computes elements of F_2^m for RS code */
    fq_vec_rand_distinct_2(this->alpha, this->n1, *(this->ctx), state);

    
    fq_poly_clear(P, *(this->ctx_q));
    _fq_vec_clear(tmp, n, *(this->ctx_q));
}
