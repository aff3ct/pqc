#include "HQC_secret_key.hpp"

#include "../tools.hpp"

#include <flint/fq_poly.h>

#include <iostream>

HQC_secret_key::HQC_secret_key(int n, fq_ctx_t* ctx_q) {
    this->n = n;
    this->ctx_q = ctx_q;

    fq_poly_init((this->x), *(this->ctx_q)); // init x
    fq_poly_zero((this->x), *(this->ctx_q)); // put it to 0

    fq_poly_init((this->y), *(this->ctx_q)); // init y
    fq_poly_zero((this->y), *(this->ctx_q)); // put it to 0
}


HQC_secret_key::~HQC_secret_key() {
    fq_poly_clear((this->x), *(this->ctx_q));
    fq_poly_clear((this->y), *(this->ctx_q));
}


int
HQC_secret_key::get_n() const {
    return this->n;
}


fq_ctx_t*
HQC_secret_key::get_ctx_q() const {
    return this->ctx_q;
}


void
HQC_secret_key::keygen(const int w) {
    int e[this->n];
    
    for (int i = 0 ; i < n; ++i) {
	e[i] = 0;
    }

    fq_struct* tmp = _fq_vec_init(this->n, *(this->ctx_q));

    // x
    hqc_gen_e(e, this->n, w);
    _int_vec_2_fq(tmp, e, this->n, *(this->ctx_q));
    fq_poly_set_coeffs(this->x, tmp, n, *(this->ctx_q));

    // y
    for (int i = 0 ; i < n; ++i) {
	e[i] = 0;
    }
    hqc_gen_e(e, this->n, w);
    _int_vec_2_fq(tmp, e, this->n, *(this->ctx_q));
    fq_poly_set_coeffs(this->y, tmp, n, *(this->ctx_q));

    _fq_vec_clear(tmp, n, *(this->ctx_q));
}
