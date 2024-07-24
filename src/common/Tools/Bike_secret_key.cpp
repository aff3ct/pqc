#include "Bike_secret_key.hpp"

#include "tools.hpp"

#include <flint/fq_poly.h>

#include <iostream>

Bike_secret_key::Bike_secret_key(int r, fq_ctx_t* ctx_q) {
    this->r = r;
    this->ctx_q = ctx_q;

    fq_poly_init((this->h0), *(this->ctx_q)); // init h0
    // _fq_poly_set_length(this->h0, this->r, *(this->ctx_q));
    fq_poly_zero((this->h0), *(this->ctx_q)); // put it to 0

    fq_poly_init((this->h1), *(this->ctx_q)); // init h1
    // _fq_poly_set_length(this->h1, this->r, *(this->ctx_q));
    fq_poly_zero((this->h1), *(this->ctx_q)); // put it to 0
}


Bike_secret_key::~Bike_secret_key() {
    fq_poly_clear((this->h0), *(this->ctx_q));
    fq_poly_clear((this->h1), *(this->ctx_q));
}


int
Bike_secret_key::get_r() const {
    return this->r;
}


fq_ctx_t*
Bike_secret_key::get_ctx_q() const {
    return this->ctx_q;
}


void
Bike_secret_key::keygen(const int w) {
    int e[this->r];
    
    for (int i = 0 ; i < r; ++i) {
	e[i] = 0;
    }

    fq_struct* tmp = _fq_vec_init(this->r, *(this->ctx_q));

    // h0
    bike_gen_e(e, this->r, w/2);
    _int_vec_2_fq(tmp, e, this->r, *(this->ctx_q));
    fq_poly_set_coeffs(this->h0, tmp, r, *(this->ctx_q));

    // h1
    for (int i = 0 ; i < r; ++i) {
	e[i] = 0;
    }
    bike_gen_e(e, this->r, w/2);
    _int_vec_2_fq(tmp, e, this->r, *(this->ctx_q));
    fq_poly_set_coeffs(this->h1, tmp, r, *(this->ctx_q));

    _fq_vec_clear(tmp, r, *(this->ctx_q));
}
