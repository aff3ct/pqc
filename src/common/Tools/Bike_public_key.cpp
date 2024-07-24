#include "tools.hpp"
#include "Bike_public_key.hpp"

#include <iostream>
#include <stdlib.h>
using namespace std;

Bike_public_key::Bike_public_key(int r, fq_ctx_t* ctx_q) {
    this->r = r;
    this->ctx_q = ctx_q;

    fq_poly_init((this->h), *(this->ctx_q)); // init h
    fq_poly_zero((this->h), *(this->ctx_q)); // put it to 0
}

Bike_public_key::~Bike_public_key() {
    fq_poly_clear(this->h, *ctx_q);
}


int
Bike_public_key:: get_r() const {
    return this->r;
}


fq_ctx_t*
Bike_public_key:: get_ctx_q() const {
    return this->ctx_q;
}


int
Bike_public_key::keygen(const Bike_secret_key& sk) {
    fq_poly_t U, V, G, P;
    fq_poly_init(U, *(this->ctx_q)); // init U
    fq_poly_init(V, *(this->ctx_q)); // init V
    fq_poly_init(G, *(this->ctx_q)); // init G
    fq_poly_init(P, *(this->ctx_q)); // init P = X^r-1

    fq_poly_zero(P, *(this->ctx_q));
    fq_t tmp; fq_init(tmp, *(this->ctx_q));
    fq_one(tmp, *(this->ctx_q));
    fq_poly_set_coeff(P, this->r, tmp, *(this->ctx_q));
    fq_poly_set_coeff(P, 0, tmp, *(this->ctx_q));

       
    fq_poly_xgcd(G, U, V, sk.h0, P, *(this->ctx_q));

    fq_poly_print_pretty(G, "t", *(this->ctx_q));
    printf("\n");
    
    int b = fq_poly_is_one(G, *(this->ctx_q));

    if (b == 1) {
	fq_poly_mul(this->h, sk.h1, U, *(this->ctx_q));
    }

    fq_poly_clear(U, *(this->ctx_q));
    fq_poly_clear(V, *(this->ctx_q));
    fq_poly_clear(G, *(this->ctx_q));
    fq_poly_clear(P, *(this->ctx_q));
    
    return b;
}
