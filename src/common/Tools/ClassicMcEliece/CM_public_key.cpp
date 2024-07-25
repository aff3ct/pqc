#include "../tools.hpp"
#include "CM_public_key.hpp"

#include <iostream>
#include <stdlib.h>
using namespace std;

CM_public_key::CM_public_key(int n, int m, int d, fq_ctx_t* ctx_q) {
    this->n = n;
    this->m = m;
    this->d = d;
    this->ctx_q = ctx_q;

    fq_mat_init(this->T, d * m, n - d * m, *ctx_q);
    fq_mat_zero(this->T, *ctx_q);
}

CM_public_key::~CM_public_key() {
    fq_mat_clear(this->T, *ctx_q);
}


int
CM_public_key:: get_n() const {
    return this->n;
}

int
CM_public_key:: get_m() const {
    return this->m;
}

int
CM_public_key:: get_d() const {
    return this->d;
}

fq_ctx_t*
CM_public_key:: get_ctx_q() const {
    return this->ctx_q;
}


int
CM_public_key::keygen(const CM_secret_key& sk, const fq_ctx_t& ctx) {
    fq_mat_t H, HH, T, I;

    // fq_ctx_t* ctx = sk.get_ctx();    
    /* compute parity check matrix over F_qᵐ */
    fq_mat_init(H, this->d, this->n, ctx);
    Goppa_parity_check(H, sk.get_alpha(), sk.g, ctx);
    
    /* expand parity check matrix */
    fq_mat_init(HH, (this->m)*(this->d), this->n, *(this->ctx_q));
    fq_matrix_expand(HH, H, ctx, *(this->ctx_q));
    
    int r = fq_mat_rref(HH, HH, *(this->ctx_q));
    
    fq_mat_window_init(I, HH, 0, 0, this->d * this->m,
		       (this->d) * (this->m), *(this->ctx_q));
    
    int b = fq_mat_is_one(I, *(this->ctx_q));
    
    /* int ptt = fq_mat_print_pretty(T, ctx_q); */
    if (b != 0) {
	fq_mat_window_init(T, HH, 0, (this->d) * (this->m), (this->d) * (this->m), (this->n),
			   *(this->ctx_q));
	fq_mat_set(this->T, T, *(this->ctx_q));
	fq_mat_window_clear(T, *(this->ctx_q));
    }

    /* clear matrices */
    fq_mat_clear(H, ctx);
    fq_mat_clear(HH, *(this->ctx_q));
    fq_mat_window_clear(I, *(this->ctx_q));
        
    return b;
}

// fq_mat_t
// CM_public_key:: get_T() const {
//     return (this->T);
// } 

// fq_poly_t&
// CM_public_key:: get_g() const {
//     return g;
// }
