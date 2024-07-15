#include "tools.hpp"
#include "CM_public_key.hpp"

#include <iostream>
#include <stdlib.h>
using namespace std;

CM_public_key::CM_public_key(int n, int m, int t, fq_ctx_t* ctx_q) {
    this->n = n;
    this->m = m;
    this->t = t;
    this->ctx_q = ctx_q;

    fq_mat_init(this->T, t * m, n - t * m, *ctx_q);
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
CM_public_key:: get_t() const {
    return this->t;
}

int
CM_public_key::keygen(const CM_secret_key& sk, const fq_ctx_t& ctx) {
    fq_mat_t H, HH, T, I;

    // fq_ctx_t* ctx = sk.get_ctx();    
    /* compute parity check matrix over F_qᵐ */
    fq_mat_init(H, this->t, this->n, ctx);
    Goppa_parity_check(H, sk.get_alpha(), sk.g, ctx);
    
    /* expand parity check matrix */
    fq_mat_init(HH, (this->m)*(this->t), this->n, *(this->ctx_q));
    fq_matrix_expand(HH, H, ctx, *(this->ctx_q));
    
    int r = fq_mat_rref(HH, HH, *(this->ctx_q));
    
    fq_mat_window_init(I, HH, 0, 0, this->t * this->m,
		       (this->t) * (this->m), *(this->ctx_q));
    
    int b = fq_mat_is_one(I, *(this->ctx_q));
    
    /* printf("%d\n", b); */
    /* int ptt = fq_mat_print_pretty(T, ctx_q); */
    if (b) {
	fq_mat_window_init(T, HH, 0, (this->t) * (this->m),
			   (this->t) * (this->m), (this->n),
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
