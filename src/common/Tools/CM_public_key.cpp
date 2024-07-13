#include "tools.hpp"
#include "CM_public_key.hpp"


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

// fq_mat_t
// CM_public_key:: get_T() const {
//     return (this->T);
// } 

// fq_poly_t&
// CM_public_key:: get_g() const {
//     return g;
// }
