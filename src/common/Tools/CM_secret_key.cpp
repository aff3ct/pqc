#include "CM_secret_key.hpp"

#include "tools.hpp"

CM_secret_key::CM_secret_key(int m, fq_ctx_t* ctx) {
    n = m;
    this->ctx = ctx;

    this->alpha = _fq_vec_init(n, *ctx); // init alpha
    _fq_vec_zero(alpha, n, *ctx);	 // put it to 0
    
    fq_poly_init((this->g), *ctx); // init g
    fq_poly_zero((this->g), *ctx); // put it to 0
}

CM_secret_key::~CM_secret_key() {
    _fq_vec_clear(this->alpha, this->n, *(this->ctx));
    fq_poly_clear((this->g), *(this->ctx));
}


int
CM_secret_key:: get_n() const {
    return this->n;
}

fq_struct*
CM_secret_key:: get_alpha() const {
    return this->alpha;
}

// fq_poly_t&
// CM_secret_key:: get_g() const {
//     return g;
// }
