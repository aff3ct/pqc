#include "CM_secret_key.hpp"

#include "tools.hpp"

#include <iostream>

CM_secret_key::CM_secret_key(int m, fq_ctx_t* ctx) {
    n = m;
    this->ctx = ctx;

    this->alpha = _fq_vec_init(n, *ctx); // init alpha
    _fq_vec_zero(alpha, n, *ctx);	 // put it to 0
    
    fq_poly_init((this->g), *ctx); // init g
    fq_poly_zero((this->g), *ctx); // put it to 0
}

CM_secret_key::~CM_secret_key() {
    std::cerr << "In SK destructor! \n" << std::endl;
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

fq_ctx_t*
CM_secret_key:: get_ctx() const {
    return this->ctx;
}

void
CM_secret_key:: keygen(const int d, flint_rand_t state) {
    /* n random elements */
    fq_vec_rand_distinct_2(this->alpha, this->n, *(this->ctx), state);
    std::cout << "after computing alpha !" <<     std::endl;
    /* irr pol over Fq (="ctx") without any root in alpha */
    cm_fq_poly_irr_pol(this->g, d, this->alpha, this->n, *(this->ctx),
		       state);
    std::cout << "after generating pol !" <<     std::endl;
}


// fq_poly_t&
// CM_secret_key:: get_g() const {
//     return g;
// }
