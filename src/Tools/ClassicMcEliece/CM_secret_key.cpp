#include <iostream>
#include "Tools/tools.hpp"
#include "Tools/ClassicMcEliece/CM_secret_key.hpp"

CM_secret_key::CM_secret_key(int m) {
    n = m;
    this->alpha.SetMaxLength(n);
    this->alpha.SetLength(n);

}

CM_secret_key::~CM_secret_key() {

}

int
CM_secret_key:: get_n() const {
    return this->n;
}


vec_GF2E CM_secret_key:: get_alpha() const {
    return this->alpha;
}


void
CM_secret_key:: keygen(const int d, flint_rand_t state) {
    /* n random elements */
    vec_rand_distinct_2(this->alpha, this->n,state);

    cm_poly_irr_pol(this->g,d,this->alpha,this->n);
}
