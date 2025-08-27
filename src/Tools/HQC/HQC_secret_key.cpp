#include <iostream>
#include <flint/fq_poly.h>
#include "Tools/tools.hpp"
#include "Tools/HQC/HQC_secret_key.hpp"

HQC_secret_key::HQC_secret_key(int n) {
    this->n = n;

    // INIT _x,_y? maybe not necessary 
}


HQC_secret_key::~HQC_secret_key() {

}


int
HQC_secret_key::get_n() const {
    return this->n;
}



void
HQC_secret_key::keygen(const int w) {
    
    int e[this->n];
    
    for (int i = 0 ; i < this->n; ++i) {
	e[i] = 0;
    }

    // x
    hqc_gen_e(e, this->n, w);
    x = int_vec_to_GF2X(e,this->n);

    // y
    for (int i = 0 ; i < this->n; ++i) {
	e[i] = 0;
    }
    hqc_gen_e(e, this->n, w);
    y = int_vec_to_GF2X(e,this->n);
}
