#include <iostream>
#include <flint/fq_poly.h>

#include "Tools/tools.hpp"
#include "Tools/Bike/Bike_secret_key.hpp"

Bike_secret_key::Bike_secret_key(int r) {
    this->r = r;
}


Bike_secret_key::~Bike_secret_key() {

}


int
Bike_secret_key::get_r() const {
    return this->r;
}


void
Bike_secret_key::keygen(const int w) {
    int e[this->r];
    
    for (int i = 0 ; i < r; ++i) {
	e[i] = 0;
    }


    // h0
    bike_gen_e(e, this->r, w/2);
    this->h0 = int_vec_to_GF2X(e,this->r);

    // h1
    for (int i = 0 ; i < r; ++i) {
	e[i] = 0;
    }
    bike_gen_e(e, this->r, w/2);
    this->h1 = int_vec_to_GF2X(e,this->r);

}
