#include "Tools/tools.hpp"
#include "Tools/Bike/Bike_public_key.hpp"

#include <iostream>
#include <stdlib.h>
using namespace std;

Bike_public_key::Bike_public_key(int r) {
    this->r = r;

}

Bike_public_key::~Bike_public_key() {
}


int
Bike_public_key:: get_r() const {
    return this->r;
}


int
Bike_public_key::keygen(const Bike_secret_key& sk) {
    GF2X U,V,G,P;
    P = gf2x_set_cyclic(this->r);

    //GCD
    XGCD(G,U,V,sk.h0,P);

    int b = IsOne(G);

    if(b){
        this->h = sk.h1 * U;
    }
    
    return b;
}
