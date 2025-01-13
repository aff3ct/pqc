#ifndef BIKE_PUBLIC_KEY_HPP
#define BIKE_PUBLIC_KEY_HPP

#include "../tools.hpp"
#include "../codes.hpp"
#include "Bike_secret_key.hpp"

class Bike_public_key {
private:
    int r;
    fq_ctx_t* ctx_q;

public:
    fq_poly_t h;

    Bike_public_key(int r, fq_ctx_t* ctx_q);
    ~Bike_public_key();

    int get_r() const;
    fq_ctx_t* get_ctx_q() const;

    
    int keygen(const Bike_secret_key& sk);
};


#endif /* BIKE_PUBLIC_KEY_HPP */
