#ifndef BIKE_SECRET_KEY_HPP
#define BIKE_SECRET_KEY_HPP

#include "tools.hpp"
#include "codes.hpp"


class Bike_secret_key {
private:
    int r;
    fq_ctx_t* ctx_q;
    

public:
    fq_poly_t h0;
    fq_poly_t h1;

    Bike_secret_key(int r, fq_ctx_t* ctx_q);
    ~Bike_secret_key();

    int get_r() const;
    fq_ctx_t* get_ctx_q() const;


    /* Parameters as in the Bike specification document
       r is the length, i.e. the degree of X^r - 1
       w is the hamming weight of the polynomials h0, h1 */
    void keygen(const int w);
    
};


#endif /* BIKE_SECRET_KEY_HPP */
