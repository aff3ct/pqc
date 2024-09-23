#ifndef HQC_SECRET_KEY_HPP
#define HQC_SECRET_KEY_HPP

#include "../tools.hpp"
#include "../codes.hpp"


class HQC_secret_key {
private:
    int n;
    fq_ctx_t* ctx_q;

    
public:
    fq_poly_t x;
    fq_poly_t y;

    HQC_secret_key(int n, fq_ctx_t* ctx_q);
    ~HQC_secret_key();

    int get_n() const;
    fq_ctx_t* get_ctx_q() const;


    /* Parameters as in the HQC specification document
       n is the length, i.e. the degree of X^n - 1
       w is the hamming weight of the polynomials x, y */
    void keygen(const int w);
    
};


#endif /* HQC_SECRET_KEY_HPP */
