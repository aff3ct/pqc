#ifndef HQC_PUBLIC_KEY_HPP
#define HQC_PUBLIC_KEY_HPP

#include "../tools.hpp"
#include "../codes.hpp"
#include "HQC_secret_key.hpp"

class HQC_public_key {
    /* Parameters as in the HQC specification document
       n is the length, i.e. the degree of X^n - 1
       n1 is the length of the shortened RS code
       w is the hamming weight of the polynomials x, y */

    
private:
    int n;      
    int n1;
    
    fq_ctx_t* ctx;
    fq_ctx_t* ctx_q;
    
    

public:
    fq_poly_t h;
    fq_poly_t s;
    fq_struct* alpha;
    
    HQC_public_key(int n, int n1, fq_ctx_t* ctx, fq_ctx_t* ctx_q);
    ~HQC_public_key();

    int get_n() const;
    int get_n1() const;
    fq_ctx_t* get_ctx_q() const;
    fq_ctx_t* get_ctx() const;

    
    void keygen(const HQC_secret_key& sk, flint_rand_t state);
};


#endif /* HQC_PUBLIC_KEY_HPP */
