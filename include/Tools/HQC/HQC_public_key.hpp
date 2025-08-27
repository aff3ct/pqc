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
    

public:    
    GF2X h;
    GF2X s;
    vec_GF2E alpha;
    
    HQC_public_key(int n, int n1);
    ~HQC_public_key();

    int get_n() const;
    int get_n1() const;
    
    void keygen(const HQC_secret_key& sk, flint_rand_t state);
};


#endif /* HQC_PUBLIC_KEY_HPP */
