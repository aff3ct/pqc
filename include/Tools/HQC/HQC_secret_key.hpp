#ifndef HQC_SECRET_KEY_HPP
#define HQC_SECRET_KEY_HPP

#include "../tools.hpp"
#include "../codes.hpp"


class HQC_secret_key {
    /* Parameters as in the HQC specification document
       n is the length, i.e. the degree of X^n - 1
       w is the hamming weight of the polynomials x, y */

private:
    int n;
    
public:
    GF2X x;
    GF2X y;

    HQC_secret_key(int n);
    ~HQC_secret_key();

    int get_n() const;

    void keygen(const int w);
    
};

#endif /* HQC_SECRET_KEY_HPP */
