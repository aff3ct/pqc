#ifndef HQC_PUBLIC_KEY_HPP
#define HQC_PUBLIC_KEY_HPP

#include "../tools.hpp"
#include "../codes.hpp"
#include "HQC_secret_key.hpp"

class HQC_public_key {
private:
    int r;
    fq_ctx_t* ctx_q;

public:
    fq_poly_t h;

    HQC_public_key(int r, fq_ctx_t* ctx_q);
    ~HQC_public_key();

    int get_r() const;
    fq_ctx_t* get_ctx_q() const;

    
    int keygen(const HQC_secret_key& sk);
};


#endif /* HQC_PUBLIC_KEY_HPP */
