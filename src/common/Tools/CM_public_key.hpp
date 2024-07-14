#ifndef CM_PUBLIC_KEY_HPP
#define CM_PUBLIC_KEY_HPP

#include "tools.hpp"
#include "codes.hpp"
#include "CM_secret_key.hpp"

class CM_public_key {
private:
    int n;
    int m;
    int t;
    fq_ctx_t* ctx_q;

public:
    fq_mat_t T;

    CM_public_key(int n, int m, int t, fq_ctx_t* ctx_q);
    ~CM_public_key();

    int get_n() const;
    int get_m() const;
    int get_t() const;


    
    int keygen(const CM_secret_key sk, const fq_ctx_t ctx);
    // fq_mat_t get_T() const;
};


#endif /* CM_PUBLIC_KEY_HPP */
