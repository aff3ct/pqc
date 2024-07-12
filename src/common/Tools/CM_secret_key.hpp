#ifndef CM_SECRET_KEY_HPP
#define CM_SECRET_KEY_HPP

#include "tools.hpp"
#include "codes.hpp"


class CM_secret_key {
private:
    int n;
    fq_ctx_t* ctx;
    fq_struct* alpha;

public:
    CM_secret_key(int m, fq_ctx_t* ctx);
    ~CM_secret_key();

    int get_n() const;
    fq_struct* get_alpha() const;

    fq_poly_t g;
    // fq_poly_t& get_g() const;
};


#endif /* CM_SECRET_KEY_HPP */
