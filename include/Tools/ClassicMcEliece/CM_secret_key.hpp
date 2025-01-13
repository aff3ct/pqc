#ifndef CM_SECRET_KEY_HPP
#define CM_SECRET_KEY_HPP

#include "../tools.hpp"
#include "../codes.hpp"


class CM_secret_key {
private:
    int n;
    fq_ctx_t* ctx;
    fq_struct* alpha;

public:
    fq_poly_t g;

    CM_secret_key(int m, fq_ctx_t* ctx);
    ~CM_secret_key();

    int get_n() const;
    fq_struct* get_alpha() const;
    fq_ctx_t* get_ctx() const;


  /* Parameters as in the Classic McEliece specification document
     alpha and g are assumed to be initialised
     Note that the parameter m as in the CM spec. is encapsulated within "ctx" here.
     n is the length
     d is the degree of the polynomial g defining the Goppa code */
    void keygen(const int d, flint_rand_t state);
    // fq_poly_t& get_g() const;
};


#endif /* CM_SECRET_KEY_HPP */
