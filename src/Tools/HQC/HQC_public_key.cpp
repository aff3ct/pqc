#include <iostream>
#include <stdlib.h>

#include "Tools/tools.hpp"
#include "Tools/HQC/HQC_public_key.hpp"

using namespace std;

HQC_public_key::HQC_public_key(int n, int n1)
{
    this->n = n;
    this->n1 = n1;

    this->alpha.SetMaxLength(n1);
    this->alpha.SetLength(n1);

    this->h.SetMaxLength(n);
    this->h.SetLength(n);

    this->s.SetMaxLength(n);
    this->s.SetLength(n);
}

HQC_public_key::~HQC_public_key()
{
}

int HQC_public_key::get_n() const
{
    return this->n;
}

int HQC_public_key::get_n1() const
{
    return this->n1;
}

void HQC_public_key::keygen(const HQC_secret_key &sk, flint_rand_t state)
{
    /* init P and set it to X^n-1 */
    GF2X P = gf2x_set_cyclic(this->n);

    int e[this->n];
    for (int i = 0; i < n; ++i)
    {
        e[i] = 0;
    }

    /* Computes h as a random element of F_2[X] / (X^n-1) */
    random_bits(e, this->n);

    vec_rand_distinct_2(this->alpha, this->n1, state);

    for (int i = 0; i < this->n; i++)
    {
        SetCoeff(h, i, e[i]);
    }

    s = MulMod(h, sk.y, P);
    s += sk.x;
}
