#include "Tools/HQC/HQC_keygen.hpp"

/* random generation of secret key as in HQC */
void
HQC_keygen_naive(HQC_secret_key& SK, HQC_public_key& PK, const int w,
		 flint_rand_t state) {
        
    SK.keygen(w);

    PK.keygen(SK, state);

}
