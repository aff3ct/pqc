#include "HQC_secret_key.hpp"
#include "HQC_public_key.hpp"


void
HQC_keygen_naive(HQC_secret_key& SK, HQC_public_key& PK, const int w,
		 flint_rand_t state);
