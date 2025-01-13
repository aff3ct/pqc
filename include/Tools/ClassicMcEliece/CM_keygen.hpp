#include "CM_secret_key.hpp"
#include "CM_public_key.hpp"


void CM_keygen_naive(CM_secret_key& SK, CM_public_key& PK, const int n, const int t, const fq_ctx_t& ctx,
		     flint_rand_t state);
