#include "CM_keygen.hpp"

/* random generation of secret key as in CM */
void
CM_keygen_naive(CM_secret_key& SK, CM_public_key& PK, const int n, const int t, const fq_ctx_t& ctx,
		flint_rand_t state) {

    int res = 0;

    while (!res) {
	SK.keygen(t, state);
	res = PK.keygen(SK, ctx);
    }
}
