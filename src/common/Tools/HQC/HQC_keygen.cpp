#include "HQC_keygen.hpp"

/* random generation of secret key as in HQC */
void
HQC_keygen_naive(HQC_secret_key& SK, HQC_public_key& PK, const int w)  {
    
    int res = 0;
    
    while (!res) {
	SK.keygen(w);
	res = PK.keygen(SK);
	
    }
}
