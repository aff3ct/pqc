#include "Tools/Bike/Bike_keygen.hpp"

/* random generation of secret key as in Bike */
void
Bike_keygen_naive(Bike_secret_key& SK, Bike_public_key& PK, const int w)  {
    
    int res = 0;
    
    while (!res) {
	SK.keygen(w);
	res = PK.keygen(SK);
	
    }
}
