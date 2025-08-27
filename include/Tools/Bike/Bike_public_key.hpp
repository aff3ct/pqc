#ifndef BIKE_PUBLIC_KEY_HPP
#define BIKE_PUBLIC_KEY_HPP

#include "../tools.hpp"
#include "../codes.hpp"
#include "Bike_secret_key.hpp"

class Bike_public_key {
private:
    int r;

public:
    GF2X h;

    Bike_public_key(int r);
    ~Bike_public_key();

    int get_r() const;
    
    int keygen(const Bike_secret_key& sk);
};


#endif /* BIKE_PUBLIC_KEY_HPP */
