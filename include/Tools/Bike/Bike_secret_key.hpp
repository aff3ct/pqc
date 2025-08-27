#ifndef BIKE_SECRET_KEY_HPP
#define BIKE_SECRET_KEY_HPP

#include "../tools.hpp"
#include "../codes.hpp"


class Bike_secret_key {
private:
    int r;

    
public:
    GF2X h0;
    GF2X h1;

    Bike_secret_key(int r);
    ~Bike_secret_key();

    int get_r() const;


    /* Parameters as in the Bike specification document
       r is the length, i.e. the degree of X^r - 1
       w is the hamming weight of the polynomials h0, h1 */
    void keygen(const int w);
    
};


#endif /* BIKE_SECRET_KEY_HPP */
