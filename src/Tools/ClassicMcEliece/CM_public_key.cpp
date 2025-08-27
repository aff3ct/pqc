#include <iostream>
#include <stdlib.h>
#include "Tools/tools.hpp"
#include "Tools/ClassicMcEliece/CM_public_key.hpp"

using namespace std;

CM_public_key::CM_public_key(int n, int m, int d) {
    this->n = n;
    this->m = m;
    this->d = d;
}

CM_public_key::~CM_public_key() {
}

int
CM_public_key:: get_n() const {
    return this->n;
}

int
CM_public_key:: get_m() const {
    return this->m;
}

int
CM_public_key:: get_d() const {
    return this->d;
}


int
CM_public_key::keygen(const CM_secret_key& sk) {
mat_GF2E H;
mat_GF2 HH, I, tmp;

H.SetDims(this->d, this->n);
int dim1 = (this->m) * (this->d);
int dim2 = this->n - dim1;

HH.SetDims(dim1, this->n);
I.SetDims(dim1, dim1);
T.SetDims(dim1, dim2);

H = Goppa_parity_check(sk.get_alpha(), sk.g, this->d, this->n);

HH = matrix_expand(H);

int r = gauss(HH);  /* Row echelon */
rref(HH);           /* Recuded row echelon*/

// for (int i = 0; i < dim1; i++) {
//     for (int j = 0; j < dim1; j++) {
//         I[i][j] = HH[i][j];
//     }
// }

for (int i = 0; i < dim1; i++) {
        I[i] = VectorCopy(HH[i], dim1); // truncates to dim1
}


int b = IsIdent(I, dim1);

if (b != 0) {
    for (int i = 0; i < dim1; i++) {
        for (int j = dim1; j < (this->n); j++) {
            T[i][j - dim1] = HH[i][j];
        }
    }

    // // Temporary variables for transpose
    // mat_GF2 HH_tmp = transpose(HH);
    // mat_GF2 T_tmp = transpose(T);

    // // Copy rows
    // for (int j = dim1; j < (this->n); j++) {
    //     T_tmp[j - dim1] = HH_tmp[j];
    // }
    
    // T = transpose(T_tmp);
}

return b;
}