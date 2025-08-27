#include <iostream>
#include <random>
#include <unistd.h> 

#include <immintrin.h>  // For _mm_popcnt_u64 
#include <bitset>       // For std::popcount 

#include "Tools/tools.hpp"

using namespace std;
using namespace NTL;
/* **************************************************************************** */
/*                              MISCELLANEOUS                                   */
/* **************************************************************************** */

/* Write the positions of non-zero bits in the pos vector*/
void extract_position(uint64_t* poly_in, int* pos ,int length){
    int i = 0;
    for (int word_idx = 0; word_idx < length; ++word_idx) {
        uint64_t word = poly_in[word_idx];
        int base_idx = word_idx * 64;  

        while (word) {  
            int bit_pos = __builtin_ctzll(word); 
            pos[i]= (base_idx + bit_pos); 
            word &= (word - 1);  
            i++;
        }
    }
}


/* Stacks an int vector containing binary values into a 64 bit-chunks*/
void bit_to_uint(uint64_t *output, int *input, int length)
{
    int size_int = (length + 63) / 64;
    for (int i = 0; i < size_int; i++)
    {
        output[i] = 0;
    }

    for (int i = 0; i < length; i++)
    {
        output[i / 64] |= ((uint64_t)input[i] << (i % 64));
    }
}

/* return integer r with "len" bits such that 2 is primitive modulo r */
int random_suitable_integer(const int len)
{
    int res = 0;
    FLINT_TEST_INIT(state);
    fq_nmod_t a;
    fq_nmod_ctx_t ctx;
    while (1)
    {
        res = n_randprime(state, len, 0);
        fq_nmod_ctx_init_ui(ctx, res, 1, "x");
        fq_nmod_init(a, ctx);
        fq_nmod_set_ui(a, 2, ctx);
        if (fq_nmod_is_primitive(a, ctx))
            break;
    }
    flint_randclear(state);
    fq_nmod_clear(a, ctx);
    fq_nmod_ctx_clear(ctx);
    return res;
}

/**
 * Computes the threshold for the BGF algo.
 * Use hard-coded values from Bike spec.
 * TODO: find a better definition
 */
int compute_threshold(const int w, const int ind, const int r)
{
    if (r == 12323)
    {
        return std::max(0.0069722 * w + 13.530, 36.);
    }
    else if (r == 24659)
    {
        return std::max(0.005265 * w + 15.2588, 52.);
    }
    else
    {
        return std::max(0.00402312 * w + 17.8785, 69.);
    }
}

/* **************************************************************************** */
/*                      RANDOM INDICES / PERMUTATIONS                           */
/* **************************************************************************** */

/* naive fisher-yates algorithm */
void fisher_yates(int *perm, const int n)
{
    int j;

    std::random_device rd;
    std::mt19937 rand_gen(rd());

    for (int i = 0; i < n; ++i)
    {
        std::uniform_int_distribution<std::mt19937::result_type> dis(0, i);
        j = dis(rand_gen);
        if (i != j)
        {
            perm[i] = perm[j];
        }
        perm[j] = i;
    }
}

slong *
random_indices(const slong len, flint_rand_t state)
{
    slong *perm = _perm_init(len);
    _perm_randtest(perm, len, state);
    return (perm);
}

/* check repetitions in a naive way */
int int_check_repeat(const int *a, const int len)
{
    for (int i = 0; i < len; i++)
    {
        for (int j = i + 1; j < len; j++)
        {
            if (a[i] == a[j])
                return 1;
        }
    }
    return 0;
}

/**
 * Computes t random and pairwise distinct indexes in [1,n] as in Classic
 * McEliece.
 * Draws random τ > t elts in [1, N[ and takes the t first ones that are in [1,n].
 * Returns the number of elts in [1,n] effectively computed.
 * Typically in CM, N = 2^m.
 */
int cm_random_indices(int *res, const int n, const int t, const int N, const int tau)
{
    int ind = 0;
    int count_loop = 0;
    int a;

    // rand with mersenne
    std::random_device rd;
    std::mt19937 rand_gen(rd());
    std::uniform_int_distribution<std::mt19937::result_type> dis(0, N - 1);

    while (ind < t && count_loop < tau)
    {
        a = dis(rand_gen);
        if (a < n)
        {
            res[ind] = a;
            ind++;
        }
        count_loop++;
    }
    return ind;
}

void random_bits(int *res, const int len)
{
    int ind = 0;
    int count_loop = 0;
    int a;

    // rand with mersenne
    std::random_device rd;
    std::mt19937 rand_gen(rd());
    std::uniform_int_distribution<std::mt19937::result_type> dis(0, 1);

    for (int i = 0; i < len; ++i)
    {
        res[i] = dis(rand_gen);
    }
}

void random_bytes(int *res, const int len)
{
    int ind = 0;
    int count_loop = 0;
    int a;

    // rand with mersenne
    std::random_device rd;
    std::mt19937 rand_gen(rd());
    std::uniform_int_distribution<std::mt19937::result_type> dis(0, 255);

    for (int i = 0; i < len; ++i)
    {
        res[i] = dis(rand_gen);
    }
}

/* **************************************************************************** */
/*                              ERROR GENERATION                                */
/* **************************************************************************** */

/**
 * Computes e as in Classic McEliece, more or less
 */
void cm_gen_e(int *e, const int n, const int t, const int N, const int tau)
{
    int len, repet = 1;
    int inds[t];
    while (repet)
    {
        len = cm_random_indices(inds, n, t, N, tau); /* compute random set of indices */
        if (len == t)
        { /* check if we obtained enough indices */
            repet = int_check_repeat(inds, t);
        }
    }
    for (int i = 0; i < t; i++)
    {
        e[inds[i]] = 1;
    }
}

/**
 * Generates fixed weight vector as in Bike, i.e. using
 * Fisher-Yates
 */
void bike_gen_e(int *e, const int len, const int t)
{
    int perm[len];
    for (int i = 0; i < len; ++i)
    {
        perm[i] = 0;
    }
    fisher_yates(perm, len);
    for (int i = 0; i < t; ++i)
    {
        e[perm[i]] = 1;
    }
}

/**
 * Generates fixed weight vector as in HQC, i.e. using
 * Fisher-Yates
 */
void hqc_gen_e(int *e, const int len, const int t)
{
    int perm[len];
    for (int i = 0; i < len; ++i)
    {
        e[i] = 0;
        perm[i] = 0;
    }
    fisher_yates(perm, len);
    for (int i = 0; i < t; ++i)
    {
        e[perm[i]] = 1;
    }
}

/* **************************************************************************** */
/*                              FINITE FIELDS                                   */
/* **************************************************************************** */
/**
 * Generates a random vector whose coeffs are in a finite field and all distincts
 * WORK ONLY for *BINARY* fields
 */
void vec_rand_distinct_2(vec_GF2E &res, const int len,
                              flint_rand_t state)
{
    slong d = 1 << GF2E::degree();
    slong *inds = random_indices(d, state); /* here is the binary stuff */

    int inds2[len];
    for (int i = 0; i < len; i++)
    {
        inds2[i] = (int)inds[i];
    }
    GF2EX tmp = int_vec_to_GF2EX(inds2, len);
    res = VectorCopy(tmp, len);
}


/* **************************************************************************** */
/*                           POLYNOMIAL MANIPULATION                            */
/* **************************************************************************** */
/**
 * Check if f has a root in alpha.
 */
int gf2ex_eval_zero(const GF2EX &f, const vec_GF2E &alpha, const int len)
{
    int b = 0;
    vec_GF2E result = eval(f, alpha);
    for (int i = 0; i < len; i++)
    {
        b = IsZero(result[i]);
        if (b)
            break;
    }
    return (b);
}

/**
 * Compute an irreducible polynomial of degree "deg" as in CM ie so that it
 * generates a Goppa code together with the roots alpha Γ(alpha, res)
 */

void cm_poly_irr_pol(GF2EX &res, const int deg, const vec_GF2E &alpha, const int len)
{
    GF2EX g = BuildIrred_GF2EX(deg);
    res = BuildRandomIrred(g);
    while (gf2ex_eval_zero(res, alpha, len))
    {
        printf("a");
        res = BuildRandomIrred(g);
    }
}

/**
 * Computes early abort extended gcd of a and b in finite field
 */
void XGCD_abort(GF2EX& d, GF2EX& s, GF2EX& t, const GF2EX& a, const GF2EX& b, long k) {
    GF2EX u0, u1, v0, v1, d0, d1, q, r, tmp;
    // Initialisation
    clear(u1); set(u0); // u0 = 1, u1 = 0
    clear(v0); set(v1); // v0 = 0, v1 = 1
    d0 = a;
    d1 = b;
    int cpt = 0;

    while (deg(d0) > k && !IsZero(d1)) {
        DivRem(q, r, d0, d1); // q = d0 / d1, r = d0 % d1

        cpt++;
        d0 = d1;
        d1 = r;

        //tmp = u1 * q;
        tmp = u0 - u1 * q;// u0 - tmp;
        u0 = u1;
        u1 = tmp;

        //tmp = v1 * q;
        tmp = v0 - v1 * q; //v0 - tmp;
        v0 = v1;
        v1 = tmp;
    }

    // Results
    d = d0;
    s = u0;
    t = v0;
}


/* **************************************************************************** */
/*                                 MATRICES                                     */
/* **************************************************************************** */

/**
 * Computes the reduced row echelon form from the row echelon form
 */
void rref(mat_GF2& A) {
    long m = A.NumRows();
    long n = A.NumCols();
    
    long lead = 0;
    for (long r = 0; r < m; r++) {
        if (lead >= n) return;
        long i = r;
        while (A[i][lead] == 0) { 
            i++;
            if (i == m) {
                i = r;
                lead++;
                if (lead == n) return;
            }
        }
        swap(A[i], A[r]);
        for (long j = 0; j < m; j++) {
            if (j != r && A[j][lead] == 1) A[j] += A[r]; 
        }
        lead++;
    }
}

/**
 * Expansion of a matrix over F_q^m to a matrix over FF_q
 * NEED TO OBTAIN OUTPUT AS A MATRIX IN FF_q
 */


/* **************************************************************************** */
/*                                   BIKE                                       */
/* **************************************************************************** */

/**
 * Counter function as in Bike documentation.
 * Computes the number of agreeing bits in vector v and jth column vector of H.
 */
int ctrv2(const GF2X &sp, const int *pos, const int j, const int weight, const int r)
{

    int count = 0;
    int k = 0;
    const int degree = deg(sp);

    for (int i = 0; i < weight; ++i)
    {
        k = (j + pos[i]) % r;
        // printf("%d ", k);

        if (k < degree && rep(sp[k]))
        {
            count++;
        }
    }
    /* printf("\n"); */
    return count;
}


/** 
 * Optimized version of ctrv2
 */
int ctrv2_v2(unsigned char *spBytes, const int *pos, const int j, const int weight, const int r)
{

    int count = 0;
    int k = 0;

    for (int i = 0; i < weight; ++i)
    {
        k = j + pos[i]; // Could be pre-computed
        if (k >= r) k -= r;  // Since k< 2r-2

        count += (spBytes[k >> 3] >> (k & 7)) & 1;

    }
    /* printf("\n"); */
    return count;
}


/**
 * Computes black and gray positions as in Bike specification doc.
 */
void BFIterv2(GF2X &e0, GF2X &e1, int *black, int *gray, const GF2X &sp,
                   const int *pos0, const int *pos1, const int weight, const int r, const int T,
                   const int tau)
{

    int i, j, count;

    int n = 2 * r;

    // Not degree, size
    int deg0 = NumBits(e0);  // deg(e0) + 1;
    int deg1 = NumBits(e1);  // deg(e1) + 1;
    int degsp = NumBits(sp); // deg(sp) + 1;

    // Conversion to bytes for faster counting (ctrv2)
    long numB = NumBytes(sp);
    unsigned char spBytes[numB];
    BytesFromGF2X(spBytes, sp, numB);

    // GF2 one; one = 1;
    GF2 tmp;

    /* initialisation black and gray vectors */
    for (i = 0; i < n; ++i)
    {
        black[i] = 0;
        gray[i] = 0;
    }
    
    if (degsp == 0) {
    std::fill(black, black + n, 0);
    std::fill(gray, gray + n, 0);
} else {
    // Première boucle pour j < r
    for (j = 0; j < r; ++j) {
        count = ctrv2_v2(spBytes, pos0, j, weight, r);

        if (count >= T) {
            if (j < deg0) {
                SetCoeff(e0, j, e0[j] + 1);
            } else {
                SetCoeff(e0, j, 1);
            }
            black[j] = 1;
        } else if (count >= T - tau) {
            gray[j] = 1;
        }
    }

    // Deuxième boucle pour j >= r
    for (j = r; j < n; ++j) {
        count = ctrv2_v2(spBytes, pos1, j - r, weight, r);

        if (count >= T) {
            if ((j - r) < deg1) {
                SetCoeff(e1, j - r, e1[j - r] + 1);
            } else {
                SetCoeff(e1, j - r, 1);
            }
            black[j] = 1;
        } else if (count >= T - tau) {
            gray[j] = 1;
        }
    }
}

    e0.normalize();
    e1.normalize();
}

/**
 * Modify e using black or gray positions as in Bike specification doc.
 */
void BFMaskedIterv2(GF2X &e0, GF2X &e1, const GF2X &sp, const int *pos0,
                         const int *pos1, const int weight, const int r, const int T,
                         const int *mask)
{
    int j, count;

    int n = 2 * r;

    int deg0 = NumBits(e0); // = deg(e0) + 1;
    int deg1 = NumBits(e1); // = deg(e1) + 1;

    GF2 tmp, tmp1;

    // Conversion into bytes (speeds up ctrv2)
    long numB = NumBytes(sp);
    unsigned char spBytes[numB];
    BytesFromGF2X(spBytes, sp, numB);

    for (j = 0; j < n; ++j)
    {
        if (j < r)
        {
            count = ctrv2_v2(spBytes, pos0, j, weight, r);
        }
        else
        {
            count = ctrv2_v2(spBytes, pos1, j - r, weight, r);
        }
        if (count >= T)
        {
            tmp = mask[j];

            if (j < r)
            {
                if (j < deg0)// && deg0 > 0)
                    SetCoeff(e0, j, e0[j]+tmp);
                else
                    SetCoeff(e0, j, tmp);

                /****** */
                // GF2 tmp_e0;
                // if (j < deg0)// && deg0 > 0)
                //     tmp_e0 = e0[j];
                // else
                //     tmp_e0 = 0;
                // tmp1 = tmp_e0 + tmp;
                // SetCoeff(e0, j, tmp1);
            }
            else
            {
                if ((j - r) < deg1)// && deg1 > 0)
                    SetCoeff(e1, j - r, e1[j-r]+tmp);
                else
                    SetCoeff(e1, j - r, tmp);

                // /****** */
                // GF2 tmp_e1;
                // if ((j - r) < deg1)// && deg1 > 0)
                //     tmp_e1 = e1[j - r];
                // else
                //     tmp_e1 = 0;
                // tmp1 = tmp_e1 + tmp;
                // SetCoeff(e1, j - r, tmp1);
            }
        }
    }
    e0.normalize();
    e1.normalize();
}

/**
 * Short: Computes CM parameters as in the specification document.
 */
void CM_params(int &m, int &n, int &t, GF2X &P, const int level)
{
    if (level == 1)
    { /* mceliece348864 */
        m = 12;
        n = 3488;
        t = 64;
        // Context for GF2E
        SetCoeff(P, 12); // X^12
        SetCoeff(P, 7);  // X^7
        SetCoeff(P, 6);  // X^6
        SetCoeff(P, 5);  // X^5
        SetCoeff(P, 3);  // X^3
        SetCoeff(P, 1);  // X
        SetCoeff(P, 0);  // 1
    }
    else if (level == 3)
    { /* mceliece6688128 */
        m = 13;
        n = 6688;
        t = 128;
        // Context for GF2E
        SetCoeff(P, 13); // X^13
        SetCoeff(P, 4);  // X^4
        SetCoeff(P, 3);  // X^3
        SetCoeff(P, 1);  // X^1
        SetCoeff(P, 0);  // X^0
    }
    else if (level == 5)
    { /* mceliece8192128 */
        m = 13;
        n = 8192;
        t = 128;
        // Context for GF2E
        SetCoeff(P, 13); // X^13
        SetCoeff(P, 4);  // X^4
        SetCoeff(P, 3);  // X^3
        SetCoeff(P, 1);  // X^1
        SetCoeff(P, 0);  // X^0
    }
}

/**
 * Short: Computes BIKE parameters as in the specification document.
 */
void Bike_params(int &r, int &weight, int &error_weight, const int level)
{
    if (level == 1)
    {
        r = 12323;
        weight = 142;
        error_weight = 134;
    }
    else if (level == 3)
    {
        r = 24659;
        weight = 206;
        error_weight = 199;
    }
    else if (level == 5)
    {
        r = 40973;
        weight = 274;
        error_weight = 264;
    }
}

/**
 * Short: Computes BGF parameters NbIter and tau as in specification document.
 */
void BGF_params(int &NbIter, int &tau, const int level)
{
    if (level == 1)
    {
        NbIter = 5;
        tau = 3;
    }
    else if (level == 3)
    {
        NbIter = 5;
        tau = 3;
    }
    else if (level == 5)
    {
        NbIter = 5;
        tau = 3;
    }
}

/**
 * Short: Computes HQC parameters `almost` as in the specification document.
 */
void HQC_params(int &k, int &n1, int &n2, int &r, int &n, int &w, int &we, int &wr,
                const int level)
{
    if (level == 1)
    {
        k = 16;
        n1 = 46;
        n2 = 128;
        r = 3;
        n = 17669;
        w = 66;
        we = 75;
        wr = 75;
    }
    else if (level == 3)
    {
        k = 24;
        n1 = 56;
        n2 = 128;
        r = 5;
        n = 35851;
        w = 100;
        we = 114;
        wr = 114;
    }
    else if (level == 5)
    {
        k = 32;
        n1 = 90;
        n2 = 128;
        r = 5;
        n = 57637;
        w = 131;
        we = 149;
        wr = 149;
    }
}


/* **************************************************************************** */
/*                                 MATRICES                                     */
/* **************************************************************************** */

/**
 * Horizontal concatenation C = [M1 M2]
 */ 
mat_GF2 concat_hor_mat_GF2(const mat_GF2& M1, const mat_GF2& M2){
    if (M1.NumRows() != M2.NumRows()) {
        throw std::invalid_argument("Unequal row number! Please check matrix dimensions");
    }

    long rows = M1.NumRows();
    long colsM1 = M1.NumCols();
    long colsM2 = M2.NumCols();

    // Output matrix
    mat_GF2 C;
    C.SetDims(rows, colsM1 + colsM2);

    // FIlling left part
    for (long i = 0; i < rows; i++) {
        for (long j = 0; j < colsM1; j++) {
            C[i][j] = M1[i][j];
        }
    }

    // Filling right section
    for (long i = 0; i < rows; i++) {
        for (long j = 0; j < colsM2; j++) {
            C[i][colsM1 + j] = M2[i][j];
        }
    }

    return C;
}

/**
 * Expansion of a matrix over F_q^m to a matrix over FF_q
 * mat_GF2E is used since mat_GF2X does not exist
 */
mat_GF2 matrix_expand(mat_GF2E& M){
    mat_GF2 exp_M;
    long row = M.NumRows();
    long col = M.NumCols();
    long m = GF2E::degree(); // Extension degree 
    exp_M.SetDims(row*m, col);
    clear(exp_M);

    for(int i = 0; i < row ; i++){
        for(int j = 0; j < col; j++){
            GF2X tmp = rep(M[i][j]); // GF2X representation of GF2E 
            int length = NumBits(tmp);
            for (int n = 0; n < m; n++) {  
                if (n < length) {  
                    exp_M[i * m + n][j] = coeff(tmp, n);  
                } else {
                    exp_M[i * m + n][j] = 0;  
                }
            }
        }
    }
    return exp_M;
}



/**
 * GF2X to int representation 
 * Mainly used for binary vector in McEliece?
 * rep(poly[0]) for binary vector
 */
int gf2x_to_int(const GF2X& poly) {
    int result = 0;
    for (long i = 0; i <= deg(poly); ++i) {
        if (IsOne(coeff(poly, i))) {
            result |= (1 << i); 
        }
    }
    return result;
}


/**
 * GF2 to int array conversion
 */
void vec_gf2_to_int(const Vec<GF2>& in_vec, int* result) {
    long size = in_vec.length();
    for (long i = 0; i < size; i++) {
        result[i] = rep(in_vec[i]); 
    }
}

/**
 * Used for binary converison...  only for binary values (maybe not optimal)
 */
void vec_gf2e_to_int(const Vec<GF2E>& in_vec, int* result) {
    long size = in_vec.length();
    for (long i = 0; i < size; i++) {
        result[i] = IsOne(in_vec[i]); 
    }
}

/**
 * Various conversions from int vec to polynomials
 */
GF2EX uint8_vec_to_GF2EX(int* in_vec, long len){
    GF2EX out_vec(INIT_SIZE,len); // Allocates memory
    out_vec.SetLength(len);       // Sets the number of coefficient (sets coeff >=len to 0)

    for(int i = 0; i< len; i++){
        GF2X tmp_coeff;
        tmp_coeff.SetMaxLength(8);
        tmp_coeff.SetLength(8);
        for(int j = 0; j<8;j++){
            tmp_coeff[j] = (in_vec[i] >> j) & 1;
        }
        tmp_coeff.normalize();  // Must be called after setting the coefficients manually
        out_vec[i] = to_GF2E(tmp_coeff);
    }
    return out_vec;
}

GF2EX uint8_vec_to_GF2EX_V2(unsigned char* in_vec, long len){
    GF2EX poly(INIT_SIZE,len);
    GF2X coeff;
    for(int i = 0; i<len;i++){
        poly[i] = to_GF2E(GF2XFromBytes(in_vec, 1));
    }
    return poly;
}

GF2EX int_vec_to_GF2EX(int* in_vec, long len){
    int ncoeff = GF2E::degree() + 1;
    GF2EX out_vec;
    out_vec.SetMaxLength(len);
    out_vec.SetLength(len);

    for(int i = 0; i< len; i++){
        GF2X tmp_coeff;
        tmp_coeff.SetMaxLength(ncoeff);
        tmp_coeff.SetLength(ncoeff);
        for(int j = 0; j<ncoeff;j++){
            tmp_coeff[j] = (in_vec[i] >> j) & 1;
        }
        tmp_coeff.normalize();  // Must be called after setting the coefficients manually
        out_vec[i] = to_GF2E(tmp_coeff);
    }
    return out_vec;
}


GF2X int_vec_to_GF2X(int* in_vec,long len){
    GF2X out_vec;
    out_vec.SetMaxLength(len);
    out_vec.SetLength(len);

    for(int i = 0; i< len; i++){
        out_vec[i] = in_vec[i];
    }
    out_vec.normalize();
    return out_vec;
}

/**
 * Sets d to "cyclic polynomial generator", i.e. X^d - 1 (over FF_2)
 */
GF2X gf2x_set_cyclic(const int d) {
    GF2X P;
    P.SetMaxLength(d);
    P.SetLength(d);
    P[d] = 1;
    P[0] = 1;
    return P;
}


/**
 * Converts each coeff of the GF2EX polynomial into an integer array 
 */ 
void GF2EX_to_int(const GF2EX& in_pol, int* result, int len){
    int m = GF2E::degree();
    int d = deg(in_pol);
    for(int i = 0; i <= d; i++){
        GF2X tmp = rep(in_pol[i]);
        int tmp_deg = (deg(tmp)+1);
        int tmp_val = 0;
        for(int k = 0; k < tmp_deg ; k++){
            tmp_val += (IsOne(tmp[k])) << k;
        }
        result[i]  = tmp_val;
    }
    
    // Needed for polynomials with highest coeff to 0: 0X^n + 1X^(n-1) + ... + 1
    for(int i = d+1; i < len ; i ++){
        result[i] = 0;
    }
}


/**
 * Converts each coeff of the GF2EX polynomial into an integer array
 */
void GF2X_to_int(const GF2X& in_pol, int* result, int len){
    int d = deg(in_pol);
    for(int i = 0; i <= len; i++){
        result[i] = rep(in_pol[i]);
    }
    
    // Needed for polynomials with highest coeff to 0: 0X^n + 1X^(n-1) + ... + 1
    for(int i = d+1; i < len ; i ++){
        result[i] = 0;
    }
}