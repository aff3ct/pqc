#include <iostream>
#include "Tools/codes.hpp"

#include <NTL/GF2X.h>
#include <NTL/GF2EX.h>
#include <vector>

// #include <chrono>

/* **************************************************************************** */
/*                                   DECODERS                                   */
/* **************************************************************************** */

using namespace std;
/**
 * Decoding algorithm for Reed-Solomon codes
 * Decode received message c wrt to RS_k(alpha)
 */
int RS_decoder(GF2EX &res, const Vec<GF2E> &c, const Vec<GF2E> &alpha, int n, int k)
{
    GF2EX A,P, tmp, u, v, D;

    A.SetMaxLength(2 * n + 1);  /* Oversized to fit the derivative */ 
    A.SetLength(1);             /* A has a starting length of 1 (A = 1) */ 
    A[0] = GF2E(1);             /* Setting the first coefficient to 1 -> A = 1 */ 

    tmp.SetMaxLength(2);         
    tmp.SetLength(2);           
    tmp[1] = GF2E(1);           /* setting tmp = X */     

    P.SetMaxLength(n);
    P.SetLength(n);

    for (int i = 0; i < n; ++i)
    {
        tmp[0] = -alpha[i];     /* tmp = X - alpha[i] */
        A *= tmp;   
    }
    P = interpolate(alpha, c);

    int d = k + (n - k + 1) / 2 - 1; /* computes k + ceil((n-k)/2) -1 */
    XGCD_abort(D, u, v, A, P, d);

    int b = divide(res, D, v);

    return b && (deg(res) <= k);
}

/**
 * Generalized RS decoder
 */
int GRS_decoder(Vec<GF2E> &res, const Vec<GF2E> &c, const Vec<GF2E> &alpha, const Vec<GF2E> &beta, int n, int k)
{
    vec_GF2E c1(INIT_SIZE, n);
    c1.SetLength(n);

    GF2EX tmp_rs_res(INIT_SIZE, n); /* Temporary variable to store the RS result */
    tmp_rs_res.SetLength(1);

    for (int i = 0; i < n; i++)
    {
        c1[i] = c[i] * inv(beta[i]);
    }

    int b = RS_decoder(tmp_rs_res, c1, alpha, n, k);

    if (b != 0)
    {
        clear(res);
        vec_GF2E tmp_eval = eval(tmp_rs_res,alpha);

        for (int ii = 0; ii < n; ii++)
        {
            if(!IsZero(tmp_eval[ii])){
                res[ii] = beta[ii] * tmp_eval[ii];
            }
        }
    }

    return b;
}

/* Decoding algorithm for *BINARY* Goppa codes */
int Goppa_decoder_bin(Vec<GF2E> &res, const Vec<GF2E> &c, const Vec<GF2E> &alpha,
                           const GF2EX &g, const int n, const int k)
{

    GF2EX A(INIT_SIZE, n + 2);
    A.SetLength(1); // A.SetLength(n+2);
    A[0] = GF2E(1); /* A = 1 */

    GF2EX tmp_poly(INIT_SIZE, 2);
    tmp_poly.SetLength(2);
    tmp_poly[1] = GF2E(1); /* 1*x */

    for (int i = 0; i < n; ++i)
    {
        tmp_poly[0] = -alpha[i]; /* tmp_poly = X + 1*alpha[i] */ 
        A *= tmp_poly;
    }
    diff(A, A);

    vec_GF2E beta(INIT_SIZE, n);
    beta.SetLength(n);

    vec_GF2E eval_g = eval(g, alpha);
    vec_GF2E eval_A = eval(A, alpha);
    for (int i = 0; i < n; i++)
    {
        beta[i] = sqr(eval_g[i]) / eval_A[i];
    }

    int b = GRS_decoder(res, c, alpha, beta, n, k);

    return b;
}

/* **************************************************************************** */
/*                                   MISC                                       */
/* **************************************************************************** */

/* Computes the parity check matrix of a Goppa code */
mat_GF2E Goppa_parity_check(const vec_GF2E &alpha, const GF2EX &g, int row, int col)
{
    mat_GF2E M;
    M.SetDims(row, col);

    vec_GF2E eval_g = eval(g, alpha);
    //vec_GF2E beta(INIT_SIZE, col);

    for (int i = 0; i < col; i++)
    {
        //beta[i] = inv(eval_g[i]);
        //M[0][i] = beta[i];
        M[0][i] = inv(eval_g[i]);
    }

    /* computes the following lines by multiplying by alpha_i */
    for (int i = 1; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            M[i][j] = M[i - 1][j] * alpha[j];
        }
    }

    return M;
}

/* **************************************************************************** */
/*                                   ENCODERS                                   */
/* **************************************************************************** */
/**
 * Encoding for duplicated Reed-Muller codes of parameters (1, m)
 * IN: alpha is a vector over FF_2
 * OUT: vector of 1 + sum_{0 < i < m} alpha_i *  X_i evaluated in all elements of FF_2^n
 * duplicated r times
 */
void RM_encoding_duplicated(vec_GF2 &res, const mat_GF2 &G, const vec_GF2 &alpha, int m, int r)
{
    int len = 1 << m; // Code length

    /* RM encoding with generator matrix */
    vec_GF2 tmp_res(INIT_SIZE, len);

    mul(tmp_res, alpha, G);
    // tmp_res.SetLength(len); 

    res = tmp_res;
   
    for (int i = 0; i < (r-1); ++i)
    {
        res.append(tmp_res);
    }

}

/**
 * Encodes with concatenated Reed-Solomon and Reed-Muller codes.
 * Here `f` is a polynomial in F_2^m[X], `alpha` is a vector of elements in F_2^m,
 * `n` is the length of the RS code, we are using RM(1, m-1) duplicated `r` times
 * and `ctx` / `ctx_q` is F_2^m / F_2.
 */
void RS_RM_concatenated_encoding(vec_GF2 &out, const GF2EX &f,
                                 const vec_GF2E &alpha, const int n, const int r, int m)
{
    int i, j;
    int len = (1 << (m - 1)); /* length of RM(1, m-1) code */
    int _m = m-1;

    vec_GF2E rs_res = eval(f, alpha); // Reed-Solomon encoding

    // Generator matrix
    mat_GF2 G;
    G.SetDims(m, len);

    for (long i = 0; i < len; i++) {
        long tmp_i = i;
        for (long j = 0; tmp_i > 0; j++, tmp_i >>= 1) { 
            if (tmp_i & 1) { 
                G[j + 1][i] = 1; 
            }
        }
        G[0][i] = 1;
    }

    /* Temporary variables */
    vec_GF2 tmp_alpha(INIT_SIZE,m);
    tmp_alpha.SetLength(m);
    vec_GF2 tmp_res(INIT_SIZE, r*len);
    tmp_res.SetLength(r*len);

    /* First iteration of the decoding loop */
    tmp_alpha = conv<vec_GF2>(rep(rs_res[0]));
    tmp_alpha.SetLength(m);

    RM_encoding_duplicated(tmp_res, G, tmp_alpha, m - 1, r);
    out = tmp_res;
    clear(tmp_res);

    /* Rest of the loop */
    for (i = 1; i < n; ++i)
    {
        tmp_alpha = conv<vec_GF2>(rep(rs_res[i]));
        tmp_alpha.SetLength(m); // length ajustement to m for matrix-vector mul

        RM_encoding_duplicated(tmp_res, G, tmp_alpha, m - 1, r);
        
        out.append(tmp_res);
        clear(tmp_res);
    }
}

/**
 * Encoding as in CM
 * Compute the syndrome associated s = H * e associated to e.
 */
void CM_encoding(vec_GF2 &res, const mat_GF2 &T, const vec_GF2 &e)
{

    int n = T.NumRows();
    int m = T.NumCols();
    
    vec_GF2 e0 = VectorCopy(e,n);
    vec_GF2 e1(INIT_SIZE,m);
    e1.SetLength(m);

    e1 = VectorCopy(shift(e,-n),m);

    mul(res, T,e1);
    res += e0; 

}

/**
 * Encoding as in Bike
 * Compute s = e0 + e1 * h modulo X^r - 1
*/
void Bike_encoding(GF2X &res, const int *e, const GF2X &h, const int r)
{
    /* init P to X^r -1 */
    GF2X P = gf2x_set_cyclic(r);
    GF2X e0(INIT_SIZE, r);
    GF2X e1(INIT_SIZE, r);

    /* Set e0 to e[0..r-1] and e1 to e[r..2r-1] */
    for (int i = 0; i < r; ++i)
    {
        SetCoeff(e0, i, e[i]);
        SetCoeff(e1, i, e[i + r]);
    }

    /* Set coeff > degree to 0 */
    e0.normalize();
    e1.normalize();

    res = MulMod(e1, h, P) + e0;
}

/**
 * Encoding as in HQC
 */
void HQC_encoding(GF2X &res1, GF2X &res2, const GF2EX &m, const GF2X &h,
                       const GF2X &s, const vec_GF2E &alpha, const int n, const int k,
                       const int n1, const int r, const int we, const int wr,
                       int _m, GF2X &P)
{
    int len = (1 << (_m - 1)); /* length of RM(1, m-1) code */
    int out_size = r * len * n1;

    int v1[n];
    int v2[n];
    int v3[n];

    /* computes r1, r2 with fixed weight `wr` */
    hqc_gen_e(v1, n, wr); /* r1 */
    hqc_gen_e(v2, n, wr); /* r2 */
    hqc_gen_e(v3, n, wr); /* tmp pol */

    GF2X r1(INIT_SIZE,n);
    r1.SetLength(n);

    GF2X r2(INIT_SIZE,n);
    r2.SetLength(n);

    GF2X tmp_pol(INIT_SIZE,n);
    tmp_pol.SetLength(n);

    for (int i = 0; i < n; i++)
    {
        r1[i] = v1[i];
        r2[i] = v2[i];
        tmp_pol = v3[i];
    }
    r1.normalize();
    r2.normalize();
    tmp_pol.normalize();

    vec_GF2 out(INIT_SIZE,out_size);
    out.SetLength(out_size);
    RS_RM_concatenated_encoding(out, m, alpha, n1, r, _m);

    GF2X pol_out2 = conv<GF2X>(out);

    res1 = r1 + MulMod(h, r2, P);
    res2 = MulMod(s, r2, P) + pol_out2 + tmp_pol;

    res2 = trunc(res2, out_size);
}

/* **************************************************************************** */
/*                                   DECODERS                                   */
/* **************************************************************************** */
/* decode syndrome as in CM */
int CM_syndrome_decoding_bin(Vec<GF2E> &res, const Vec<GF2E> &s, const Vec<GF2E> &alpha, const GF2EX &g,
                                  const int len, const int t)
{
    int d = deg(g);
    int m = GF2E::degree();

    vec_GF2E c = s;

    int b = Goppa_decoder_bin(res, c, alpha, g, len, len - 2 * t);

    return b;
}


/**
 * Decoding as in Bike : use the Black-Gray-Flip (BGF) algorithm
*/
int Bike_decoding_v2(vec_GF2 &res, const GF2X &s, const GF2X &h0,
                          const GF2X &h1, const int r, const int weight,
                          const int NbIter, const int tau)
{
    int i, T, w;
    int n = 2 * r;
    int d = (weight + 1) / 2 + 1;
    int b = 0;

    /* P = X^r - 1 */
    GF2X P = gf2x_set_cyclic(r);

    /* Size variables */
    int size0 = NumBits(h0); //deg(h0) + 1;
    int size1 = NumBits(h1); //deg(h1) + 1;

    int size_int_0 = (size0+63)/64;
    int size_byte_0 = NumBytes(h0); 

    int size_int_1 = (size1+63)/64;
    int size_byte_1 = NumBytes(h1); 
    uint8_t h0_byte[size_byte_0];
    uint8_t h1_byte[size_byte_1];

    /* Conversion into 64-bits array for faster computation using extract_position */ 
    // Init??
    uint64_t h0_int[size_int_0];
    uint64_t h1_int[size_int_1];

    BytesFromGF2X(h0_byte, h0, size_byte_0);
    memcpy(h0_int, h0_byte, size_byte_0); 
    BytesFromGF2X(h1_byte, h1, size_byte_1);
    memcpy(h1_int, h1_byte, size_byte_1); 

    /* position of non-zero coefficients of h0, h1 (secret key) */
    int pos0[weight];
    int pos1[weight];
    /* Extract the position of 1 bit within h0 and h1 */
    extract_position(h0_int,pos0,size_int_0);
    extract_position(h1_int,pos1,size_int_1);

    GF2X e0(INIT_SIZE, r); 
    GF2X e1(INIT_SIZE, r); 
    GF2X tmp(INIT_SIZE, r); 
    GF2X synd; // Syndrome

    int black[n];
    int gray[n];

    for (i = 0; i < NbIter; ++i)
    {
        /* Computes s' = s + e*h */
        tmp = MulMod(e0, h0, P) + MulMod(e1, h1, P) + s;

        /* Hamming weight of s' + threshold */
        w = NTL::weight(tmp);
        T = compute_threshold(w, i, r);

        BFIterv2(e0, e1, black, gray, tmp, pos0, pos1, weight, r, T, tau);

        if (i == 0)
        {
            BFMaskedIterv2(e0, e1, tmp, pos0, pos1, weight, r, d, black);
            BFMaskedIterv2(e0, e1, tmp, pos0, pos1, weight, r, d, gray);
        }
    }

    /* Computes the syndrome */
    synd = MulMod(e0, h0, P) + MulMod(e1, h1, P);

    if (synd == s)
    {
        vec_GF2 tmp0 = conv<vec_GF2>(e0);
        tmp0.SetMaxLength(r); tmp0.SetLength(r);
        vec_GF2 tmp1 = conv<vec_GF2>(e1);
        tmp1.SetMaxLength(r); tmp1.SetLength(r);
        append(tmp0,tmp1);

        res = tmp0;
        
        b = 1;
    }

    return b;
}

/**
 * Reed-Muller decoding using Hadamard transform
 * Adapted from Ref. Implementation ~ https://pqc-hqc.org/implementation.html
 */
void RM_decoding_duplicated(GF2X &res, const vec_GF2 &c, const int m, const int r)
{
    int len = 1 << m;
    std::vector<int> c_prime(len, 0);
    std::vector<int> d(len, 0);

    for (int i = 0; i < r; i++)
    {
        int idx = i * len;
        for (int j = 0; j < len; j++)
        {
            c_prime[j] += rep(c[idx + j]);
        }
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < len / 2; j++)
        {
            d[j] = c_prime[2 * j] + c_prime[2 * j + 1];
            d[j + len / 2] = c_prime[2 * j] - c_prime[2 * j + 1];
        }
        std::swap(d, c_prime);
    }

    c_prime[0] -= (len / 2) * r;

    int32_t peak_abs_value = 0, peak_value = 0, peak_pos = 0;

    for (int i = 0; i < len; i++)
    {
        int32_t t = c_prime[i];
        int32_t absolute = std::abs(t);

        if (absolute > peak_abs_value)
        {
            peak_value = t;
            peak_pos = i;
            peak_abs_value = absolute;
        }
    }

    peak_pos |= len * (peak_value > 0);
    uint8_t lower_8_bits = peak_pos & 0xFF;
    lower_8_bits = (lower_8_bits << 1 ) | ((lower_8_bits>>7)&1);

    GF2XFromBytes(res,&lower_8_bits,1);

}

/**
 * Encodes with concatenated Reed-Solomon and Reed-Muller codes.
 * Here `f` is a polynomial in F_2^m[X], `alpha` is a vector of elements in F_2^m,
 * `n` is the length of the RS code, we are using RM(1, m-1) duplicated `r` times
 */
void RS_RM_concatenated_decoding(GF2EX &res, const vec_GF2 &c,
                                      const vec_GF2E &alpha, const int n, const int k, const int r)
{

    int i, j;
    int m = GF2E::degree();
    int len = (1 << (m - 1));
    int rlen = r * len;

    /* Temporary vectors */
    vec_GF2E tmp_c1(INIT_SIZE,n);
    tmp_c1.SetLength(n);    
    vec_GF2 tmp_c2(INIT_SIZE, rlen);
    tmp_c2.SetLength(rlen);

    GF2X tmp_res(INIT_SIZE,m);
    tmp_res.SetLength(m);

    /* will decode each block of r * 2^len bits */
    for (i = 0; i < n; i++)
    {
        /* Decode a block of r*len bits using RM(1, m-1) code */
        int idx = i * rlen;

        tmp_c2 = VectorCopy(shift(c,-idx),rlen);
        
        clear(tmp_res); // Useless
        RM_decoding_duplicated(tmp_res, tmp_c2, m - 1, r);
        tmp_c1[i] = to_GF2E(tmp_res);
    }

    RS_decoder(res, tmp_c1, alpha, n, k);
}

/* Decoding as in HQC*/
void HQC_decoding(GF2EX &res, const GF2X &u, const GF2X &v, const GF2X &y,
                       const vec_GF2E &alpha, const int n, const int n1, const int k, const int r)
{

    int _m = GF2E::degree();
    int len = (1 << (_m - 1)); /* length of RM(1, m-1) code (n2)*/

    GF2X P = gf2x_set_cyclic(n);

    GF2X tmp_poly;
    tmp_poly = v + MulMod(u, y, P);

    tmp_poly = trunc(tmp_poly, r * len * n1);

    vec_GF2 tmp_vec(INIT_SIZE,r * len * n1);
    tmp_vec = conv<vec_GF2>(tmp_poly);
    RS_RM_concatenated_decoding(res, tmp_vec, alpha, n1, k, r);
}

