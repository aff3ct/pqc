#ifndef CODES_H
#define CODES_H

#include "tools.hpp"

// ENCODER
void RM_encoding_duplicated(vec_GF2 &res, const mat_GF2 &G, const vec_GF2 &alpha, int m, int r);

void RS_RM_concatenated_encoding(vec_GF2 &out, const GF2EX &f,
								 const vec_GF2E &alpha, const int n, const int r, int m);

void CM_encoding(vec_GF2 &res, const mat_GF2 &T, const vec_GF2 &e);

void Bike_encoding(GF2X &res, const int *e, const GF2X &h, const int r);

void HQC_encoding(GF2X &res1, GF2X &res2, const GF2EX &m, const GF2X &h,
					   const GF2X &s, const vec_GF2E &alpha, const int n, const int k,
					   const int n1, const int r, const int we, const int wr, int _m, GF2X &P);

// DECODER
void RM_decoding_duplicated(GF2X &res, const vec_GF2 &c, const int m, const int r);

void RS_RM_concatenated_decoding(GF2EX &res, const vec_GF2 &c,
									  const vec_GF2E &alpha, const int n, const int k, const int r);

int RS_decoder(GF2EX &res, const Vec<GF2E> &c, const Vec<GF2E> &alpha, int n, int k);

int GRS_decoder(Vec<GF2E> &res, const Vec<GF2E> &c, const Vec<GF2E> &alpha, const Vec<GF2E> &beta, int n, int k);

int Goppa_decoder_bin(Vec<GF2E> &res, const Vec<GF2E> &c, const Vec<GF2E> &alpha,
						   const GF2EX &g, const int n, const int k);

mat_GF2E Goppa_parity_check(const vec_GF2E &alpha, const GF2EX &g, int row, int col);

void HQC_decoding(GF2EX &res, const GF2X &u, const GF2X &v, const GF2X &y,
					   const vec_GF2E &alpha, const int n, const int n1, const int k, const int r);

int Bike_decoding_v2(vec_GF2 &res, const GF2X &s, const GF2X &h0,
						  const GF2X &h1, const int r, const int weight,
						  const int NbIter, const int tau);

int CM_syndrome_decoding_bin(Vec<GF2E> &res, const Vec<GF2E> &s, const Vec<GF2E> &alpha, const GF2EX &g,
								  const int len, const int t);

#endif // CODES_H
