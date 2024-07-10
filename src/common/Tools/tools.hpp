#ifndef TOOLS_H
#define TOOLS_H


void fisher_yates(int* perm, const int n);

int int_check_repeat(const int* a, const int len);

int cm_random_indices(int* res, const int n, const int t, const int N, const int tau);

void cm_gen_e(int* e, const int n, const int t, const int N, const int tau);


#endif // TOOLS_H
