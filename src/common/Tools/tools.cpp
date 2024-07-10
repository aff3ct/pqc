#include <iostream>
#include <random>

#include "tools.hpp"

using namespace std;



void
fisher_yates(int* perm, const int n) {
  int j;
  
  std::random_device rd;
  std::mt19937 rand_gen(rd());
  
  
  for (int i = 0; i < n; ++i) {
      std::uniform_int_distribution<std::mt19937::result_type> dis(0, i);
      j = dis(rand_gen);
      if (i != j) {
	  perm[i] = perm[j];
      }
      perm[j] = i;
  }
}


/* check repetitions in a naive way */
int
int_check_repeat(const int* a, const int len) {
  for (int i = 0; i < len; i++){
    for (int j = i+1; j < len; j++) {
      if (a[i] == a[j]) return 1;
    }
  }
  return 0;
}



/*  computes t random and pairwise distinct indexes in [1,n] as in Classic
    McEliece
    draw random τ > t elts in [1, N[ and take the t first ones that are in [1,n]
    return the number of elts in [1,n] effectively computed
    Typically in CM, N = 2^m
*/
int
cm_random_indices(int* res, const int n, const int t, const int N, const int tau) {
    int ind = 0;
    int count_loop = 0;
    int a;

    // rand with mersenne
    std::random_device rd;
    std::mt19937 rand_gen(rd());
    std::uniform_int_distribution<std::mt19937::result_type> dis(0, N-1);
    
    while (ind < t && count_loop < tau) {
	a = dis(rand_gen);
	if (a < n) {
	    res[ind] = a;
	    ind++;
	}
	count_loop++;
    }
    return ind;
}


/* computes e as in Classic McEliece, more or less */
void
cm_gen_e(int* e, const int n, const int t, const int N, const int tau) {
    int len, repet = 1;
    int inds[t];
    while (repet) {
	len = cm_random_indices(inds, n, t, N, tau); /* compute random set of indices */
	if (len == t) {				 /* check if we obtained enough indices */
	    repet = int_check_repeat(inds, t);
	}
    }
    for (int i = 0; i <t; i++) {
	e[inds[i]] = 1;
    }
}


// int main(int argc, char *argv[])
// {
//   int n = 20;
//   int t = 5;
//   int tau = 10;
//   int N = 50;	  // essentially 2^m in CM
  
//   int perm[n];


//   int e[n];
//   for (int i = 0; i < n; i++) e[i] = 0;
  
//   fisher_yates(perm, n);

//   // for (int i = 0 ; i < n; ++i) {
//   //     cout << perm[i] << endl;
//   // }
  
//   cm_gen_e(e, n, t,  N,  tau);

//   for (int i = 0 ; i < n; ++i) {
//       cout << e[i] << endl;
//   }
  
  
//   return 0;
// }
