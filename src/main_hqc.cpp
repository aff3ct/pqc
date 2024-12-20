#include <stdlib.h>
#include <iostream>


#include <vector>
#include <streampu.hpp>

#include <flint/flint.h>
#include <flint/fmpz.h>		/* large integers */
#include <flint/fq.h>		/* finite fields */

#include "flint/fq_poly.h"	/* pol. in finite fields */
#include "flint/fmpz_poly.h"	/* pol. in integers */
#include "flint/fmpz_vec.h"	/* vectors integers */
#include "flint/fq_vec.h"	/* vectors finite fields */
#include "flint/perm.h"		/* permutations */
#include "flint/fq_mat.h"	/* matrix / finite fields */


#include "common/Tools/tools.hpp"
#include "common/Tools/codes.hpp"


#include "Tools/HQC/HQC_secret_key.hpp"
#include "Tools/HQC/HQC_public_key.hpp"
#include "Tools/HQC/HQC_keygen.hpp"


#include "Modules/RandomVector/RandomVector.hpp"
#include "Modules/Comparator/Comparator.hpp"
#include "Modules/SyndComparator/SyndComparator.hpp"
#include "Modules/HQC_Encoder/HQC_Encoder.hpp"
#include "Modules/HQC_Decoder/HQC_Decoder.hpp"


using namespace spu;
using namespace spu::module;
using namespace std;


int main(int argc, char** argv, char** env) {
    /* ************************************************************************* */
    /*                                 TESTS FOR HQC                             */
    /* ************************************************************************* */


    FLINT_TEST_INIT(state);	/* for randomness */

    
    fmpz_t p;
    slong m;
    fq_ctx_t ctx, ctx_q;
    
    fmpz_init_set_ui(p, 2);
    m = 8;
    fq_ctx_init_conway(ctx, p, m, "x");
    fq_ctx_init_conway(ctx_q, p, 1, "α");

    int level = 1;

    /* int k = 16;			/\* degree for RS code *\/ */
    /* int n1 = 46;		/\* length for RS code *\/ */
    /* int n2 = 1 << (m-1);	/\* length of RM(1, m) code *\/ */
    /* int r = 3;			/\* number of times the RM code is duplicated *\/ /
    /* int n = 17669;  */
    /*     int weight = 66; */
    /* int we = 75; */
    /* int wr = 75;     */


    int k, n1, n2, r, n, w, we, wr;
    HQC_params(k, n1, n2, r, n, w, we, wr, level);
    
    int len = r * n1 * n2;


    
    
    const int FRAME_SIZE = k;
    const int GLOBAL_LENGTH = n;
    const int WEIGHT = w;
    const int ERROR_WEIGHT = we;
    const int N1 = n1;
    const int CODE_LENGTH = len;
    const int R = r;
    
    HQC_secret_key SK = HQC_secret_key(GLOBAL_LENGTH, &ctx_q);
    HQC_public_key PK = HQC_public_key(GLOBAL_LENGTH, N1, &ctx, &ctx_q);
    HQC_keygen_naive(SK, PK, WEIGHT, state);
    
    module::Initializer   <int> initializer(FRAME_SIZE);
    module::Finalizer     <int> finalizer(FRAME_SIZE);
    
    module::Comparator comp(FRAME_SIZE);
    module::RandomVector random_vector(FRAME_SIZE);
    module::HQC_Encoder hqc_encode(PK, FRAME_SIZE, CODE_LENGTH, R, ERROR_WEIGHT);
    module::HQC_Decoder hqc_decode(SK, PK, FRAME_SIZE, CODE_LENGTH, R);
    

    
    initializer   ["initialize::out" ] = random_vector   ["random_vector::input"];
    random_vector   ["random_vector::output" ] = hqc_encode   ["hqc_encoder::input"];
    hqc_encode ["hqc_encoder::output1"] = hqc_decode ["hqc_decoder::input1"];
    hqc_encode ["hqc_encoder::output2"] = hqc_decode ["hqc_decoder::input2"];

    hqc_decode ["hqc_decoder::output"] = comp ["compare::input1"];
    random_vector   ["random_vector::output" ] = comp ["compare::input2"];
    comp ["compare::output"] = finalizer ["finalize::in"  ];


    std::vector<runtime::Task*> first = {&initializer("initialize")};

    runtime::Sequence seq(first);

    std::ofstream file("graph.dot");
    seq.export_dot(file);

    for (auto lt : seq.get_tasks_per_types())
        for (auto t : lt)
	    {
		t->set_stats(true);
		t->set_debug(true);
	    }

    seq.exec_seq();

    
    




    

    
    /* ************************************************************************* */
    /* ************************************************************************* */
   
}
