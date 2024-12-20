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


#include "Tools/ClassicMcEliece/CM_secret_key.hpp"
#include "Tools/ClassicMcEliece/CM_public_key.hpp"
#include "Tools/ClassicMcEliece/CM_keygen.hpp"



#include "Modules/RandomVector/RandomVector.hpp"
#include "Modules/Comparator/Comparator.hpp"
#include "Modules/SyndComparator/SyndComparator.hpp"
#include "Modules/CM_RandomFixedWeight/CM_RandomFixedWeight.hpp"
#include "Modules/CM_Encoder/CM_Encoder.hpp"
#include "Modules/CM_Decoder/CM_Decoder.hpp"


using namespace spu;
using namespace spu::module;
using namespace std;


int main(int argc, char** argv, char** env) {

    FLINT_TEST_INIT(state);

    fmpz_t p;
    fq_ctx_t ctx, ctx_q;
    
    fmpz_init_set_ui(p, 2);


    
    
    
    /* ************************************************************************* */
    /*                           TESTS FOR CLASSIC MCELIECE                      */
    /* ************************************************************************* */
    int m, n, t, tau;

    int level = 1;

    CM_params(m, n, t, level);
    if ((1 << m) == n) {
	tau = t;
    } else {
	tau = 2*t;
    }
    

    
    const int FRAME_SIZE = n;
    const int DEG = t;
    const int WEIGHT = t;
    const int N = 1 << m;
    const int TAU = tau;
    const int OUTPUT_SIZE = DEG * m;


    fq_ctx_init_conway(ctx, p, m, "x");
    fq_ctx_init_conway(ctx_q, p, 1, "α");
	
    CM_secret_key SK = CM_secret_key(FRAME_SIZE, &ctx);
    CM_public_key PK = CM_public_key(FRAME_SIZE, m, DEG, &ctx_q);

    
    CM_keygen_naive(SK, PK, FRAME_SIZE, DEG, ctx, state);
    
    module::Initializer   <int> initializer(FRAME_SIZE);
    module::Incrementer   <int> incr1(FRAME_SIZE);
    module::Finalizer     <int> finalizer(FRAME_SIZE);

    module::Comparator comp(FRAME_SIZE);
    module::CM_RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT, N, TAU);
    module::CM_Encoder cm_encode(FRAME_SIZE, OUTPUT_SIZE, PK);
    module::CM_Decoder cm_decode(FRAME_SIZE, OUTPUT_SIZE, WEIGHT, SK);

    
    initializer   ["initialize::out" ] = randfixed   ["random_fixed_weight::input"];
    randfixed   ["random_fixed_weight::output" ] = cm_encode   ["cm_encoder::input"];
    cm_encode ["cm_encoder::output"] = cm_decode ["cm_decoder::input"];
    cm_decode ["cm_decoder::output"] = comp ["compare::input1"];
    randfixed   ["random_fixed_weight::output" ]= comp ["compare::input2"];
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
