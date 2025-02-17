#include <stdlib.h>
#include <iostream>


#include <vector>
#include <streampu.hpp>

#include <flint/flint.h>
#include <flint/fmpz.h>		/* large integers */
#include <flint/fq.h>		/* finite fields */

#include <flint/fq_poly.h>	/* pol. in finite fields */
#include <flint/fmpz_poly.h>	/* pol. in integers */
#include <flint/fmpz_vec.h>	/* vectors integers */
#include <flint/fq_vec.h>	/* vectors finite fields */
#include <flint/perm.h>		/* permutations */
#include <flint/fq_mat.h>	/* matrix / finite fields */

#include <pqc.hpp>

using namespace spu;
using namespace spu::module;
using namespace std;


int main(int argc, char** argv, char** env) {

    fmpz_t p;
    slong m;
    fq_ctx_t ctx, ctx_q;
    
    fmpz_init_set_ui(p, 2);
    m = 8;
    fq_ctx_init_conway(ctx, p, m, "x");
    fq_ctx_init_conway(ctx_q, p, 1, "α");

    FLINT_TEST_INIT(state);
    



    /* ************************************************************************* */
    /*                                TESTS FOR BIKE                             */
    /* ************************************************************************* */

    // !! needs r such that 2 is primitive mod r !!
#if 1
    int r = random_suitable_integer(6);
#else
    int r = 0;
#endif

    int level = 1;

    std::cout << "BIKE Level: " << level << std::endl;
    
    int weight = 0;
    int error_weight = 0;
    int NbIter = 0;
    int tau = 0;
    
    Bike_params(r, weight, error_weight, level);
    BGF_params(NbIter, tau, level);
    
    
    int len = 2*r;
    int n = len;

    const int FRAME_SIZE = len;
    const int SYND_SIZE = r;
    const int WEIGHT = weight;
    const int ERROR_WEIGHT = error_weight;
    const int NBITER = NbIter;
    const int TAU = tau;

    
    Bike_secret_key SK = Bike_secret_key(SYND_SIZE, &ctx_q);
    Bike_public_key PK = Bike_public_key(SYND_SIZE, &ctx_q);
    Bike_keygen_naive(SK, PK, weight);

    
    module::Initializer   <int> initializer(FRAME_SIZE);
    module::Incrementer   <int> incr1(FRAME_SIZE);
    module::Finalizer     <int> finalizer(FRAME_SIZE);
    
    module::Comparator comp(FRAME_SIZE);
    module::Bike_RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT);
    module::Bike_Encoder bike_encode(PK);
    module::Bike_Decoder bike_decode(WEIGHT, NBITER, TAU, SK);

    
    initializer   ["initialize::out" ] = randfixed   ["random_fixed_weight::input"];
    randfixed   ["random_fixed_weight::output" ] = bike_encode   ["bike_encoder::input"];
    bike_encode ["bike_encoder::output"] = bike_decode ["bike_decoder::input"];
    bike_decode ["bike_decoder::output"] = comp ["compare::input1"];
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
		t->set_debug(false);
	    }

    for(int i = 0; i<10; i++){
        auto t_start = std::chrono::steady_clock::now();

        seq.exec_seq();

        std::chrono::nanoseconds duration = std::chrono::steady_clock::now() - t_start;

        auto elapsed_time = duration.count() / 1000.f / 1000.f;
        std::cout << "Sequence elapsed time: " << elapsed_time << " ms" << std::endl;
    }
    

    tools::Stats::show(seq.get_modules_per_types(), true, false);
    
    
    

    
    /* ************************************************************************* */       
    /* ************************************************************************* */
    

    
    
    
    
}
