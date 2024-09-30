#include <stdlib.h>
#include <iostream>
// #include <cstdlib>
// #include <verilated.h>
// #include <verilated_vcd_c.h>

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

// #include "SerialPort.hpp"

#include "common/Tools/tools.hpp"
#include "common/Tools/codes.hpp"


#include "Tools/ClassicMcEliece/CM_secret_key.hpp"
#include "Tools/ClassicMcEliece/CM_public_key.hpp"
#include "Tools/ClassicMcEliece/CM_keygen.hpp"


#include "Tools/Bike/Bike_secret_key.hpp"
#include "Tools/Bike/Bike_public_key.hpp"
#include "Tools/Bike/Bike_keygen.hpp"



// #include "VTop_Level.h"
// #include "VerilatorSimulation.hpp"
// #include "MySource.hpp"
#include "Modules/Comparator/Comparator.hpp"
#include "Modules/SyndComparator/SyndComparator.hpp"
#include "Modules/CM_RandomFixedWeight/CM_RandomFixedWeight.hpp"
#include "Modules/Bike_RandomFixedWeight/Bike_RandomFixedWeight.hpp"
#include "Modules/CM_Encoder/CM_Encoder.hpp"
#include "Modules/CM_Decoder/CM_Decoder.hpp"
#include "Modules/Bike_Encoder/Bike_Encoder.hpp"
#include "Modules/Bike_Decoder/Bike_Decoder.hpp"


using namespace spu;
using namespace spu::module;
using namespace std;

// #define VERIF_START_TIME 7
// vluint64_t sim_time = 0;
// vluint64_t posedge_cnt = 0;

// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <vector>
// #include <random>

int main(int argc, char** argv, char** env) {

    fmpz_t p;
    slong m;
    fq_ctx_t ctx, ctx_q;
    
    fmpz_init_set_ui(p, 2);
    m = 7;
    fq_ctx_init_conway(ctx, p, m, "x");
    fq_ctx_init_conway(ctx_q, p, 1, "α");

    FLINT_TEST_INIT(state);
    

    /* ************************************************************************* */
    /*                                TESTS FOR BIKE                             */
    /* ************************************************************************* */

    int len = (1 << m) - 1;
    int k = len - 10;
    int t = (len-k)/2;
    
    fq_struct* alpha = _fq_vec_init(len, ctx);
    fq_struct* codeword = _fq_vec_init(len, ctx);
    fq_struct* c = _fq_vec_init(len, ctx);
    fq_struct* cs = _fq_vec_init(20, ctx);

    fq_vec_rand_distinct_2(alpha, len, ctx, state);


    


    /* fq_poly_t message;     fq_poly_t m1; */
    /* fq_poly_init(message, ctx);     fq_poly_init(m1, ctx); */

    /* /\* in flint : length of poly is degree+1 *\/ */
    /* fq_poly_randtest(message, state, k+1, ctx); */


    /* RS_encoding(codeword, message, alpha, len, ctx); */
    
    /* RS_decoder(m1, codeword, alpha, len, k, ctx); */

    
    /* printf("the decoding is correct: %d\n", fq_poly_equal(m1, message, ctx)); */

    
    // /* ************************************************************************* */
    // /*                                TESTS FOR BIKE                             */
    // /* ************************************************************************* */

    // /* !! needs r such that 2 is primitive mod r !!  */
    // /* int r = random_suitable_integer(6); */


    // int level = 5;

    
    // int r = 0;
    // int weight = 0;
    // int error_weight = 0;
    // int NbIter = 0;
    // int tau = 0;
    
    // Bike_params(r, weight, error_weight, level);
    // BGF_params(NbIter, tau, level);


    // cout << r << endl;
    // cout << NbIter << endl;
    
    
    // int len = 2*r;
    // int n = len;

    // const int FRAME_SIZE = len;
    // const int SYND_SIZE = r;
    // const int WEIGHT = weight;
    // const int ERROR_WEIGHT = error_weight;
    // const int NBITER = NbIter;
    // const int TAU = tau;

    
    // Bike_secret_key SK = Bike_secret_key(SYND_SIZE, &ctx_q);
    // Bike_public_key PK = Bike_public_key(SYND_SIZE, &ctx_q);
    // Bike_keygen_naive(SK, PK, weight);

    
    // module::Initializer   <int> initializer(FRAME_SIZE);
    // module::Incrementer   <int> incr1(FRAME_SIZE);
    // module::Finalizer     <int> finalizer(FRAME_SIZE);
    
    // module::Comparator comp(FRAME_SIZE);
    // module::Bike_RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT);
    // module::Bike_Encoder bike_encode(PK);
    // module::Bike_Decoder bike_decode(WEIGHT, NBITER, TAU, SK);

    
    // initializer   ["initialize::out" ] = randfixed   ["random_fixed_weight::input"];
    // randfixed   ["random_fixed_weight::output" ] = bike_encode   ["bike_encoder::input"];
    // bike_encode ["bike_encoder::output"] = bike_decode ["bike_decoder::input"];
    // bike_decode ["bike_decoder::output"] = comp ["compare::input1"];
    // randfixed   ["random_fixed_weight::output" ]= comp ["compare::input2"];
    // comp ["compare::output"] = finalizer ["finalize::in"  ];


    // std::vector<runtime::Task*> first = {&initializer("initialize")};

    // runtime::Sequence seq(first);

    // std::ofstream file("graph.dot");
    // seq.export_dot(file);

    // for (auto lt : seq.get_tasks_per_types())
    //     for (auto t : lt)
    // 	    {
    // 		t->set_stats(true);
    // 		t->set_debug(true);
    // 	    }

    // seq.exec_seq();
    
    
    

    
    /* ************************************************************************* */       
    /* ************************************************************************* */
    

    
    
    
    /* ************************************************************************* */
    /*                           TESTS FOR CLASSIC MCELIECE                      */
    /* ************************************************************************* */
    /* int deg = 15;  */
    /* slong len = 128; */
    /* int t = deg; */
    /* int tau = 20; */
    
    /* // AFTER THAT :     TESTS WITH MODULES */
    
    /* const int FRAME_SIZE = len; */
    /* const int DEG = deg; */
    /* const int WEIGHT = t; */
    /* const int TAU = tau; */
    /* const int N = 1 << m; */
    /* const int OUTPUT_SIZE = DEG * m; */

	
    /* CM_secret_key SK = CM_secret_key(FRAME_SIZE, &ctx); */
    /* CM_public_key PK = CM_public_key(FRAME_SIZE, m, DEG, &ctx_q); */

    
    /* CM_keygen_naive(SK, PK, FRAME_SIZE, DEG, ctx, state); */
    
    /* module::Initializer   <int> initializer(FRAME_SIZE); */
    /* module::Incrementer   <int> incr1(FRAME_SIZE); */
    /* module::Finalizer     <int> finalizer(FRAME_SIZE); */
    /* // module::Finalizer     <int> finalizer(FRAME_SIZE); */

    /* // // module::Finalizer     <int> finalizer_hw(FRAME_SIZE); */
    /* // module::MySource    my_source(FRAME_SIZE); */

    /* module::Comparator comp(FRAME_SIZE); */
    /* // module::CM_RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT, N, TAU); */
    /* module::Bike_RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT); */
    /* module::CM_Encoder cm_encode(FRAME_SIZE, OUTPUT_SIZE, PK); */
    /* module::CM_Decoder cm_decode(FRAME_SIZE, OUTPUT_SIZE, WEIGHT, SK); */


    /* // module::SyndComparator comp(FRAME_SIZE, OUTPUT_SIZE); */
    
    /* // module::Comparator comp_fpga(FRAME_SIZE); */
    /* // VerilatorSimulation sim(FRAME_SIZE); */
    /* // // SerialPort serial("/dev/tty.usbserial-210292ABF7641", 115200, FRAME_SIZE); */

    

    /* // initializer   ["initialize::out" ] = incr1           ["increment::in"]; */
    /* // // my_source   ["generate::output" ] = sim             ["simulate::input"]; */
    /* // // // my_source   ["generate::output" ] = serial          ["write::input"]; */

    /* // // sim         ["simulate::output" ] = comp_sim            ["compare::input2"]; */
    /* // // comp_sim    ["compare::output"  ] = finalizer_sw        ["finalize::in"]; */
    
    /* // comp ["compare::input1"] = incr1     ["increment::out"]; */
    /* // comp ["compare::input2"] = incr1     ["increment::out"]; */

    /* // randfixed["random_fixed_weight::input"] = comp["compare::output"]; */
    
    
    /* // randfixed["random_fixed_weight::output"] = finalizer ["finalize::in"  ]; */
    /* // comp ["compare::output"] = finalizer ["finalize::in"  ]; */

    /* // cout << "done \n" <<endl; */
    
    /* initializer   ["initialize::out" ] = randfixed   ["random_fixed_weight::input"]; */
    /* randfixed   ["random_fixed_weight::output" ] = cm_encode   ["cm_encoder::input"]; */
    /* cm_encode ["cm_encoder::output"] = cm_decode ["cm_decoder::input"]; */
    /* cm_decode ["cm_decoder::output"] = comp ["compare::input1"]; */
    /* randfixed   ["random_fixed_weight::output" ]= comp ["compare::input2"]; */
    /* comp ["compare::output"] = finalizer ["finalize::in"  ]; */


    /* // randfixed   ["random_fixed_weight::output" ] = comp ["compare::input1" ]; */
    /* // cm_decode ["cm_decoder::output"] = comp ["compare::input2" ]; */
    /* // cm_encode ["cm_encoder::output"] = comp ["compare::input3" ]; */

    /* // comp ["compare::output"] = finalizer ["finalize::in"]; */
    
    /* // randfixed["random_fixed_weight::output"] = finalizer ["finalize::in"  ]; */


    
    /* // // incr1       ["increment::out" ]   = comp_fpga            ["compare::input1"]; */
    /* // // serial      ["write::output"   ]   = comp_fpga            ["compare::input2"]; */
    /* // // comp_fpga   ["compare::output"  ] = finalizer_hw        ["finalize::in"]; */

    /* std::vector<runtime::Task*> first = {&initializer("initialize")}; */

    /* runtime::Sequence seq(first); */

    /* std::ofstream file("graph.dot"); */
    /* seq.export_dot(file); */

    /* for (auto lt : seq.get_tasks_per_types()) */
    /*     for (auto t : lt) */
    /* 	    { */
    /* 		t->set_stats(true); */
    /* 		t->set_debug(true); */
    /* 	    } */

    /* seq.exec_seq(); */
    /* // seq.exec_seq(); */



    

    
    /* ************************************************************************* */
    /* ************************************************************************* */
    
}
