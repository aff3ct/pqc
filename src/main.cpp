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
#include "Modules/Bike_Encoder/Bike_Encoder.hpp"
#include "Modules/CM_Decoder/CM_Decoder.hpp"



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

    /* !! needs r such that 2 is primitive mod r !!  */

    int r = random_suitable_integer(9);
    
    
    // int r = 101;
    int len = 2*r;
    int n = len;
    int weight = 10;
    
    
    Bike_secret_key SK = Bike_secret_key(r, &ctx_q);
    Bike_public_key PK = Bike_public_key(r, &ctx_q);
        
    Bike_keygen_naive(SK, PK, weight);
    
    fq_mat_t H0, H1;
    fq_mat_init(H0, r, r, ctx_q);
    fq_mat_init(H1, r, r, ctx_q);
    fq_poly_t P;
    fq_poly_init(P, ctx_q);
    fq_poly_set_cyclic(P, r, ctx_q);


    fq_mult_matrix(H0, SK.h0, P, ctx_q);
    fq_mult_matrix(H1, SK.h1, P, ctx_q);

    fq_mat_print_pretty(H0, ctx_q);

    fq_mat_t H;
    fq_mat_init(H, r, n, ctx_q);
    fq_mat_concat_horizontal(H, H0, H1, ctx_q);
    fq_mat_print_pretty(H, ctx_q);

    int e[n];
    for (int i = 0; i < n; ++i) {
	e[i] = 0;
    }
    bike_gen_e(e, n, weight);

    fq_struct* ee = _fq_vec_init(n, ctx_q);
    _int_vec_2_fq(ee, e, n, ctx_q);
    

    fq_struct* s = _fq_vec_init(r, ctx_q);
    fq_poly_t sp; fq_poly_init(sp, ctx_q);
    fq_poly_zero(sp, ctx_q);
    
    Bike_encoding(sp, e, PK.h,  r, ctx_q);

    fq_struct* ss = _fq_vec_init(r, ctx_q);

    
    // fq_poly_set_coeffs(sp, s, r, ctx_q);

    fq_poly_mul(sp, sp, SK.h0, ctx_q);
    for (int k = 0; k < r; k++) {
	fq_poly_get_coeff(&ss[k], sp, k, ctx_q);
    }
    
    
    fq_struct* res = _fq_vec_init(n, ctx_q);
    
    // fq_mat_mul_vec(s, H, ee, n, ctx_q);

    int b = Bike_decoding(res, ss, H, weight, 5, 3, ctx_q);

    printf("%d \n", b);

    
    _fq_vec_print(res, n, ctx_q);
    printf("\n");

    for (int ii = 0; ii < n; ii++) {
	fq_print_pretty(&ee[ii], ctx_q);
	    printf(" ");
    }
    printf("\n");
    // _fq_vec_print(ee, n, ctx_q);
    
    
    
    // const int FRAME_SIZE = len;
    // const int WEIGHT = weight;

    // module::Initializer   <int> initializer(FRAME_SIZE);
    // module::Incrementer   <int> incr1(FRAME_SIZE);
    // module::Finalizer     <int> finalizer(r);

    // module::Bike_RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT);
    // module::Bike_Encoder bike_encode(PK);


        
    // initializer   ["initialize::out" ] = randfixed   ["random_fixed_weight::input"];
    // randfixed   ["random_fixed_weight::output" ] = bike_encode   ["bike_encoder::input"];
    // bike_encode ["bike_encoder::output"] = finalizer ["finalize::in"  ];


    // // randfixed   ["random_fixed_weight::output" ] = comp ["compare::input1" ];
    // // cm_decode ["cm_decoder::output"] = comp ["compare::input2" ];
    // // cm_encode ["cm_encoder::output"] = comp ["compare::input3" ];

    // // comp ["compare::output"] = finalizer ["finalize::in"];
    
    // // randfixed["random_fixed_weight::output"] = finalizer ["finalize::in"  ];


    
    // // // incr1       ["increment::out" ]   = comp_fpga            ["compare::input1"];
    // // // serial      ["write::output"   ]   = comp_fpga            ["compare::input2"];
    // // // comp_fpga   ["compare::output"  ] = finalizer_hw        ["finalize::in"];

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
    // // seq.exec_seq();


    
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



    
    /* fq_ctx_clear(ctx); */
    /* fq_ctx_clear(ctx_q); */




    /* int a = _fq_vec_print(SK.get_alpha(), len, ctx);  */

    /* cout << " "  << endl; */

    /* fq_poly_t& h  = SK.g; */
    /* // h = SK.get_g(); */
    /* fq_poly_print_pretty(h, "X",  ctx); */
    /* cout << " "  << endl; */

    /* fq_poly_one(h, ctx); */
    /* fq_poly_print_pretty(h, "X",  ctx); */
    /* cout << " "  << endl; */
    
    /* fq_poly_print_pretty(SK.g, "X",  ctx); */
    /* cout << " "  << endl;  */

    
    /* ************************************************************************* */
    /* ************************************************************************* */
    
}
