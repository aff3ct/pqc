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
    m = 8;
    fq_ctx_init_conway(ctx, p, m, "x");
    fq_ctx_init_conway(ctx_q, p, 1, "α");

    FLINT_TEST_INIT(state);
    


    /* ************************************************************************* */
    /*                                 TESTS FOR HQC                             */
    /* ************************************************************************* */

    /* int n1 = 46;		/\* length for RS code *\/ */
    /* int k = 16;			/\* degree for RS code *\/ */
    /* int r = 3;			/\* number of times the RM code is duplicated *\/ */
    /* int n2 = 1 << (m-1);	/\* length of RM(1, m) code *\/ */
    /* int len = r * n1 * n2; */
    /* int n = 17669; */
    
    /* cout << "length is: " << len << endl; */
    
    /* int weight = 66; */
    /* int we = 75; */
    /* int wr = 75; */
        
    /* fq_struct* alpha = _fq_vec_init(n1, ctx); */
    /* fq_struct* codeword = _fq_vec_init(len, ctx_q); */
    /* fq_struct* codeword1 = _fq_vec_init(len, ctx_q); */
    /* fq_struct* codeword2 = _fq_vec_init(len, ctx_q); */
    
    /* HQC_secret_key SK = HQC_secret_key(n, &ctx_q); */
    /* HQC_public_key PK = HQC_public_key(n, n1, &ctx, &ctx_q); */
    /* HQC_keygen_naive(SK, PK, weight, state); */


    /* _fq_vec_zero(codeword, len, ctx_q); */
    /* _fq_vec_zero(codeword1, len, ctx_q); */
    /* _fq_vec_zero(codeword2, len, ctx_q); */
    
    /* fq_vec_rand_distinct_2(alpha, n1, ctx, state); */

    /* fq_poly_t message;     fq_poly_t m1; */
    /* fq_poly_init(message, ctx_q);     fq_poly_init(m1, ctx_q); */

    
    /* fq_poly_t res; fq_poly_init(res, ctx_q); */
    
    /* /\* in flint : length of poly is degree+1 *\/ */
    /* fq_poly_randtest(message, state, k, ctx_q); */

    /* RS_RM_concatenated_encoding(codeword, message, PK.alpha, n1, ctx, r, ctx_q); */
    
    /* int ee[len]; */
    /* hqc_gen_e(ee, len, 7000); */
    /* fq_t tmp; fq_init(tmp, ctx_q); */
    /* for (int i = 0; i < len; i++) { */
    /* 	fq_set_ui(tmp, ee[i], ctx_q); */
    /* 	fq_add(&codeword1[i], &codeword[i], tmp, ctx_q); */
    /* } */
    
    /* _fq_vec_add(codeword2, codeword, codeword1, len, ctx_q); */
    /* // _fq_vec_print_pretty(codeword2, len, ctx_q);     */
    
    /* cout << "hamming distance: " << hamming_distance(codeword, codeword1, len, ctx_q) << endl; */
    /* cout << "hamming weight: " << hamming_weight(codeword2, len, ctx_q) << endl; */
    
    /* RS_RM_concatenated_decoding(res, codeword1, PK.alpha, n1, k, ctx, r, ctx_q); */

    /* cout << "decoding successful ?  " << fq_poly_equal(message, res, ctx) << endl; */

    /* fq_poly_t res1; fq_poly_init(res1, ctx_q); fq_poly_zero(res1, ctx_q); */
    /* fq_poly_t res2; fq_poly_init(res2, ctx_q); fq_poly_zero(res2, ctx_q); */
    /* fq_poly_t res3; fq_poly_init(res3, ctx_q); fq_poly_zero(res3, ctx_q); */
    /* fq_struct* res4 = _fq_vec_init(len, ctx_q);  */
    
    /* HQC_encoding(res1, res2, message, PK.h, PK.s, PK.alpha, n, k, n1, r,  we, wr, ctx, */
    /* 		 ctx_q); */
    
    /* HQC_decoding(res3, res1, res2, SK.y, PK.alpha, n, n1, k, r, ctx, ctx_q); */

    /* cout << "decoding working ?    " << fq_poly_equal(message, res3, ctx_q) << endl; */


    /* fq_poly_print_pretty(message, "T", ctx_q); cout << endl; */
    /* fq_poly_print_pretty(res3, "T", ctx_q); */
    
    int n1 = 46;		/* length for RS code */
    int k = 16;			/* degree for RS code */
    int r = 3;			/* number of times the RM code is duplicated */
    int n2 = 1 << (m-1);	/* length of RM(1, m) code */
    int len = r * n1 * n2;
    int n = 17669;
    
    cout << "length is: " << len << endl;
    
    int weight = 66;
    int we = 75;
    int wr = 75;
    
    const int FRAME_SIZE = k;
    const int GLOBAL_LENGTH = n;
    const int WEIGHT = weight;
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

    
    

    // /* ************************************************************************* */
    // /*                                TESTS FOR BIKE                             */
    // /* ************************************************************************* */

    /* !! needs r such that 2 is primitive mod r !!  */
    /* int r = random_suitable_integer(6); */


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
    /* Bike_public_key PK = Bike_public_key(SYND_SIZE, &ctx_q); */
    /* Bike_keygen_naive(SK, PK, weight); */

    
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
