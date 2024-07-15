#include <stdlib.h>
#include <iostream>
// #include <cstdlib>
// #include <verilated.h>
// #include <verilated_vcd_c.h>

#include <vector>
#include <streampu.hpp>


#include <flint.h>
#include <fmpz.h>		/* large integers */
#include <fq.h>		/* finite fields */

// #include flint/fq_poly.h"	/* pol. in finite fields */
// #include "flint/fmpz_poly.h"	/* pol. in integers */
// #include "flint/fmpz_vec.h"	/* vectors integers */
// #include "flint/fq_vec.h"	/* vectors finite fields */
// #include "flint/perm.h"		/* permutations */
// #include "flint/fq_mat.h"	/* matrix / finite fields */



// #include "VTop_Level.h"
// #include "VerilatorSimulation.hpp"
// #include "MySource.hpp"
#include "Modules/Comparator/Comparator.hpp"
#include "Modules/RandomFixedWeight/RandomFixedWeight.hpp"

// #include "SerialPort.hpp"

#include "Tools/tools.hpp"
#include "Tools/codes.hpp"
#include "Tools/CM_secret_key.hpp"
#include "Tools/CM_public_key.hpp"

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
    slong d;
    fq_ctx_t ctx, ctx_q;
    
    fmpz_init_set_ui(p, 2);
    d = 7;
    fq_ctx_init_conway(ctx, p, d, "X");
    fq_ctx_init_conway(ctx_q, p, 1, "Y");

    int deg = 15; 
    slong len = 128;
    int t = deg/2;
    int tau = 20;



    

    // int a = _fq_vec_print(SK.get_alpha(), len, ctx); 

    // cout << " "  << endl;

    // fq_poly_t& h  = SK.g;
    // // h = SK.get_g();
    // fq_poly_print_pretty(h, "X",  ctx);
    // cout << " "  << endl;

    // fq_poly_one(h, ctx);
    // fq_poly_print_pretty(h, "X",  ctx);
    // cout << " "  << endl;
    
    // fq_poly_print_pretty(SK.g, "X",  ctx);
    // cout << " "  << endl;


    FLINT_TEST_INIT(state);

    CM_secret_key SK = CM_secret_key(len, &ctx);
    CM_public_key PK = CM_public_key(len, d, t, &ctx_q);
    
    int res = 0;

    
    while (!res) {
	SK.keygen(t, state);
   
	cout << "after keygen" << endl;

	// fq_struct* beta = _fq_vec_init(len, ctx);
	// _fq_vec_set(beta, SK.get_alpha(), len, ctx);
	// a = _fq_vec_print(SK.get_alpha(), len, ctx);
	// for (int i = 0; i < len; i++) {
	//     fq_print_pretty(&beta[i], ctx);
	//     cout  <<  endl;
	// }
    
	// cout << " "  << endl;
    
	// fq_poly_print_pretty(SK.g, "Y",  ctx);
	// cout << " "  << endl;

    
	
	// fq_mat_t& M = PK.T;

	res = PK.keygen(SK, ctx);
	
	cout << "keygen success ? \n" << res << endl;
    }
    

    // AFTER THAT :     TESTS WITH MODULES
    
    // const int FRAME_SIZE = 20;
    // const int WEIGHT = 10;
    // const int TAU = 20;
    // const int N = 25;

    // module::Initializer   <int> initializer(FRAME_SIZE);
    // module::Incrementer   <int> incr1(FRAME_SIZE);
    // module::Finalizer     <int> finalizer(FRAME_SIZE);

    // // // module::Finalizer     <int> finalizer_hw(FRAME_SIZE);
    // // module::MySource    my_source(FRAME_SIZE);

    // module::Comparator comp(FRAME_SIZE);
    // module::RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT, N, TAU);
    
    // // module::Comparator comp_fpga(FRAME_SIZE);
    // // VerilatorSimulation sim(FRAME_SIZE);
    // // // SerialPort serial("/dev/tty.usbserial-210292ABF7641", 115200, FRAME_SIZE);

    

    // initializer   ["initialize::out" ] = incr1           ["increment::in"];
    // // my_source   ["generate::output" ] = sim             ["simulate::input"];
    // // // my_source   ["generate::output" ] = serial          ["write::input"];

    // // sim         ["simulate::output" ] = comp_sim            ["compare::input2"];
    // // comp_sim    ["compare::output"  ] = finalizer_sw        ["finalize::in"];
    
    // comp ["compare::input1"] = incr1     ["increment::out"];
    // comp ["compare::input2"] = incr1     ["increment::out"];

    // randfixed["random_fixed_weight::input"] = comp["compare::output"];
    
    
    // randfixed["random_fixed_weight::output"] = finalizer ["finalize::in"  ];

    // // // incr1       ["increment::out" ]   = comp_fpga            ["compare::input1"];
    // // // serial      ["write::output"   ]   = comp_fpga            ["compare::input2"];
    // // // comp_fpga   ["compare::output"  ] = finalizer_hw        ["finalize::in"];

    // std::vector<runtime::Task*> first = {&initializer("initialize")};

    // runtime::Sequence seq(first);

    // std::ofstream file("graph.dot");
    // seq.export_dot(file);

    // for (auto lt : seq.get_tasks_per_types())
    //     for (auto t : lt)
    //     {
    //         t->set_stats(true);
    //         t->set_debug(true);
    //     }

    // seq.exec_seq();
    // // seq.exec_seq();

}
