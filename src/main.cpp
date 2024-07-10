#include <stdlib.h>
#include <iostream>
// #include <cstdlib>
// #include <verilated.h>
// #include <verilated_vcd_c.h>
// #include "VTop_Level.h"

// #include "VerilatorSimulation.hpp"
// #include "MySource.hpp"
#include "Comparator/Comparator.hpp"
#include "RandomFixedWeight/RandomFixedWeight.hpp"

// #include "SerialPort.hpp"

#include <vector>
#include <streampu.hpp>


#include <flint.h>
#include <fmpz.h>		/* large integers */

using namespace spu;
using namespace spu::module;

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
    fmpz_init_set_ui(p, 2);
    flint_printf("test");
  
    const int FRAME_SIZE = 20;
    const int WEIGHT = 10;
    const int TAU = 20;
    const int N = 25;

    module::Initializer   <int> initializer(FRAME_SIZE);
    module::Incrementer   <int> incr1(FRAME_SIZE);
    module::Finalizer     <int> finalizer(FRAME_SIZE);

    // // module::Finalizer     <int> finalizer_hw(FRAME_SIZE);
    // module::MySource    my_source(FRAME_SIZE);

    module::Comparator comp(FRAME_SIZE);
    module::RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT, N, TAU);
    
    // module::Comparator comp_fpga(FRAME_SIZE);
    // VerilatorSimulation sim(FRAME_SIZE);
    // // SerialPort serial("/dev/tty.usbserial-210292ABF7641", 115200, FRAME_SIZE);

    

    initializer   ["initialize::out" ] = incr1           ["increment::in"];
    // my_source   ["generate::output" ] = sim             ["simulate::input"];
    // // my_source   ["generate::output" ] = serial          ["write::input"];

    // sim         ["simulate::output" ] = comp_sim            ["compare::input2"];
    // comp_sim    ["compare::output"  ] = finalizer_sw        ["finalize::in"];
    
    comp ["compare::input1"] = incr1     ["increment::out"];
    comp ["compare::input2"] = incr1     ["increment::out"];

    randfixed["random_fixed_weight::input"] = comp["compare::output"];
    
    
    randfixed["random_fixed_weight::output"] = finalizer ["finalize::in"  ];

    // // incr1       ["increment::out" ]   = comp_fpga            ["compare::input1"];
    // // serial      ["write::output"   ]   = comp_fpga            ["compare::input2"];
    // // comp_fpga   ["compare::output"  ] = finalizer_hw        ["finalize::in"];

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
    // seq.exec_seq();

}
