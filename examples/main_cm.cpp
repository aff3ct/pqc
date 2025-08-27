#include <stdlib.h>
#include <iostream>
#include <getopt.h>

#include <vector>
#include <streampu.hpp>

#include <flint/flint.h>

#include <pqc.hpp>

using namespace spu;
using namespace spu::module;
using namespace std;

int main(int argc, char **argv, char **env)
{
    option longopts[] =
    { { "level", required_argument, NULL, 'l' },
      { "throughput", no_argument, NULL, 'b' },
        { "help", no_argument, NULL, 'h' },
        { NULL, 0, NULL, 0 } };


    FLINT_TEST_INIT(state);

    /* ************************************************************************* */
    /*                           TESTS FOR CLASSIC MCELIECE                      */
    /* ************************************************************************* */
    int m, n, t, tau;
    GF2X P; // Context mod polynomial

    bool throughput = false;
    int level = 1;

    while (1)
    {
        const int opt = getopt_long(argc, argv, "bl:h", longopts, 0);
        if (opt == -1) break;
        switch (opt)
        {
            case 'b':
                throughput = true;
                break;
            case 'l':
                level = atoi(optarg);
                switch (level)
                {
                        case 1:
                        case 3:
                        case 5:
                            break;
                        default:
                            std::cerr << "Invalid safety level argument: " << level << ", expected 1, 3 or 5" << std::endl;
                            exit(1);
                }
                break;
            case 'h':
                std::cerr << "usage: " << argv[0] << " [options]" << std::endl;
                std::cerr << std::endl;
                std::cerr << "  -l, --level            "
                          << "Safety level (1, 3, 5)   "
                          << "[" << level << "]" << std::endl;
                std::cerr << "  -h, --help             "
                          << "This help                "
                          << std::endl;
                exit(0);
                break;
            default:
                break;
        }
    }

    std::cout << "CM Level: " << level << std::endl;

    CM_params(m, n, t, P, level);
    if ((1 << m) == n)
    {
        tau = t;
    }
    else
    {
        tau = 2 * t;
    }

    const int FRAME_SIZE = n;
    const int DEG = t;
    const int WEIGHT = t;
    const int N = 1 << m;
    const int TAU = tau;
    const int OUTPUT_SIZE = DEG * m;

    //  Init GF2E
    GF2E::init(P);

    CM_secret_key SK = CM_secret_key(FRAME_SIZE);
    CM_public_key PK = CM_public_key(FRAME_SIZE, m, DEG);
    CM_keygen_naive(SK, PK, FRAME_SIZE, DEG, state);
    
    module::Initializer<int> initializer(FRAME_SIZE);
    module::Incrementer<int> incr1(FRAME_SIZE);
    module::Finalizer<int> finalizer(FRAME_SIZE);

    module::Comparator comp(FRAME_SIZE);
    module::CM_RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT, N, TAU);
    module::CM_Encoder cm_encode(FRAME_SIZE, OUTPUT_SIZE, PK);
    module::CM_Decoder cm_decode(FRAME_SIZE, OUTPUT_SIZE, WEIGHT, SK);

    initializer["initialize::out"] = randfixed["random_fixed_weight::input"];
    randfixed["random_fixed_weight::output"] = cm_encode["cm_encoder::input"];
    cm_encode["cm_encoder::output"] = cm_decode["cm_decoder::input"];
    cm_decode["cm_decoder::output"] = comp["compare::input1"];
    randfixed["random_fixed_weight::output"] = comp["compare::input2"];
    comp["compare::output"] = finalizer["finalize::in"];

    std::vector<runtime::Task *> first = {&initializer("initialize")};

    runtime::Sequence seq(first);

    std::ofstream file("graph.dot");
    seq.export_dot(file);

    for (auto lt : seq.get_tasks_per_types())
        for (auto t : lt)
        {
            t->set_stats(true);
            t->set_debug(false);
        }

    seq.exec_seq();

    const std::vector<int> &final_data = finalizer.get_final_data()[0];
    int error_sum = 0;

    for (const int &val : final_data)
        error_sum += val;

    cout << "Error sum : " << error_sum << endl;

    tools::Stats::show(seq.get_modules_per_types(), true, throughput);

    /* ************************************************************************* */
    /* ************************************************************************* */
}
