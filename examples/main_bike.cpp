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

    slong m = 8;

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


    std::cout << "BIKE Level: " << level << std::endl;

    int weight = 0;
    int error_weight = 0;
    int NbIter = 0;
    int tau = 0;

    Bike_params(r, weight, error_weight, level);
    BGF_params(NbIter, tau, level);

    int len = 2 * r;
    int n = len;

    const int FRAME_SIZE = len;
    const int SYND_SIZE = r;
    const int WEIGHT = weight;
    const int ERROR_WEIGHT = error_weight;
    const int NBITER = NbIter;
    const int TAU = tau;

    Bike_secret_key SK = Bike_secret_key(SYND_SIZE);
    Bike_public_key PK = Bike_public_key(SYND_SIZE);
    Bike_keygen_naive(SK, PK, weight);

    module::Initializer<int> initializer(FRAME_SIZE);
    module::Incrementer<int> incr1(FRAME_SIZE);
    module::Finalizer<int> finalizer(FRAME_SIZE);

    module::Comparator comp(FRAME_SIZE);
    module::Bike_RandomFixedWeight randfixed(FRAME_SIZE, WEIGHT);
    module::Bike_Encoder bike_encode(PK);
    module::Bike_Decoder bike_decode(WEIGHT, NBITER, TAU, SK);

    initializer["initialize::out"] = randfixed["random_fixed_weight::input"];
    randfixed["random_fixed_weight::output"] = bike_encode["bike_encoder::input"];
    bike_encode["bike_encoder::output"] = bike_decode["bike_decoder::input"];
    bike_decode["bike_decoder::output"] = comp["compare::input1"];
    randfixed["random_fixed_weight::output"] = comp["compare::input2"];
    comp["compare::output"] = finalizer["finalize::in"];

    std::vector<runtime::Task *> first = {&initializer("initialize")};

    runtime::Sequence seq(first);

    std::ofstream file("graph.dot");
    seq.export_dot(file);
    int error = 0;
    for (auto lt : seq.get_tasks_per_types())
        for (auto t : lt)
        {
            t->set_stats(true);
            t->set_debug(false);
        }

    for (int i = 0; i < 10; i++)
    {
        auto t_start = std::chrono::steady_clock::now();

        seq.exec_seq();

        std::chrono::nanoseconds duration = std::chrono::steady_clock::now() - t_start;

        auto elapsed_time = duration.count() / 1000.f / 1000.f;
        const std::vector<int> &final_data = finalizer.get_final_data()[0];

        for (const int &val : final_data)
        {
            //std::cout << val << ",";
            if (val != 0)
            {
                error++;
            }
        }
        //std::cout << std::endl;

    }
    std::cout << "Error sum  " << error << std::endl;

    tools::Stats::show(seq.get_modules_per_types(), true, throughput);

    /* ************************************************************************* */
    /* ************************************************************************* */
}
