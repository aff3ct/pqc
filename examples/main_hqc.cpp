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



    /* ************************************************************************* */
    /*                                 TESTS FOR HQC                             */
    /* ************************************************************************* */

    FLINT_TEST_INIT(state); /* for randomness */

    slong m = 8;

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
                std::cerr << "  -b, --throughput        "
                          << "Show throughput stats     "
                          << std::endl;
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


   std::cout << "HQC Level: " << level << std::endl;

    GF2X P; // Context mod polynomial

    // Context for GF2E ~ Constant for HQC
    // X^8 + X^4 + X^3 X^2 + 1
    SetCoeff(P, 8);
    SetCoeff(P, 4);
    SetCoeff(P, 3);
    SetCoeff(P, 2);
    SetCoeff(P, 0);
    GF2E::init(P);

    int k_bytes, n1, n2, r, n, w, we, wr;
    HQC_params(k_bytes, n1, n2, r, n, w, we, wr, level);

    int len = r * n1 * n2;

    const int FRAME_SIZE = k_bytes;
    const int GLOBAL_LENGTH = n;
    const int WEIGHT = w;
    const int ERROR_WEIGHT = we;
    const int N1 = n1;
    const int CODE_LENGTH = len;
    const int R = r;

    HQC_secret_key SK = HQC_secret_key(GLOBAL_LENGTH);
    HQC_public_key PK = HQC_public_key(GLOBAL_LENGTH, N1);
    HQC_keygen_naive(SK, PK, WEIGHT, state);

    module::Initializer<int> initializer(FRAME_SIZE);
    module::Finalizer<int> finalizer(FRAME_SIZE);

    module::Comparator comp(FRAME_SIZE);
    module::RandomVector random_vector(FRAME_SIZE);
    module::HQC_Encoder hqc_encode(PK, FRAME_SIZE, CODE_LENGTH, R, ERROR_WEIGHT);
    module::HQC_Decoder hqc_decode(SK, PK, FRAME_SIZE, CODE_LENGTH, R);

    initializer["initialize::out"] = random_vector["random_vector::input"];
    random_vector["random_vector::output"] = hqc_encode["hqc_encoder::input"];
    hqc_encode["hqc_encoder::output1"] = hqc_decode["hqc_decoder::input1"];
    hqc_encode["hqc_encoder::output2"] = hqc_decode["hqc_decoder::input2"];

    hqc_decode["hqc_decoder::output"] = comp["compare::input1"];
    random_vector["random_vector::output"] = comp["compare::input2"];
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

    int error = 0;
    for (int i = 0; i < 1; i++)
    {

        seq.exec_seq();

        const std::vector<int> &final_data = finalizer.get_final_data()[0];
        for (const int &val : final_data)
        {
            if (val != 0)
            {
                error++;
            }
        }
    }

    printf("Error:%i\n", error);

    tools::Stats::show(seq.get_modules_per_types(), true, throughput);

    /* ************************************************************************* */
    /* ************************************************************************* */
}
