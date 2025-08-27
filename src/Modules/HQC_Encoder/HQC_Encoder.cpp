#include <iostream>
#include "Modules/HQC_Encoder/HQC_Encoder.hpp"

using namespace spu;
using namespace spu::module;
using namespace std;

HQC_Encoder::HQC_Encoder(HQC_public_key &PK, int k, int len, int r, int w) : Stateful(),
                                                                             frame_size(k),
                                                                             output_size1(PK.get_n()),
                                                                             output_size2(len),
                                                                             n1(PK.get_n1()),
                                                                             r(r),
                                                                             w(w)
{

    int output_size1_64t = (((output_size1) / (8)) + ((output_size1) % (8) == 0 ? 0 : 1));
    int output_size2_64t = (((output_size2) / (8)) + ((output_size2) % (8) == 0 ? 0 : 1));

    this->set_name("HQC_Encoder");
    this->set_short_name("HQC_Encoder");

    auto &t = create_task("hqc_encoder");
    auto input = create_socket_in<int>(t, "input", frame_size);
    auto output1 = create_socket_out<uint64_t>(t, "output1", output_size1_64t);
    auto output2 = create_socket_out<uint64_t>(t, "output2", output_size2_64t);

    this->create_codelet(t, [input, output1, output2, &PK](Module &m, runtime::Task &t, const size_t frame_id) -> int
                         {
	static_cast<HQC_Encoder&>(m).hqc_encoder(static_cast<int*>(t[input].get_dataptr()),
						 static_cast<uint64_t*>(t[output1].get_dataptr()),
						 static_cast<uint64_t*>(t[output2].get_dataptr()),
						 static_cast<HQC_public_key&>(PK),
						 frame_id);
	return 0; });
}

HQC_Encoder::~HQC_Encoder()
{
}

void HQC_Encoder::hqc_encoder(int *input, uint64_t *output1, uint64_t *output2, const HQC_public_key &PK,
                              const int frame_id)
{

    int k = (this->frame_size); /* Input size in bytes*/
    int m = GF2E::degree();

    /* Variables to store the results */
    GF2X res1(INIT_SIZE, this->output_size1);
    res1.SetLength(this->output_size1);
    
    GF2X res2(INIT_SIZE, this->output_size2);
    res2.SetLength(this->output_size2);

    GF2EX tmp_pol_m = uint8_vec_to_GF2EX(input, k); /* Conversion to FF_(2^8) polynomial with input bytes as coefficients */ 
    GF2X P = gf2x_set_cyclic(this->output_size1);

    /* Encoding */
    HQC_encoding(res1, res2, tmp_pol_m, PK.h, PK.s, PK.alpha, this->output_size1,
                      this->output_size2, this->n1, this->r, this->w, this->w, m, P);
    
    /* Conversion to int vectors */
    uint8_t output1_bytes[this->output_size1];
    BytesFromGF2X(output1_bytes, res1, this->output_size1);
    uint8_t output2_bytes[this->output_size2];
    BytesFromGF2X(output2_bytes, res2, this->output_size2);

    memcpy(output1, output1_bytes, this->output_size1);
    memcpy(output2, output2_bytes, this->output_size2);
}
