// MandlebrotRender.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <cstdio>
#include <vector>
#include <cinttypes>
#include "../include/BigFloat.h"

template<int exp_l, int mnt_l>
void test(double x, double y) {
    BigFloat<exp_l, mnt_l> big_x(x);
    BigFloat <exp_l, mnt_l> big_y(y);

    BigFloat<exp_l, mnt_l> big_product = big_x * big_y;
    BigFloat<exp_l, mnt_l> big_sum = big_x + big_y;
    BigFloat<exp_l, mnt_l> big_diff = big_x - big_y;
    BigFloat<exp_l, mnt_l> big_quotient = big_x / big_y;

    printf(
        "\n=========\nTESTING FOR EXP_L: %d, MNT_L: %d, X=%.16f, Y=%.16f\nbig_x: %.16f\nbig_y: %.16f\nbig_product: %.16f, actual: %.16f\nbig_sum: %.16f, acual: %.16f\nbig_diff: %.16f, actual: %.16f\nbig_quotient: %.16f, actual: %.16f\n=========",
        exp_l, mnt_l, x, y, big_x.to_double(), big_y.to_double(), big_product.to_double(), x * y, big_sum.to_double(), x + y, big_diff.to_double(), x - y, big_quotient.to_double(), x / y
    );
}

template<int exp_l, int mnt_l>
void validate_on_many() {
    std::vector<double> some_random_doubles = { 0, 1, 3, 123123412, -0.000012312, -132123312.2, 44123123, -3.33333, 0.1, -64, 89, 3432413212.1, -0.332132, 0.5531431, -55123.32313254, 5012312.23123, 2, -4, 7, 8 };
    for (double x : some_random_doubles) {
        for (double y : some_random_doubles) {
            BigFloat<exp_l, mnt_l> big_x(x);
            BigFloat <exp_l, mnt_l> big_y(y);

            BigFloat<exp_l, mnt_l> big_product = big_x * big_y;
            BigFloat<exp_l, mnt_l> big_sum = big_x + big_y;
            BigFloat<exp_l, mnt_l> big_diff = big_x - big_y;
            BigFloat<exp_l, mnt_l> big_quotient = y == 0? BigFloat<exp_l, mnt_l>::zero() : big_x / big_y;

            double eps = 0.00000001; //eventually shouldn't need
            if (abs(big_product.to_double() - (x * y)) > eps) printf("ERROR ON MULT %.20f * %.20f. big_calc: %.20f, actual: %.20f\n", x, y, big_product.to_double(), x * y);
            if (abs(big_sum.to_double() - (x + y)) > eps) printf("ERROR ON ADD %.20f + %.20f. big_calc: %.20f, actual: %.20f\n", x, y, big_sum.to_double(), x + y);
            if (abs(big_diff.to_double() - (x - y)) > eps) printf("ERROR ON SUB %.20f - %.20f. big_calc: %.20f, actual: %.20f\n", x, y, big_diff.to_double(), x - y);
            //if (y != 0 && abs(big_quotient.to_double() - (x / y)) > eps) printf("ERROR ON DIV %.14f / %.14f. big_calc: %.14f, actual: %.14f\n", x, y, big_quotient.to_double(), x / y);
        }
    }
}

//int main() {
//    ////test<3, 3>(0.11111111111111111111, 0.11111111111111111111);
//    //test<3, 3>(-1.5, -0.75);
//    ///*test<2, 2>(7, 7);
//    //test<1, 1>(0.5, 0.5);
//    //test<1, 2>(64, 128);
//    //test<3, 3>(1.5, 0.75);
//    //test<3, 3>(-1.5, 0.75);
//    //test<3, 3>(1.5, -0.75);
//    //test<3, 3>(-1.5, -0.75);*/
//    //test<3, 6>(3.1415926535, 312123.319412);
//    //double w = 1.0 / 0.0003154;
//    //float a = 1.0f / float(0.0003154);
//    //BigFloat<3, 3> u(1.0 / 0.0003154);
//    //BigFloat<3, 3> b(1.0f / float(0.0003154));
//    //test<3, 6>(0.00001, 0.0003154);
//
//
//    ///*validate_on_many<1, 2>();
//    //validate_on_many<2, 2>();
//    //validate_on_many<3, 2>();
//    //validate_on_many<2, 3>();
//    //validate_on_many<15, 15>();*/
//    ////validate_on_many<1, 100>();
//
//    ////printf("%.20f, %.20f\n", 0.5, BigFloat<2, 2>(0.5).to_double());
//
//    ////{
//    ////    float one = 1.0f;
//    ////    uint32_t bits = 0x3FFFFFFF;
//    ////    float f = *(float*)&bits;
//    ////    float inv = 1.0f / f;
//    ////    uint32_t inv_bits = *(uint32_t*)&inv;
//    ////    int inv_exponent = (0x7F100000 & inv_bits) >> 23 - 127;
//
//    ////    printf("one: %x, f: %f, inv: %f, inv_exponent: %d\n", *(uint32_t*)&one, f, inv, inv_exponent);
//    ////}
//
//    struct scientific_actual_pair {
//        std::string scientific;
//        double actual;
//    };
//
//    std::vector<scientific_actual_pair> test_pairs = {
//        scientific_actual_pair{"1.5e+1", 15 },
//        scientific_actual_pair{"0.15e-0", 0.15 },
//        scientific_actual_pair{"1.0e+1", 10 },
//        scientific_actual_pair{"1", 1},
//        scientific_actual_pair{"1.0", 1 },
//        scientific_actual_pair{"1.0e-1", 0.1 },
//        scientific_actual_pair{"1.5", 1.5 },
//        scientific_actual_pair{"31.415926535e-1", 3.1415926535 },
//        scientific_actual_pair{"31.415926535e+3", 31415.926535 },
//    };
//
//    for (scientific_actual_pair p : test_pairs) {
//        printf("scientific: %s, actual: %.16f, bigfloat: %s\n", p.scientific.c_str(), p.actual, BigFloat<1, 2>(p.scientific).to_scientific_notation().c_str());
//    }
//}
