#include "../include/BigFloat.h"

template<int length>
bool operator>(const BitString<length>& a, uint32_t x) {
	if (a.bits[0] == x) return length > 1 && a.bits[1] != 0;
	return a.bits[0] > x;
}

template<int length>
bool operator<(const BitString<length>& a, uint32_t x) {
	return a.bits[0] < x;
}

FloatingPointBinaryString FloatingPointDecimalString::to_binary(int max_bits_right_of_fp) const {
	FloatingPointBinaryString out;

	FloatingPointDecimalString x{ left_of_fp, {} };
	while (!x.left_of_fp_zero()) {
		out.left_of_fp.insert(out.left_of_fp.begin(), x.left_of_fp.back() % 2 != 0);
		x = x.div2();
	}

	x = FloatingPointDecimalString{ {}, right_of_fp };
	while (!x.right_of_fp_zero() && out.right_of_fp.size() < max_bits_right_of_fp) {
		//multiply x by 2, checking if it overflowed a 1 to he left of fp
		x = x.times_digit(2);
		bool had_carry = !x.left_of_fp_zero();
		x.left_of_fp = {};
		//check if carry out increased the size
		out.right_of_fp.insert(out.right_of_fp.end(), had_carry);
	}

	return out;
}

template<int exp_len, int mant_len>
BigFloat<mant_len> operator*(double d, BigFloat<mant_len> b) {
	if (d == 2) return b.times2();
	return BigFloat<exp_len, mant_len>(d) * b;
}