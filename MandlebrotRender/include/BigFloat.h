#include <cinttypes>
#include <vector>
#include <string>
#include <array>
#include <regex>

template<int length>
class BitString {
public:

	BitString() {}
	BitString(const std::array<uint32_t, length> b)
		: bits(b)
	{}
	BitString(int32_t c) {
		bits[length - 1] = c;
		for (int i = 0; i < length - 2; i++) {
			uint32_t v = c < 0 ? 0xFFFFFFFF : 0; //negative/positive with twos complement

		}
	}

	std::array<uint32_t, length> bits; //index 0 has most significant bits, index length-1 has least signficant bits

	bool add_with_carryover(int target_index, uint32_t x) {
		while (x != 0 && target_index >= 0) {
			//adding a bit in the 33rd spot that can be carried over during subtraction
			uint64_t sum = uint64_t(bits[target_index]) + uint64_t(x);
			bits[target_index] = uint32_t(sum);

			//carry over can only be a single 1
			x = (sum & 0x0000000100000000) != 0 ? 1 : 0;
			target_index--;
		}

		return x != 0;
	}

	bool operator+=(const BitString& b) {
		bool has_carry = false;
		for (int i = length - 1; i >= 0; i--) {
			has_carry = add_with_carryover(i, b.bits[i]);
		}
		return has_carry;
	}

	bool operator-=(const BitString& b) {
		return *this += -b;
	}

	bool operator+=(int32_t x) {
		if (x >= 0) return add_with_carryover(length - 1, x);
		else {
			//to properly represent negative numbers with twos compliment, need to have the leading Fs
			bool already_negative = is_negative();
			add_with_carryover(length - 1, x);
			for (int i = length - 2; i >= 0; i--) add_with_carryover(i, 0xFFFFFFFF);
			return already_negative && !is_negative(); //check if underflow occured
		}
	}

	bool operator-=(int32_t x) {
		return *this += -x;
	}

	BitString operator-() const {
		BitString out;
		for (int i = 0; i < length; i++) out.bits[i] = ~bits[i];
		out.add_with_carryover(length - 1, 1);
		return out;
	}

	int32_t sub_capped_at_int32_t_size(const BitString& b) const {
		BitString difference = *this;
		difference -= b;

		bool diff_negative = difference.is_negative();
		for (int i = 0; i < length - 1; i++) {
			if ((diff_negative && difference.bits[i] != 0xFFFFFFFF) || (!diff_negative && difference.bits[i] != 0)) {
				//difference greater than size of int32_t, return max negative/positive value
				return diff_negative ? 0x80000000 : 0x7FFFFFFF;
			}
		}
		return difference.bits[length - 1];
	}

	bool is_negative() const {
		return (bits[0] & 0x80000000) != 0;
	}

	void shift_right_with(uint32_t shift_amount, uint32_t shift_in = 0) {
		if (shift_amount == 0) return;

		int shift_indx_gap = shift_amount / 32;
		//simple case
		if (shift_amount % 32 == 0) {
			for (int i = length - 1; i >= 0; i--) {
				if (i >= shift_indx_gap)         bits[i] = bits[i - shift_indx_gap];
				else if (i == shift_indx_gap - 1) bits[i] = shift_in;
				else                              bits[i] = 0;
			}
			return;
		}

		uint32_t shift_in_lower_mask = 0xFFFFFFFF << (shift_amount % 32);
		uint32_t shift_in_upper_mask = ~shift_in_lower_mask;

		for (int i = length - 1; i >= 0; i--) {
			uint32_t upper_shifted_in_bits;
			uint32_t lower_shifted_in_bits;

			if (i >= shift_indx_gap + 1)  upper_shifted_in_bits = (shift_in_upper_mask & bits[i - shift_indx_gap - 1]) << (32 - (shift_amount % 32));
			else if (i == shift_indx_gap) upper_shifted_in_bits = (shift_in_upper_mask & shift_in) << (32 - (shift_amount % 32));
			else                          upper_shifted_in_bits = 0;

			if (i >= shift_indx_gap)          lower_shifted_in_bits = (shift_in_lower_mask & bits[i - shift_indx_gap]) >> (shift_amount % 32);
			else if (i == shift_indx_gap - 1) lower_shifted_in_bits = (shift_in_lower_mask & shift_in) >> (shift_amount % 32);
			else                              lower_shifted_in_bits = 0;

			bits[i] = upper_shifted_in_bits | lower_shifted_in_bits;
		}
	}

	void operator<<=(uint32_t shift_amount) {
		if (shift_amount == 0) return;

		int shift_indx_gap = shift_amount / 32;
		//simple case
		if (shift_amount % 32 == 0) {

			for (int i = 0; i < length - shift_indx_gap; i++) {
				if (i >= shift_indx_gap)         bits[i] = bits[i + shift_indx_gap];
				else                              bits[i] = 0;
			}
			return;
		}

		uint32_t shift_in_lower_mask = 0xFFFFFFFF >> (shift_amount % 32);
		uint32_t shift_in_upper_mask = ~shift_in_lower_mask;

		for (int i = 0; i < length - shift_indx_gap; i++) {
			uint32_t upper_shifted_in_bits;
			uint32_t lower_shifted_in_bits;

			if (i < length - shift_indx_gap - 1) upper_shifted_in_bits = (shift_in_upper_mask & bits[i + shift_indx_gap + 1]) >> (32 - (shift_amount % 32));
			else                                 upper_shifted_in_bits = 0;

			if (i < length - shift_indx_gap)  lower_shifted_in_bits = (shift_in_lower_mask & bits[i + shift_indx_gap]) << (shift_amount % 32);
			else                              lower_shifted_in_bits = 0;

			bits[i] = upper_shifted_in_bits | lower_shifted_in_bits;
		}
	}

	int32_t get_leftmost_set_bit() const {
		int i = 0;
		while (i < length && bits[i] == 0) i++;
		if (i == length) return -1;

		uint32_t chunk = bits[i];
		int chunk_pos = 0;
		while (chunk >>= 1) chunk_pos++;

		return (length - 1 - i) * 32 + chunk_pos;
	}

	bool is_zero() const {
		for (int i = 0; i < length; i++) if (bits[i] != 0) return false;
		return true;
	}

	void set_to_min_negative() {
		bits[0] = 0x80000000;
		for (int i = 1; i < length; i++) bits[i] = 0;
	}

	void set_to_max_positive() {
		bits[0] = 0x7FFFFFFF;
		for (int i = 1; i < length; i++) bits[i] = 0xFFFFFFFF;
	}

	bool operator<(const BitString<length>& b) const {
		int compare_indx = 0;
		//find most signficant uint32_t where this and b differ
		while (compare_indx < length && bits[compare_indx] == b.bits[compare_indx]) compare_indx++;

		if (compare_indx == length) return false;
		return bits[compare_indx] < b.bits[compare_indx];
	}

	bool operator>(const BitString<length>& b) const {
		int compare_indx = 0;
		//find most signficant uint32_t where this and b differ
		while (compare_indx < length&& bits[compare_indx] == b.bits[compare_indx]) compare_indx++;

		if (compare_indx == length) return false;
		return bits[compare_indx] > b.bits[compare_indx];
	}
};

template<int length>
bool operator>(const BitString<length>& a, uint32_t x);

template<int length>
bool operator<(const BitString<length>& a, uint32_t x);

struct DecimalString {
	std::vector<uint8_t> digits; //index 0 is leftmost, size - 1 rightmost

	bool is_zero() const {
		for (uint8_t i : digits) if (i != 0) return false;
		return true;
	}

	bool is_even() const {
		return digits.empty() || digits.back() % 2 == 0;
	}

	DecimalString divby2() const {
		DecimalString out{ {} };
		if (is_zero()) return out;

		//long division
		uint32_t x = 0;
		for (uint8_t d : digits) {
			x *= 10;
			x += d;

			if (x > 1) {
				out.digits.push_back(x / 2);
				x %= 2;
			}
			else {
				if (!out.is_zero()) out.digits.push_back(0); //avoid leading zeroes
			}
		}

		return out;
	}

	DecimalString multby2() const {
		DecimalString out{ {} };
		if (is_zero()) return out;

		uint8_t carry = 0;
		for (int i = digits.size() - 1; i >= 0; i--) {
			uint8_t prod = 2 * digits[i] + carry;
			out.digits.insert(out.digits.begin(), prod % 10);
			carry = prod / 10;
		}
		if (carry != 0) {
			out.digits.insert(out.digits.begin(), carry);
		}

		return out;
	}

	DecimalString add1() const {
		DecimalString out{ {} };

		uint8_t carry = 1;
		for (int i = digits.size() - 1; i >= 0; i--) {
			uint8_t sum = digits[i] + carry;
			out.digits.insert(out.digits.begin(), sum % 10);
			carry = sum / 10;
		}
		if (carry != 0) {
			out.digits.insert(out.digits.begin(), carry);
		}

		return out;
	}

	DecimalString multby2_keep_same_size(bool* had_carry_out) const {
		DecimalString out{ {} };
		if (is_zero()) return *this;

		uint8_t carry = 0;
		for (int i = digits.size() - 1; i >= 0; i--) {
			uint8_t prod = 2 * digits[i] + carry;
			out.digits.insert(out.digits.begin(), prod % 10);
			carry = prod / 10;
		}
		if (carry != 0) {
			*had_carry_out = true;
		}
		else {
			*had_carry_out = false;
		}

		return out;
	}

	std::vector<bool> to_binary() const {
		std::vector<bool> bits;

		DecimalString x{ digits };
		do {
			bits.push_back(!x.is_even());
			x = x.divby2();
		} while (!x.is_zero());

		return bits;
	}

	//instead of treating as an integer, treat as a number to the right of the floaing point
	std::vector<bool> convert_fraction_to_binary(int max_output_size) const {
		std::vector<bool> bits;

		DecimalString x{ digits };
		while (!x.is_zero() && bits.size() < max_output_size) {
			//multiply x by 2, checking if it overflowed a 1 to he left of fp
			bool had_carry;
			x = x.multby2_keep_same_size(&had_carry);
			//check if carry out increased the size
			bits.insert(bits.end(), had_carry);
		}

		return bits;
	}

	uint32_t to_int() const {
		uint32_t out = 0;
		double dec = 1;
		for (int i = digits.size() - 1; i >= 0; i--) {
			out += dec * digits[i];
			dec *= 10;
		}
		return out;
	}
};

struct FloatingPointBinaryString;

struct FloatingPointDecimalString {
	std::vector<uint8_t> left_of_fp;
	std::vector<uint8_t> right_of_fp;

	FloatingPointDecimalString operator+(const FloatingPointDecimalString& f) const {
		FloatingPointDecimalString out;

		uint8_t carry = 0;
		out.right_of_fp = std::vector<uint8_t>(std::max<int>(right_of_fp.size(), f.right_of_fp.size()));

		for (int i = out.right_of_fp.size() - 1; i >= 0; i--) {
			uint8_t sum = carry;
			if (i < right_of_fp.size()) sum += right_of_fp[i];
			if (i < f.right_of_fp.size()) sum += f.right_of_fp[i];

			out.right_of_fp[i] = sum % 10;
			carry = sum / 10;
		}

		out.left_of_fp = std::vector<uint8_t>(std::max<int>(left_of_fp.size(), f.left_of_fp.size()));
		for (int i = 0; i < out.left_of_fp.size(); i++) {
			uint8_t sum = carry;
			if (i < left_of_fp.size()) sum += left_of_fp[left_of_fp.size() - 1 - i];
			if (i < f.left_of_fp.size()) sum += f.left_of_fp[f.left_of_fp.size() - 1 - i];

			out.left_of_fp[out.left_of_fp.size() - 1 - i] = sum % 10;
			carry = sum / 10;
		}
		if (carry != 0) {
			out.left_of_fp.push_back(carry);
		}

		return out;
	}

	FloatingPointDecimalString times_digit(uint8_t d) const {
		FloatingPointDecimalString out{ {}, {} };

		uint8_t carry = 0;
		out.right_of_fp = std::vector<uint8_t>(right_of_fp.size());

		for (int i = out.right_of_fp.size() - 1; i >= 0; i--) {
			uint8_t prod = d * right_of_fp[i] + carry;

			out.right_of_fp[i] = prod % 10;
			carry = prod / 10;
		}

		out.left_of_fp = std::vector<uint8_t>(left_of_fp.size());
		for (int i = out.left_of_fp.size() - 1; i >= 0; i--) {
			uint8_t prod = d * left_of_fp[i] + carry;

			out.left_of_fp[i] = prod % 10;
			carry = prod / 10;
		}
		if (carry > 0) {
			out.left_of_fp.insert(out.left_of_fp.begin(), carry);
		}

		return out;
	}

	FloatingPointDecimalString operator*(const FloatingPointDecimalString f) {
		FloatingPointDecimalString out{ {}, {} };

		for (int i = 0; i < right_of_fp.size(); i++) {
			FloatingPointDecimalString prod = f.times_digit(right_of_fp[i]);
			prod.mult_self_by_pow10(-1 - i);
			out = out + prod;
		}
		for (int i = 0; i < left_of_fp.size(); i++) {
			FloatingPointDecimalString prod = f.times_digit(left_of_fp[left_of_fp.size() - 1 - i]);
			prod.mult_self_by_pow10(i);
			out = out + prod;
		}

		return out;
	}

	FloatingPointDecimalString div2() const {
		FloatingPointDecimalString out{ {}, {} };

		uint8_t sub_term = 0;
		for (int i = 0; i < left_of_fp.size(); i++) {
			sub_term *= 10;
			sub_term += left_of_fp[i];

			if (sub_term >= 2) {
				out.left_of_fp.push_back(sub_term / 2);
				sub_term %= 2;
			}
			else if (!out.left_of_fp.empty()) {
				out.left_of_fp.push_back(0);
			}
		}
		int i = 0;
		while (i < right_of_fp.size() || sub_term != 0) {
			sub_term *= 10;
			if (i < right_of_fp.size()) {
				sub_term += right_of_fp[i];
			}

			if (sub_term >= 2) {
				out.right_of_fp.push_back(sub_term / 2);
				sub_term %= 2;
			}
			else {
				out.right_of_fp.push_back(0);
			}

			i++;
		}

		return out;
	}

	FloatingPointBinaryString to_binary(int max_bits_right_of_fp) const;

	void mult_self_by_pow10(int i) {
		while (i > 0) {
			uint8_t shift_in_digit = 0;
			if (!right_of_fp.empty()) {
				shift_in_digit = right_of_fp.front();
				right_of_fp.erase(right_of_fp.begin());
			}
			left_of_fp.push_back(shift_in_digit);
			i--;
		}
		while (i < 0) {
			uint8_t shift_in_digit = 0;
			if (!left_of_fp.empty()) {
				shift_in_digit = left_of_fp.back();
				left_of_fp.pop_back();
			}
			right_of_fp.insert(right_of_fp.begin(), shift_in_digit);
			i++;
		}
	}

	bool is_zero() const {
		return left_of_fp_zero() && right_of_fp_zero();
	}

	bool left_of_fp_zero() const {
		for (uint8_t i : left_of_fp) if (i != 0) return false;
		return true;
	}

	bool right_of_fp_zero() const {
		for (uint8_t i : right_of_fp) if (i != 0) return false;
		return true;
	}
};

struct FloatingPointBinaryString {
	std::vector<bool> left_of_fp;
	std::vector<bool> right_of_fp;


	void mult_self_by_pow2(int i) {
		while (i > 0) {
			bool shift_in_bit = false;
			if (!right_of_fp.empty()) {
				shift_in_bit = right_of_fp.front();
				right_of_fp.erase(right_of_fp.begin());
			}
			left_of_fp.push_back(shift_in_bit);
			i--;
		}
		while (i < 0) {
			uint8_t shift_in_bit = false;
			if (!left_of_fp.empty()) {
				shift_in_bit = left_of_fp.back();
				left_of_fp.pop_back();
			}
			right_of_fp.insert(right_of_fp.begin(), shift_in_bit);
			i++;
		}
	}

	FloatingPointDecimalString to_decimal() const {
		FloatingPointDecimalString fractional_part = { {}, {} };
		FloatingPointDecimalString integer_part = { {}, {} };

		for (int i = 0; i < left_of_fp.size(); i++) {
			integer_part = integer_part.times_digit(2);
			if (left_of_fp[i]) {
				integer_part = integer_part + FloatingPointDecimalString{ {1}, {} };
			}
		}

		for (int i = right_of_fp.size() - 1; i >= 0; i--) {
			fractional_part = fractional_part.div2();
			if (right_of_fp[i]) {
				fractional_part = fractional_part + FloatingPointDecimalString{ {}, {5} };
			}
		}

		return FloatingPointDecimalString{ integer_part.left_of_fp, fractional_part.right_of_fp };
	}

};

template<int length>
bool operator<(const BitString<length>& a, const BitString<length>& b);
template<int length>
bool operator>(const BitString<length>& a, const BitString<length>& b);
template<int length>
bool operator<(const BitString<length>& a, uint32_t x);
template<int length>
bool operator>(const BitString<length>& a, uint32_t x);

template <int mantissa_length>
class BigFloat {
public:
	BigFloat()
		: is_zero(true)
	{
		mantissa.bits.fill(0);
	}

	BigFloat(double d)
		: is_zero(d == 0), is_negative(d < 0)
	{
		uint64_t d_bits = *(uint64_t*)&d;

		uint64_t d_mantissa = DOUBLE_MANTISSA_MASK & d_bits;
		exponent = ((DOUBLE_EXPONENT_MASK & d_bits) >> DOUBLE_EXPONENT_OFFSET) - DOUBLE_EXPONENT_BIAS;

		mantissa.bits.fill(0);
		//store the first 32 bits of the mantissa in length - 1, and the remaining 20 bits in length - 2
		mantissa.bits[0] = d_mantissa >> FIRST_MANTISSA_CHUNK_IN_DOUBLE_OFFSET;
		if (mantissa_length > 1) {
			mantissa.bits[1] = (d_mantissa & 0xFFFFF) << NUM_BITS_CUT_OFF_FROM_SECOND_MANITSSA_CHUNK_IN_DOUBLE;
		}
	}

	BigFloat(bool is_negative, int32_t exponent, BitString<mantissa_length>&& mantissa)
		: is_negative(is_negative), is_zero(false), exponent(exponent), mantissa(std::move(mantissa))
	{}

	static bool scientific_notation_number_valid(std::string scientific_notation) {
		// Optional leading "-"
		// at least 1 digit
		// optional "." for floating point
		// 0 or more digits right of the floating point, if it exists
		// optionally e or E, followed by either + or - and then at least one digit
		std::regex scientific_notation_pattern(R"rgx(\-?[0-9]+\.?[0-9]*(([eE]\+[0-9]+)|([eE]\-[0-9]+))?)rgx");
		return std::regex_match(scientific_notation, scientific_notation_pattern);
	}

	BigFloat(std::string scientific_notation)
		: is_negative(false), is_zero(false)
	{
		if (!scientific_notation_number_valid(scientific_notation)) {
			*this = zero();
			printf("WARN: %s is not valid scientific notation format.\n", scientific_notation.c_str());
			return;
		}
		mantissa.bits.fill(0);

		std::vector<uint8_t> left_of_fp;
		std::vector<uint8_t> right_of_fp;
		int32_t exponent_val = 0;
		bool exponent_negative = false;

		enum parse_state { OPTIONAL_NEGATIVE, LEFT_OF_FLOATING_POINT, RIGHT_OF_FLOATING_POINT, EXPONENT_SIGN, EXPONENT_DIGITS };
		parse_state state = OPTIONAL_NEGATIVE;
		int i = 0;
		while (i < scientific_notation.size()) {
			char c = scientific_notation[i];
			switch (state) {
			case OPTIONAL_NEGATIVE:
				if (c == '-') {
					is_negative = true;
					i++;
				}
				state = LEFT_OF_FLOATING_POINT;
				break;
			case LEFT_OF_FLOATING_POINT:
				if (c == '.') {
					state = RIGHT_OF_FLOATING_POINT;
				}
				else if (c == 'e' || c == 'E') {
					state = EXPONENT_SIGN;
				}
				else {
					left_of_fp.push_back(c - '0');
				}
				i++;
				break;
			case RIGHT_OF_FLOATING_POINT:
				if (c == 'e' || c == 'E') {
					state = EXPONENT_SIGN;
				}
				else {
					right_of_fp.push_back(c - '0');
				}
				i++;
				break;
			case EXPONENT_SIGN:
				if (c == '-') {
					exponent_negative = true;
				}
				i++;
				state = EXPONENT_DIGITS;
				break;
			case EXPONENT_DIGITS:
				uint32_t digit = c - '0';
				exponent_val *= 10;
				exponent_val += digit;
				i++;
				break;
			}
		}

		FloatingPointDecimalString digits = { left_of_fp, right_of_fp };

		if (digits.is_zero()) {
			*this = zero();
			return;
		}

		FloatingPointBinaryString binary = digits.to_binary(mantissa_length * 160); //allow for a lot of shifting if needed

		for (int i = 0; i < binary.right_of_fp.size() && i < 32 * mantissa_length; i++) {
			uint8_t bit_value = binary.right_of_fp[i];
			mantissa.bits[i / 32] |= ((bit_value & 0x1) << (31 - (i % 32)));
		}

		uint32_t val_left_of_fp = 0;
		for (uint8_t d : digits.left_of_fp) {
			val_left_of_fp = 10 * val_left_of_fp + d;
		}

		if (val_left_of_fp > 0) {
			int required_shift_right = 0;
			uint32_t x = val_left_of_fp;
			while (x >>= 1) required_shift_right++;
			mantissa.shift_right_with(required_shift_right, val_left_of_fp);
			exponent += required_shift_right;
		}
		else {
			int required_shift_left_amount = 32 * mantissa_length - mantissa.get_leftmost_set_bit();
			mantissa <<= required_shift_left_amount;
			exponent -= required_shift_left_amount;
		}

		//very very inefficient but this code doesn't need to run somewhere performant
		BigFloat ten(10);
		for (int i = 0; i < exponent_val; i++) {
			if (exponent_negative) *this = *this / ten;
			else *this = *this * ten;
		}
	}

	static BigFloat zero() {
		return BigFloat();
	}

	double to_double() const {
		if (is_zero) return 0;

		uint64_t out_exponent = uint64_t(exponent + DOUBLE_EXPONENT_BIAS) << DOUBLE_EXPONENT_OFFSET;
		uint64_t out_sign = is_negative ? 0x8000000000000000 : 0;
		uint64_t out_mantissa_bits;
		if (mantissa_length > 1) {
			uint32_t ls_mantissa_bits = mantissa.bits[1];
			uint32_t ms_mantissa_bits = mantissa.bits[0];
			out_mantissa_bits = uint64_t(ms_mantissa_bits) << FIRST_MANTISSA_CHUNK_IN_DOUBLE_OFFSET | (uint64_t(ls_mantissa_bits) >> NUM_BITS_CUT_OFF_FROM_SECOND_MANITSSA_CHUNK_IN_DOUBLE & 0xFFFFF);
		}
		else {
			out_mantissa_bits = uint64_t(mantissa.bits[0]) << FIRST_MANTISSA_CHUNK_IN_DOUBLE_OFFSET;
		}

		uint64_t out_bits = out_sign | out_exponent | out_mantissa_bits;
		return *(double*)&out_bits;
	}

	std::string to_scientific_notation() const {
		if (is_zero) return "0.0e+0";

		FloatingPointBinaryString binary = { { true }, {} };
		for (int i = 0; i < 32 * mantissa_length; i++) {
			bool ith_bit_set = (mantissa.bits[i / 32] & (0x1 << (31 - (i % 32)))) != 0;
			binary.right_of_fp.push_back(ith_bit_set);
		}

		binary.mult_self_by_pow2(exponent);

		FloatingPointDecimalString decimal = binary.to_decimal();

		int dec_exp_number = 0;
		if (decimal.left_of_fp.size() > 0) {
			dec_exp_number = decimal.left_of_fp.size() - 1;
		}
		else {
			//set dec_exp_number to 1 greater than index of first nonzero
			while (decimal.right_of_fp[dec_exp_number++] == 0);
			dec_exp_number *= -1;
		}
		decimal.mult_self_by_pow10(-dec_exp_number);

		std::string out;
		if (is_negative) out.push_back('-');
		for (uint8_t d : decimal.left_of_fp) out.push_back('0' + d);
		out.push_back('.');
		if (decimal.right_of_fp_zero()) out.push_back('0');
		else for (uint8_t d : decimal.right_of_fp) out.push_back('0' + d);
		out.push_back('e');
		if (dec_exp_number >= 0) out.push_back('+');
		out += std::to_string(dec_exp_number);

		return out;
	}

	template<int source_mant_len>
	BigFloat(const BigFloat<source_mant_len>& b)
		: is_negative(b.is_negative), is_zero(b.is_zero), exponent(b.exponent)
	{
		for (int i = 0; i < mantissa_length; i++) {
			mantissa.bits[i] = i < source_mant_len ? b.mantissa.bits[i] : 0;
		}
	}

	BigFloat operator-() const {
		BigFloat b(*this);
		b.is_negative = !b.is_negative;
		return b;
	}

	BigFloat operator+(const BigFloat& b) const {

		if (this->is_zero) return b;
		if (b.is_zero) return *this;

		BigFloat out;
		out.is_zero = false;
		int32_t exponent_difference = exponent - b.exponent;
		const BitString<mantissa_length>* larger_exponent_mantissa = nullptr;
		if (exponent_difference > 0) {
			out.exponent = exponent;
			out.is_negative = is_negative;
			out.mantissa = b.mantissa;
			larger_exponent_mantissa = &mantissa;
		}
		else {
			out.exponent = b.exponent;
			out.is_negative = b.is_negative;
			out.mantissa = mantissa;
			larger_exponent_mantissa = &b.mantissa;
		}

		//difference in scale so large that it won't change any values
		if (abs(exponent_difference) > mantissa_length * 32) return exponent_difference > 0 ? *this : b;

		out.mantissa.shift_right_with(abs(exponent_difference), 1); //shift in 1 in front of the floating point

		if (is_negative == b.is_negative) {
			uint32_t val_left_of_fp = exponent_difference == 0 ? 2 : 1; //from 1 left of mantissa
			bool carry_left_of_fp = out.mantissa += *larger_exponent_mantissa;
			if (carry_left_of_fp) {
				val_left_of_fp++;
			}
			if (val_left_of_fp >= 2) {
				//can only be 1, 2, or 3. so can get away with just sliding right once or none
				out.mantissa.shift_right_with(1, val_left_of_fp & 1);
				out.exponent += 1;
			}
		}
		else {
			int32_t val_left_of_fp = exponent_difference == 0 ? 0 : 1; //from 1 left of mantissa
			if (out.mantissa > *larger_exponent_mantissa) {
				val_left_of_fp--;
			}
			out.mantissa -= *larger_exponent_mantissa;
			if (val_left_of_fp == -1) {
				val_left_of_fp = 0;
				out.is_negative = !out.is_negative;
			}
			else {
				//subtraction order was backwards, so need to flip
				out.mantissa = -out.mantissa;
			}
			if (val_left_of_fp == 0) {
				int leftmost_set_bit = out.mantissa.get_leftmost_set_bit();
				if (leftmost_set_bit == -1) return zero();
				int shift_left_amount = 32 * mantissa_length - leftmost_set_bit;
				out.mantissa <<= shift_left_amount;
				out.exponent -= shift_left_amount;
			}
		}

		return out;
	}

	void operator+=(const BigFloat& b) {
		*this = *this + b;
	}

	BigFloat operator-(const BigFloat& b) const {
		return *this + -b;
	}

	BigFloat operator*(const BigFloat& b) const {
		if (is_zero || b.is_zero)  return zero();

		BigFloat out(*this);
		out.is_negative = is_negative ^ b.is_negative;
		out.exponent += b.exponent;

		//detect and cap at min and max exponent values in case of over/under flow
		if (exponent < 0 && b.exponent < 0 && out.exponent > 0) out.exponent = 0x8000000;
		if (this->exponent > 0 && b.exponent > 0 && out.exponent < 0) out.exponent = 0x77777777;

		uint32_t value_left_of_floating_point = 1; //start with 1 from leading 1 x 1

		if (!this->mantissa.is_zero() && !b.mantissa.is_zero())
			for (int i = 0; i < mantissa_length; i++) {
				int max_index_within_accuracy_bounds = std::min(mantissa_length, mantissa_length - i);
				for (int j = 0; j < max_index_within_accuracy_bounds; j++) {

					uint64_t product = uint64_t(b.mantissa.bits[i]) * uint64_t(mantissa.bits[j]);
					uint32_t lower_half = product & 0x00000000FFFFFFFF;
					uint32_t upper_half = product >> 32;

					int target_index = i + j + 1;
					if (target_index < mantissa_length) {
						value_left_of_floating_point += out.mantissa.add_with_carryover(target_index, lower_half);
					}
					value_left_of_floating_point += out.mantissa.add_with_carryover(target_index - 1, upper_half);
				}
			}

		//floating point representation as an implied 1 before the mantissa (e.g. 1.<mantissa bits>).
		//we already accounted for the mutliplication for <this> times the 1 from b by initializing out with the mantissa from <this>,
		//now we need to add the mantissa from b.
		if (out.mantissa += b.mantissa) {
			value_left_of_floating_point += 1;
		}

		if (value_left_of_floating_point == 0) return out;

		uint32_t carry_out_msb_pos = 1;
		//https://stackoverflow.com/a/4970859
		uint32_t v = value_left_of_floating_point;
		while (v >>= 1) {
			carry_out_msb_pos++;
		}

		uint32_t right_shift_amount = carry_out_msb_pos - 1;
		out.mantissa.shift_right_with(right_shift_amount, value_left_of_floating_point);
		out.exponent += right_shift_amount;

		return out;
	}

	BigFloat times2() const {
		BigFloat out(*this);
		out.exponent += 1;
		return out;
	}

	BigFloat operator/(const BigFloat& b) const {
		//subtract the exponents
		if (is_zero) return zero();
		else if (b.is_zero) {
			printf("ERR: divison by zero");
			return zero();
		}

		//find inverse of a floating point approximation of the mantissa of b
		uint64_t b_mantissa_float_approx = 0x3FF0000000000000 | (uint64_t(b.mantissa.bits[0]) << FIRST_MANTISSA_CHUNK_IN_DOUBLE_OFFSET);
		if (mantissa_length > 1) b_mantissa_float_approx |= (uint64_t(b.mantissa.bits[1]) >> NUM_BITS_CUT_OFF_FROM_SECOND_MANITSSA_CHUNK_IN_DOUBLE & 0xFFFFF);
		double inverse = 1.0f / *(double*)&b_mantissa_float_approx;

		BigFloat b_inv(inverse);
		BigFloat b_mantissa(b);
		b_mantissa.is_negative = false;
		b_mantissa.exponent = 0;

		//refine inverse newton-raphson

		//computing number of iterations needed knowing that error is halfed every time,
		//and that our approximation is already accurate to the first 52 bits.
		int required_iterations = 0;
		int x = (mantissa_length * 32) >> 5;
		while (x >>= 1) required_iterations++;

		for (int i = 0; i < required_iterations; i++) {
			BigFloat apple = b_inv.times2();
			BigFloat banana = b_inv * b_inv * b_mantissa;
			b_inv = b_inv.times2() - b_inv * b_inv * b_mantissa;
		}

		b_inv.exponent -= b.exponent;
		b_inv.is_negative = b.is_negative;

		double d = b_inv.to_double();
		//compute quoient
		return *this * b_inv;
	}

	bool get_is_negative() const {
		return is_negative;
	}

	BigFloat operator*(double d) const {
		if (d == 2) return this->times2();
		return BigFloat(d) * *this;
	}

	BigFloat operator/(double d) const {
		return *this / BigFloat(d);
	}

	bool operator<(double d) const {
		return this->to_double() < d;
	}

private:
	bool is_negative;
	bool is_zero;
	int32_t exponent;
	BitString<mantissa_length> mantissa;

	//https://en.wikipedia.org/wiki/Double-precision_floating-point_format
	const static uint32_t DOUBLE_EXPONENT_BIAS = 1023;
	const static uint64_t DOUBLE_MANTISSA_MASK = 0x000FFFFFFFFFFFFF;
	const static uint64_t DOUBLE_EXPONENT_MASK = 0x7FF0000000000000;
	const static uint32_t DOUBLE_EXPONENT_OFFSET = 52;
	const static uint32_t NUM_BITS_CUT_OFF_FROM_SECOND_MANITSSA_CHUNK_IN_DOUBLE = 12;
	const static uint32_t FIRST_MANTISSA_CHUNK_IN_DOUBLE_OFFSET = DOUBLE_EXPONENT_OFFSET - 32;

	template <int other_mant_l>
	friend class BigFloat;
};