#include <assert.h>
#include <iostream>

struct LogGuard {
    LogGuard(std::ostream& os, const char* prefix, const char* file, int lineno)
        : pos(&os) {
        if (prefix) *pos << prefix;
        *pos << file << ":" << lineno << ": ";
    }
    ~LogGuard() { *pos << std::endl; }
    std::ostream* pos;
};
const char* const INFO = nullptr;
const char* const WARN = "WARN: ";
const char* const ERROR = "ERR: ";
#define LOG(level) (*LogGuard(std::cout, level, __FILE__, __LINE__).pos)
#define CHECK(expr) if (!(expr)) LOG("CHECK: ") << #expr " is false. "

template <int E, int M> struct SimFloat;
template <int E, int M> float ToInfAndNaN(const SimFloat<E, M>&) { return 0; }
template <int E, int M> SimFloat<E, M> FromInfAndNaN(float) { return SimFloat<E, M>::ZERO; }

template <int E, int M>
struct SimFloat {
    static_assert(E+M <= 15, "only support low-bit floats");
    static const int EBIAS = (1 << (E-1)) - 1;

    SimFloat(bool sign, int exponent, int mantissa, bool overflow = false, bool underflow = false)
        : sign(sign), overflow_bit(overflow), underflow_bit(underflow), exponent(exponent), mantissa(mantissa) {
        CHECK(mantissa >= 0 && mantissa < (1<<M)) << "mantissa=" << mantissa;
        CHECK(exponent >= 0 && exponent < (1<<E)) << "exponent=" << exponent << " E=" << E;
    }

    static const int EBITS = E;
    static const int MBITS = M;
    static const SimFloat<E, M> INF;
    static const SimFloat<E, M> NANVAL;
    static const SimFloat<E, M> NORMAL_MAX;
    static const SimFloat<E, M> NORMAL_MIN;
    static const SimFloat<E, M> DENORM_MAX;
    static const SimFloat<E, M> DENORM_MIN;
    static const SimFloat<E, M> ZERO;

    // bitwise conversion to float32
    float ToFloat32() const {
        if (exponent != 0) {
            const float sp = ToInfAndNaN(*this);
            if (sp != 0) return sp;
            const auto new_exponent = exponent - (1 << (E-1)) + 128;
            const auto new_mantissa = (mantissa << (23 - M));
            const uint32_t enc = (sign ? 0x80000000 : 0) | (new_exponent << 23) | new_mantissa;
            return *(const float*)&enc;
        } else { // denorms
            if (mantissa == 0) {
                return (sign ? -0 : 0);
            }
            auto new_mantissa = mantissa;
            int exponent = 1;
            if (new_mantissa < (1<<M)) {
                auto offset = __builtin_clz(new_mantissa) - 31 + M;
                new_mantissa = (new_mantissa << offset) - (1<<M);
                exponent -= offset;
            }
            new_mantissa = new_mantissa << (23 - M);
            auto new_exponent = exponent - (1 << (E-1)) + 128;
            const uint32_t enc = (sign ? 0x80000000 : 0) | (new_exponent << 23) | new_mantissa;
            return *(const float*)&enc;
        }
    }

    float ToFloat32Slow() const {
        float val;
        if (exponent == 0) {
            val = mantissa * exp2f(1 - EBIAS - M);
        } else {
            const float sp = ToInfAndNaN(*this);
            if (sp != 0) return sp;
            val = (1 + float(mantissa) / (1<<M)) * exp2f((int)exponent - EBIAS);
        }
        return sign ? -val : val;
    }

    static SimFloat FromFloat32(float fval) {
        const uint32_t fenc = *(const uint32_t*)&fval;
        const bool sign = (fenc & 0x80000000);
        if (isinf(fval)) {
            auto v = on_overflow();
            return SimFloat(sign, v.exponent, v.mantissa, true);
        }
        if (isnan(fval)) {
            if (has_nan()) {
                return SimFloat(sign, NANVAL.exponent, NANVAL.mantissa);
            } else {
                LOG(ERROR) << "Not support NaN!";
                return SimFloat(sign, 0, 0);
            }
        }
        const int exponent = (fenc >> 23) & 0xFF;
        const int mantissa = fenc & 0x7FFFFF;
        if (exponent == 0) {
            // no way to hold denorms of float32
            return SimFloat(sign, 0, 0);
        }
        int new_exponent = exponent - 127 + EBIAS;
        int new_mantissa = mantissa + (1 << 23);
        return SimFloat::_Pack(sign, new_exponent, new_mantissa, 23-M);
    }

    void PrintEncoding(std::ostream& os) const {
        os << (sign ? '1' : '0') << '_';
        for (int i = E-1; i >= 0; --i) {
            os << ((exponent & (1 << i)) ? '1' : '0');
        }
        if (M > 0) {
            os << '_';
            for (int i = M-1; i >= 0; --i) {
                os << ((mantissa & (1 << i)) ? '1' : '0');
            }
        }
        if (overflow_bit) { os << "_OF"; }
        if (underflow_bit) { os << "_UF"; }
    }

    static int RShiftMantissa(int mantissa, int shift) {
        if (shift <= 0) return (mantissa << -shift);
        if (shift >= 32) return 0;
        const bool carry = (mantissa & (1 << (shift-1)));
        mantissa >>= shift;
        if (carry) ++mantissa;
        return mantissa;
    }

    SimFloat Multiply1(const SimFloat& rhs) const {
        // TODO: Handling Inf/Nan efficiently.
        const auto new_sign = sign ^ rhs.sign;
        auto new_exponent = exponent_no_bias() + rhs.exponent_no_bias() + EBIAS;
        auto new_mantissa = fraction() * rhs.fraction();
        if (new_mantissa == 0) return SimFloat(new_sign, 0, 0);
        const auto dz = 32 - __builtin_clz(new_mantissa) - (2*M+1);
        new_exponent += dz;
        return SimFloat::_Pack(new_sign, new_exponent, new_mantissa, M+dz);
    };

    SimFloat operator*(const SimFloat& rhs) const {
        return FromFloat32(ToFloat32() * rhs.ToFloat32());
    }

    SimFloat operator+(const SimFloat& rhs) const {
        return FromFloat32(ToFloat32() + rhs.ToFloat32());
    }

    SimFloat operator-(const SimFloat& rhs) const {
        return FromFloat32(ToFloat32() - rhs.ToFloat32());
    }

    SimFloat operator/(const SimFloat& rhs) const {
        return FromFloat32(ToFloat32() / rhs.ToFloat32());
    }

    bool operator==(const SimFloat& rhs) const {
        return sign == rhs.sign && exponent == rhs.exponent && mantissa == rhs.mantissa;
    }
    bool operator!=(const SimFloat& rhs) const { return !operator==(rhs); }

    static SimFloat _Pack(bool sign, int exponent, int mantissa, int rshift) {
        bool overflow = false;
        bool underflow = false;
        if (exponent >= 1) {
            mantissa = RShiftMantissa(mantissa, rshift);
            if (mantissa >= (1<<(M+1))) {
                // the rounding in rshift could result in extra carrying.
                mantissa >>= 1;
                exponent ++;
            }
            mantissa -= (1<<M);
        } else {                
            mantissa = RShiftMantissa(mantissa, rshift+(1-exponent));
            exponent = 0;
            if (mantissa >= (1<<M)) {
                // when exponent=0, a 2M+1 bit number rshift M+1 bits could result in a M+1 bit number.
                mantissa -= (1<<M);
                exponent = 1;
            }
            if (exponent == 0 && mantissa == 0) {
                underflow = true;
            }
        }
        if (exponent > NORMAL_MAX.exponent) {
            //LOG(WARN) << "overflow exponent=" << exponent << " max=" << (1<<E)-1;
            overflow = true;
            auto v = on_overflow();
            exponent = v.exponent;
            mantissa = v.mantissa;
        }
        return SimFloat(sign, exponent, mantissa, overflow, underflow);

    }

    static void PrintStats(std::ostream& os) {
        os << type_desc() << " EBIAS=" << EBIAS
           << " NormalMax=" << NORMAL_MAX
           << " NormalMin=" << NORMAL_MIN
           << " DenormMax=" << DENORM_MAX
           << " DenormMin=" << DENORM_MIN
           << " Inf=";
        if (has_inf()) {
            os << INF;
        } else {
            os << "none";
        }
        os << " NaN=";
        if (has_nan()) {
            os << NANVAL;
        } else {
            os << "none";
        }
    }
    
    static void ListAllPositiveValues(std::ostream& os) {
        int counter = 0;
        for (int i = 0; i < (1<<E); ++i) {
            for (int j = 0; j < (1<<M); ++j) {
                SimFloat f1(0, i, j);
                os << counter++ << ") ";
                f1.PrintEncoding(os);
                auto v1 = f1.ToFloat32();
                auto v2 = f1.ToFloat32Slow();
                os << " -> " << v1;
                if (v2 != v1 && !isnan(v2) && !isnan(v1)) {
                    os << " UNMATCHED:" << v2;
                }
                os << '\n';
            }
        }
    }

    int exponent_no_bias() const {
        return exponent ? (exponent - EBIAS) : (1 - EBIAS);
    }

    int fraction() const {
        return exponent ? ((1<<M) + mantissa) : mantissa;
    }

    bool is_zero() const { return exponent == 0 && mantissa == 0; }

    static bool has_inf() { return !INF.is_zero(); }
    static SimFloat on_overflow() { return has_inf() ? INF : NORMAL_MAX; }

    static bool has_nan() { return !NANVAL.is_zero(); }

    static std::string type_desc() {
        char desc_buf_[16];
        snprintf(desc_buf_, sizeof(desc_buf_), "FP%d-E%dM%d", E+M+1, E, M);
        return desc_buf_;
    }

public:
    bool sign;
    bool overflow_bit;
    bool underflow_bit;
    int exponent;
    int mantissa;
};

template <int E, int M>
const SimFloat<E, M> SimFloat<E, M>::INF(0, 0, 0);

template <int E, int M>
const SimFloat<E, M> SimFloat<E, M>::NANVAL(0, 0, 0); // 'NAN' was occupied.

template <int E, int M>
const SimFloat<E, M> SimFloat<E, M>::NORMAL_MAX(0, (1<<E)-1, (1<<M)-1);

template <int E, int M>
const SimFloat<E, M> SimFloat<E, M>::NORMAL_MIN(0, 1, 0);

template <int E, int M>
const SimFloat<E, M> SimFloat<E, M>::DENORM_MAX(0, 0, (1<<M)-1);

template <int E, int M>
const SimFloat<E, M> SimFloat<E, M>::DENORM_MIN(0, 0, std::min(M, 1));

template <int E, int M>
const SimFloat<E, M> SimFloat<E, M>::ZERO(0, 0, 0);

// Specialize E5M2
template <>
const SimFloat<5, 2> SimFloat<5, 2>::INF(0, (1<<5)-1, 0);

template <>
const SimFloat<5, 2> SimFloat<5, 2>::NANVAL(0, (1<<5)-1, 1); 

template <>
const SimFloat<5, 2> SimFloat<5, 2>::NORMAL_MAX(0, (1<<5)-2, (1<<2)-1);

template <> float ToInfAndNaN(const SimFloat<5, 2>& sf) {
    if (sf.exponent == SimFloat<5, 2>::INF.exponent) {
        if (sf.mantissa == SimFloat<5, 2>::INF.mantissa) {
            auto inf = std::numeric_limits<float>::infinity();
            return sf.sign ? -inf : inf;
        }
        return std::numeric_limits<float>::signaling_NaN();
    }
    return 0;
}

// Specialize E4M3
template <>
const SimFloat<4, 3> SimFloat<4, 3>::NANVAL(0, (1<<4)-1, (1<<3)-1); 

template <>
const SimFloat<4, 3> SimFloat<4, 3>::NORMAL_MAX(0, (1<<4)-1, (1<<3)-2);

template <> float ToInfAndNaN(const SimFloat<4, 3>& sf) {
    if (sf.exponent == ((1<<4)-1) && sf.mantissa == ((1<<3)-1)) {
        return std::numeric_limits<float>::signaling_NaN();
    }
    return 0;
}


template <int E, int M>
std::ostream &operator<<(std::ostream &os, const SimFloat<E, M> &f) {
    f.PrintEncoding(os);
    os << "(" << f.ToFloat32() << ")";
    return os;
}

static std::pair<float, float> RandN() {
    // using box-muller method
    const float u = rand() / (float)(RAND_MAX - 1); // not strictly uniform, enough for testing.
    const float v = rand() / (float)(RAND_MAX - 1);
    const float radius = sqrtf(-2*logf(u));
    const float theta = 2*M_PI*v;
    const auto r1 = radius * cosf(theta);
    const auto r2 = radius * sinf(theta);
    return {r1, r2};
}

static std::pair<float, float> TruncatedRandN(float multiple_sd) {
    CHECK(multiple_sd >= 3) << "Too small multiple of the standard deviation (may causing many retries)";
    int max_try = 0;
    std::pair<float, float> p;
    do {
        p = RandN();
        if (abs(p.first) <= multiple_sd && abs(p.second) <= multiple_sd) {
            return p;
        }
    } while (max_try++ < 5);
    if (std::abs(p.first) > multiple_sd)
        p.first = (p.first > 0 ? multiple_sd : -multiple_sd);
    if (std::abs(p.second) > multiple_sd)
        p.second = (p.second > 0 ? multiple_sd : -multiple_sd);
    return p;
}

typedef SimFloat<2, 1> FP4_E2M1;
typedef SimFloat<3, 0> FP4_E3M0;
typedef SimFloat<3, 2> FP6_E3M2;
typedef SimFloat<2, 3> FP6_E2M3;
typedef SimFloat<5, 2> FP8_E5M2;
typedef SimFloat<4, 3> FP8_E4M3;

struct AddOp {
    template <typename T> T operator()(T a, T b) { return a + b; };
    static const char DESC = '+';
};
struct MulOp {
    template <typename T> T operator()(T a, T b) { return a * b; };
    static const char DESC = '*';
};
struct DivOp {
    template <typename T> T operator()(T a, T b) { return a / b; };
    static const char DESC = '/';
};

template <typename FPX, typename Op>
void CheckOpPrecisions(int op_run, bool show_details = false) {
    Op op;
    
    double diff_pct_sum = 0;
    int pct_count = 0;
    for (int i = 0; i < op_run; ++i) {
        auto p = TruncatedRandN(3);
        auto a1 = FPX::FromFloat32(p.first);
        auto b1 = FPX::FromFloat32(p.second);
        auto c = op(a1, b1);
        //const auto exact_c1 = a1.ToFloat32() * b1.ToFloat32();
        const auto exact_c2 = op(p.first, p.second);
        const auto diff = c.ToFloat32() - exact_c2;
        const auto diff_pct = diff ? diff * 100 / exact_c2 : 0;
        if (!c.overflow_bit && std::abs(exact_c2) >= std::max(FPX::DENORM_MIN.ToFloat32(), 0.05f)) {
            diff_pct_sum += std::abs(diff_pct);
            ++pct_count;
        }
        if (show_details) {
            LOG(INFO) << a1 << op.DESC << b1 << "=" << c << " exact=" << exact_c2 << "(" << c.ToFloat32() - exact_c2 << ", " << diff_pct << "%)";
        }
    }
    LOG(INFO) << FPX::type_desc() << "(" << op.DESC 
              << ") avg_diff_pct=" << diff_pct_sum / pct_count << "% (count=" << pct_count << ")";
}

template <typename FPX>
void CheckDotProducts(int dp_run) {
    double total_abs_diff_sum = 0;
    for (int run = 0; run < dp_run; ++run) {
        double exact_prod_sum = 0;
        double lowbit_prod_sum = 0;
        int count = 0;
        const int K = 1024;
        for (int i = 0; i < K; ++i) {
            auto p = RandN();
            auto a1 = FPX::FromFloat32(p.first);
            auto b1 = FPX::FromFloat32(p.second);
            auto c = a1 * b1;
            const auto exact_c2 = p.first * p.second;
            if (!c.overflow_bit) {
                lowbit_prod_sum += c.ToFloat32();
                exact_prod_sum += exact_c2;
                ++count;
            }
        }
        auto prod_sum_diff = lowbit_prod_sum - exact_prod_sum;
        LOG(INFO) << FPX::type_desc() << "(DP_" << run << ") K=" << K << " exact_prod_sum=" << exact_prod_sum << " lowbit_prod_sum=" << lowbit_prod_sum
                  << " diff=" << prod_sum_diff << " diff_pct=" << prod_sum_diff * 100 / exact_prod_sum << "% count=" << count;
        total_abs_diff_sum += std::abs(prod_sum_diff);
    }
    LOG(INFO) << "====== avg_abs_diff=" << total_abs_diff_sum / dp_run;
}

int Tests(int argc, char* argv[]);

int main(int argc, char* argv[]) {
    srand(time(0));

    if (argc >= 2 && strcmp(argv[1], "tests") == 0)
        return Tests(argc - 2, argv + 2);

    const int N = 1000;
    const bool show_details = (argc >= 2 && strcmp(argv[1], "verbose") == 0);

    CheckOpPrecisions<FP8_E5M2, MulOp>(N, show_details);
    CheckOpPrecisions<FP8_E5M2, AddOp>(N, show_details);
    CheckOpPrecisions<FP8_E4M3, MulOp>(N, show_details);
    CheckOpPrecisions<FP8_E4M3, AddOp>(N, show_details);
    CheckOpPrecisions<FP6_E3M2, MulOp>(N, show_details);
    CheckOpPrecisions<FP6_E3M2, AddOp>(N, show_details);
    CheckOpPrecisions<FP6_E2M3, MulOp>(N, show_details);
    CheckOpPrecisions<FP6_E2M3, AddOp>(N, show_details);
    CheckOpPrecisions<FP4_E3M0, MulOp>(N, show_details);
    CheckOpPrecisions<FP4_E3M0, AddOp>(N, show_details);
    CheckOpPrecisions<FP4_E2M1, MulOp>(N, show_details);
    CheckOpPrecisions<FP4_E2M1, AddOp>(N, show_details);

    const int DP_RUN = 20;
    CheckDotProducts<FP8_E5M2>(DP_RUN);
    CheckDotProducts<FP8_E4M3>(DP_RUN);
    CheckDotProducts<FP6_E3M2>(DP_RUN);
    CheckDotProducts<FP6_E2M3>(DP_RUN);
    CheckDotProducts<FP4_E3M0>(DP_RUN);
    CheckDotProducts<FP4_E2M1>(DP_RUN);
}

template <typename FPX, typename Op>
void TestsImpl() {
    Op op;
    FPX::PrintStats(LOG(INFO));
    FPX::ListAllPositiveValues(LOG(INFO));
    
    for (int i = 0; i < (1<<FPX::EBITS); ++i) {
        for (int j = 0; j < (1<<FPX::MBITS); ++j) {
            auto a1 = FPX(0, i, j);
            LOG(INFO) << "=== change a to " << a1 << " ===";
            for (int i2 = 0; i2 < (1<<FPX::EBITS); ++i2) {
                for (int j2 = 0; j2 < (1<<FPX::MBITS); ++j2) {
                    auto b1 = FPX(0, i2, j2);
                    auto c = op(a1, b1);
                    const auto exact_c = op(a1.ToFloat32(), b1.ToFloat32());
                    const auto diff = c.ToFloat32() - exact_c;
                    const auto diff_pct = diff ? diff * 100 / exact_c : 0;
                    LOG(INFO) << a1 << op.DESC << b1 << "=" << c << " exact=" << exact_c << "(" << c.ToFloat32() - exact_c << ", " << diff_pct << "%)";
                }
            }
        }
    }

    for (double x = -10; x < std::min(FPX::NORMAL_MAX.ToFloat32() + 5, 1000.0f); x += 0.05) {
        LOG(INFO) << x << " : " << FPX::FromFloat32(x);
    }
}

int Tests(int argc, char* argv[]) {
    if (argc <= 0) {
        TestsImpl<FP6_E3M2, MulOp>();
    } else if (strcmp(argv[0], "e3m2") == 0) {
        TestsImpl<FP6_E3M2, MulOp>();
    } else if (strcmp(argv[0], "e2m3") == 0) {
        TestsImpl<FP6_E2M3, MulOp>();
    } else if (strcmp(argv[0], "e3m0") == 0) {
        TestsImpl<FP4_E3M0, MulOp>();
    } else if (strcmp(argv[0], "e2m1") == 0) {
        TestsImpl<FP4_E2M1, MulOp>();
    } else if (strcmp(argv[0], "e5m2") == 0) {
        TestsImpl<FP8_E5M2, MulOp>();
    } else if (strcmp(argv[0], "e4m3") == 0) {
        TestsImpl<FP8_E4M3, MulOp>();
    } else {
        LOG(ERROR) << "Unknown arg=" << argv[0];
    }
    return 0;
}

