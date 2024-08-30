def print_float_encoding(E:int, M:int):
    EBIAS = 2**(E-1) - 1

    def get_float_value(sign, exponent, mantissa):
        if E >= 5:
            if exponent == 2**E-1:
                if mantissa == 0:
                    return float('-inf') if sign else float('inf')
                return float('nan')
        elif E == 4:
            if exponent == 2**E-1 and mantissa == 2**M-1:
                return float('nan')
            
        sign_value = -1 if sign else 1
        return sign_value * 2**(exponent - EBIAS) * (1 + mantissa/(2**M)) if exponent != 0 else sign_value * 2**(1-EBIAS) * (0 + mantissa/(2**M))

    index = 0
    for sign in range(2):
        for exponent in range(2**E):
            for mantissa in range(2**M):
                val = get_float_value(sign, exponent, mantissa)
                abserr = 2**(exponent-EBIAS-M-2)
                print(f"{index}) {sign}_{exponent}_{mantissa} -> {val}  //Err: {abserr} / {abserr/(val+1e-9)}")
                index += 1

print_float_encoding(3,2)
