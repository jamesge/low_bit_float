lbfc: low_bit_float_calc.cc
	clang++ -O2 low_bit_float_calc.cc -o lbfc

all: lbfc

clean:
	rm -f ./lbfc
