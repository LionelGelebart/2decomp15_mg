# 2decomp15_mg for 2decomp developers : 

multigrid evolution from the 1.5 version of 2decomp

# Install
cp src/Makefile.inc.x86 src/Makefile.inc

Adjust FFT=fftw3_f03

Adjust FFT_inc and FFT_lib

# Test example : fft_test_r2c_mg

run two fft/ifft with diffrent grid sizes

Note : tested with FFT=generic

# commit  e1d7a8ac3 : slight modification

	the pointer association also acts on public variable of 2decomp
	so that every pointers are on the correct grid
