echo running 1D FFT code
gfortran -o FFT tryfftw.f90 -L /opt/fftw-3.3.4/lib/ -lfftw3 -lm
./FFT
