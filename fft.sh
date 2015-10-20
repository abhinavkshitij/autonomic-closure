echo running 1D FFT code
gfortran -o FFT tryfftw1.f90 -L /opt/fftw-3.3.4/lib/ -lfftw3 -lm
./FFT
#cat fft_mode.csv
#python fft_plot.py
