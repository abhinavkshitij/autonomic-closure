import numpy as np
import matplotlib.pyplot as plt
import pylab

filenames = ('fft_F.csv','fft_mode.csv','fft_filter.csv','fft_B.csv')      
filenames = filenames[0:]
for f in filenames:
    print f

    data = np.loadtxt(fname=f, delimiter=',')

    fig = plt.figure(figsize=(5.0, 5.0))

    if (f == 'fft_mode.csv') or (f == 'fft_filter.csv'):
        plt.plot(data[:,1],'ok')
    else:
        plt.plot(data[:,1],'-ok')

    plt.grid(True)
    plt.show(fig)
    
