import numpy as np
import math
from utils import *


class Signal:
    def __init__(self,initial, N, snr):
        self.initial = initial
        self.N = N
        self.snr = snr

    def stft(self, g , L):
        vec = np.array(list(range(0, self.N)))
        m = np.array(list(range(0, self.N, L)))
        yy = np.zeros([len(m), self.N], dtype='complex')
        for ii in range(len(m)):
            index = np.mod(m[ii]-vec, self.N)
            xg = np.multiply(self.initial, g.initial[index])
            # print(np.fft.fft(xg))
            yy[ii, :] = np.fft.fft(xg)
        return yy

    def noise_measurement(self, g, L):
        x = self.stft(g, L)
        return noise_mea(self.snr, x)

    def compute_yc(self,g, L):
        tt = np.fft.ifft(np.abs(self.noise_measurement(g, L))**2)
        return tt

    def compute_ycl(self, g, l, L):
        return self.compute_yc(g, L)[:, l]



    def get_xl(self,l):
        xl = np.zeros(self.N, dtype='complex')
        for i in range(self.N):
            index = np.mod(i+l, self.N)
            xl[i] = np.conjugate(self.initial[i])*(self.initial[index])
        return xl

    def compute_single_xl(self, g, l, L):
        gl = g.compute_Gl(l, L)
        xl = np.dot(np.linalg.inv(gl), self.compute_ycl(g, l, L))
        return xl

    def retrieval_x(self, g, l1,l2, L):
        x_hat_phase = np.zeros(self.N, dtype='complex')
        x_mod = np.zeros(self.N, dtype='complex')
        xl1 = self.compute_single_xl(g, l1, L)
        xl2 = self.compute_single_xl(g, l2, L)
        x_hat_phase[0] = 1
        for i in range(1, self.N):
            x_hat_phase[i] = x_hat_phase[i-1]*phase(xl1)[i-1]
            x_mod[i] = xl1[i-1]*xl1[i]/xl2[i-1]
        x_mod[0] = xl1[self.N-1]*xl1[0]/xl2[self.N-1]
        return np.sqrt(x_mod) * x_hat_phase

    def using_xl(self):
        x_hat_phase = np.zeros(self.N, dtype='complex')
        x_mod = np.zeros(self.N, dtype='complex')
        xl1 = self.get_xl(1)
        xl2 = self.get_xl(2)
        x_hat_phase[0] = 1
        for i in range(1, self.N):
            x_hat_phase[i] = x_hat_phase[i-1]*phase(xl1)[i-1]
            x_mod[i] = xl1[i-1]*xl1[i]/xl2[i-1]
        x_mod[0] = xl1[self.N-1]*xl1[0]/xl2[self.N-1]
        # print(x_hat_phase)
        return np.sqrt(x_mod) * x_hat_phase

    # def compute_xl(self, L, **kwargs):
    #     xl =
    #
    #     for j1 in range(L):
    #         for j in range(L):
    #             for m in range(np.ceil(self.N/L)):
    #                 for m1 in range(np.ceil(self.N/L)):
    #                     xl(n) =np.exp(-2*math.pi*1j*(m*(m1*L-n)/N-j*n/L))
    #


class window:
    def __init__(self,initial,length):
        self.initial = initial
        self.length = length

    def compute_Gl(self,l,L):
        g_l = np.zeros([int(np.ceil(self.length/L)), self.length], dtype='complex')
        gg = np.zeros(self.length, dtype='complex')
        for i in range(self.length):
            index = np.mod(i-l, self.length)
            gg[i] = self.initial[np.mod(i, self.length)]*np.conjugate(self.initial[index])
        for m in range(int(np.ceil(self.length/L))):
            for n in range(self.length):
                g_l[m,n] = gg[int(np.mod(m*L-n, self.length))]
        return g_l



def test_single():
    N = 100
    W = 6
    snr = 5
    value = np.random.randn(N)+1j*np.random.randn(N)
    g_value = np.zeros(N)
    g_value[0:W-1] = 1
    g = window(g_value, N)
    x = Signal(value, N, snr)
    y = x.stft(g, 2)
    Yc = np.abs(y)




def Test2():
    N = 101
    W = 7
    snr = 30
    L = 1
    l = 1
    l2 =2
    value = np.random.randn(N)+1j*np.random.randn(N)+1000000
 #   value = np.random.randn(N)
    g_value = np.zeros(N)
    g_value[0:W-1] = 1
    g = window(g_value, N)
    x = Signal(value, N, snr)
    # print(x.get_xl(1))
    # print("G_l", g.compute_Gl(l,L))
    # print(np.shape(g.compute_Gl(l, L)))
    # print(np.shape(x.compute_yc(g, L)))
    print(compute_rdistance(x.initial, x.retrieval_x(g,l, l2, L)))
    print(compute_rdistance(x.initial, x.using_xl()))
    print("same: ", compute_rdistance(x.initial, x.initial))


def test_kwargs(first, *args, **kwargs):
    print('Required argument: ', first)
    for v in args:
        print('Optional argument (*args): ', v)
    for k, v in kwargs.items():
        print('Optional argument %s (*kwargs): %s' % (k, v))
#
#
def compute_y_orgianl():
    N = 100
    W = 6
    snr = 10
    L = 1
    l = 1
    l2 =2
    value = np.random.randn(N)+1j*np.random.randn(N)
    #value = np.random.randn(N)
    g_value = np.zeros(N)
    g_value[0:W-1] = 1
    g = window(g_value, N)
    x = Signal(value, N, snr)
    x.compute_yc()


def Test3():
    a = 1+1j
    b = 1+1j
    print(compute_rdistance(a,b))

if __name__ == "__main__":
     # Test1()
  #   keylist = {'121':1,'12':2}
  #     test_single()
  #  test_kwargs(1,2,4,k1=1,k2=3)
  #   Test2()
    Test2()
#