from my_test import *
import scipy.signal as signal
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#
def single_window_test11(N=101, W=7, snr=10 , L=1):
    l = 1
    l2 =2
    value = np.random.randn(N)+1j*np.random.randn(N)+0
 #   value = np.random.randn(N)
    g_value = np.zeros(N)
    g_value[0:W-1] = 1
    g = window(g_value, N)
    x = Signal(value, N, snr)
    # print(x.get_xl(1))
    # print("G_l", g.compute_Gl(l,L))
    # print(np.shape(g.compute_Gl(l, L)))
    # print(np.shape(x.compute_yc(g, L)))
    re = compute_rdistance(x.initial, x.retrieval_x(g,l, l2, L))
    print(compute_rdistance(x.initial, x.retrieval_x(g,l, l2, L)))
    print(compute_rdistance(x.initial, x.using_xl()))
    print("same: ", compute_rdistance(x.initial, x.initial))
    return re

def single_window_test(signal, g, L=1):
    l = 1
    l2 =2
    return compute_rdistance(signal.initial, signal.retrieval_x(g, l, l2, L))


def choose_signal(name='Gauss complex', N=101, snr=30, M=0):
    if name == 'Gauss complex':
        value = np.random.randn(N)+1j*np.random.randn(N)+M
        x = Signal(value, N, snr)
    elif name == 'chirp':
        t = np.arange(0, N)
        value = signal.chirp(t, 0, N, f1=1)+M
        x = Signal(value, N, snr)
    elif name == 'temper':
        pass
    return x



def choose_windows(name='Hamming', W=7,N=101):
    # Rect/Hanning/Hamming
    if name == 'Hamming':
        window1 = np.zeros(N)
        for n in range(W):
            window1[n] = 0.54 - 0.46 * np.cos(2 * np.pi * n / (W - 1))
    elif name == 'Hanning':
        window1 = np.zeros(N)
        for n in range(W):
            window1[n] = 0.5 - 0.5 * np.cos(2 * np.pi * n / (W - 1))
    elif name == 'Rect':
        window1 = np.zeros(N)
        window1[0:W-1]=1
    elif name == 'Tri':
        window1 = np.zeros(N)
        for n in range(int(np.floor(W/2))):
            window1[n] = n/W
        for n in range(int(np.floor(W/2)),W):
            window1[n]=1- n/W
    return window(window1,N)



def plot_re_snr(name='Gauss complex'):
    g_list = ['Hamming', 'Tri','Rect']
    relative_snr_all = []
    i_snr = np.arange(10,100,10)
    for g in g_list:
        relative_snr = []
        for i in i_snr:
            x = choose_signal(name, N=101, snr=i, M=10)
            re = single_window_test(x, choose_windows(g), L=1)
            relative_snr.append(re)
        relative_snr_all.append(relative_snr)
    plt.plot(i_snr, relative_snr_all[0], 'bv--',mew='4',label='Hamming window')
    plt.plot(i_snr, relative_snr_all[1], 'ro--',mew='4',label='Triangle window')
    plt.plot(i_snr, relative_snr_all[2], 'g^--',mew='4',label='Rectangle window')
    plt.title(name + ' signal')
    plt.legend(fontsize='x-large')
    plt.ylabel("Relative error")
    plt.xlabel("SNR")
    plt.savefig('image/re_snr_window_{}.pdf'.format(name))
    # plt.show()
    plt.close()

def plot_re_snr_M(name='Guass complex'):
    g = 'Hamming'
    KM= np.arange(0,80,20)
    relative_snr_all = []
    i_snr = np.arange(10,100,10)
    for K in KM:
        relative_snr = []
        for i in i_snr:
            x = choose_signal(name, N=101, snr=i,M=K)
            re = single_window_test(x, choose_windows(g), L=1)
            relative_snr.append(re)
        relative_snr_all.append(relative_snr)
    plt.plot(i_snr, relative_snr_all[0], 'bs--', mew='3',label='$\mu$=0, {} signal'.format(name))
    # plt.legend()
    plt.plot(i_snr, relative_snr_all[1], 'ro--', mew='3', label='$\mu$=20, {} signal'.format(name))
    # plt.legend()
    plt.plot(i_snr, relative_snr_all[2], 'g^--', mew='3',label='$\mu$=40, {} signal'.format(name))

    plt.plot(i_snr, relative_snr_all[3], 'yv--', mew='3', label='$\mu$=60, {} signal'.format(name))
    plt.title(name + ' signal')
    plt.legend(fontsize='x-large')
    plt.ylabel("Relative error")
    plt.xlabel("SNR")
    plt.savefig('image/re_snr_M_{}.pdf'.format(name))

    # plt.show()
    plt.close()


def plot_re_W(name='Gauss complex'):
    g_list = ['Hamming', 'Tri','Rect']
    relative_snr_all = []
    W_list = np.arange(5,100,3)
    for g in g_list:
        relative_snr = []
        for W in W_list:
            x = choose_signal(name, N=101, snr=100, M=20)
            re = single_window_test(x, choose_windows(g, W=W), L=1)
            relative_snr.append(re)
        relative_snr_all.append(relative_snr)
    plt.plot(W_list, relative_snr_all[0], 'bv--',mew='3',label='Hamming window')
    plt.plot(W_list, relative_snr_all[1], 'ro--',mew='3',label='Triangle window')
    plt.plot(W_list, relative_snr_all[2], 'g^--',mew='3',label='Rectangle window')
    plt.title(name + ' signal')
    plt.legend(fontsize='x-large')
    plt.ylabel("Relative error")
    plt.xlabel("W(length of window)")
    plt.axis([0, 80, 0, 0.2])
    plt.savefig('image/re_snr_window_W_{}.pdf'.format(name))

    # plt.show()
    plt.close()



def plot_re_N(name='chirp'):
    g_list = ['Hamming', 'Tri','Rect']
    relative_snr_all = []
    N_list = np.arange(55,101,2)
    for g in g_list:
        relative_snr = []
        for N in N_list:
            x = choose_signal(name, N=N, snr=80, M=10)
            re = single_window_test(x, choose_windows(g, W=7, N=N), L=1)
            relative_snr.append(re)
        relative_snr_all.append(relative_snr)
    plt.plot(N_list, relative_snr_all[0], 'bv--',mew='3',label='Hamming window')
    plt.plot(N_list, relative_snr_all[1], 'ro--',mew='3',label='Triangle window')
    plt.plot(N_list, relative_snr_all[2], 'g^--',mew='3',label='Rectangle window')
    plt.title(name + ' signal')
    plt.legend(fontsize='x-large')
    plt.ylabel("Relative error")
    plt.xlabel("W(length of window)")
    plt.axis([0, 80, 0, 0.2])
    plt.savefig('image/re_snr_window_n_{}.pdf'.format(name))

    # plt.show()


def Test():
    N=101
    t = np.arange(0,N)
    aa=signal.chirp(t,0,101,f1=1)
    print(len(aa))


if __name__ == "__main__":
    # plot_re_snr(name='Gauss complex')
    plot_re_snr_M(name='Gauss complex')
    # plot_re_W(name='Gauss complex')
    # plot_re_snr(name='chirp')
    # plot_re_N('chirp')
    # plot_re_snr_M(name='chirp')
    # plot_re_W(name='chirp')
    # single_window_test11()
    # Test()



