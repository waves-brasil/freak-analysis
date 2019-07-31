import numpy as np
from scipy import signal
import xarray as xr
import glob
import matplotlib.pyplot as plt
import datetime
from pylab import date2num, mean, ravel, diff, sign, conj, real
from matplotlib import mlab

def read_nc(fname):
    return xr.open_dataset(fname)

def read_uvh(fname):
    with open(fname, encoding='utf-8', errors='replace') as f:
        lines = [line.strip() for line in f.readlines()]
        return lines

def int2byf(si, dt=0.20, cut_off=1/0.03):
    '''
    This function double-integrates a signal in the time domain by an 
    equivalent product in the frequency program, and cut off energy
    of periods greater than cut-off. (frequencies below 1/cut_off)

    si      = input signal;
    dt      = sampling rate [s];
    cut_off = cut_off period [s].
    '''
    
    # Signal transformation to the frequency domain.
    si = si[:] ; lsi = len(si) ; clsi = int(np.ceil( lsi/2. ))
    f = np.arange(-clsi,clsi,1.) / lsi / dt # frequency
    i = 1j # complex number i
    
    # Signal transformation to the frequency domain.
    SI = np.fft.fft(signal.detrend( si ))
    SIshift = np.fft.fftshift( SI )
    
    # integration.
    SI[:clsi] = SIshift[:clsi] / ( i * 2. * np.pi * f[:clsi] )**2 # tricky indexation (it doesn't count clsi)
    SI[clsi]    = 0.
    SI[clsi+1:] = SIshift[clsi+1:lsi] / ( i * 2. * np.pi * f[clsi+1:] )**2
    
    # Energy elimination at frequencies smaller than 1/cut-off.
    fc = 1. / cut_off; ind = np.logical_and( f >= -fc , f <= fc )
    SI[ind] = 0.
    
    # Signal transformation back to the temporal domain.
    si = np.fft.ifft(np.fft.fftshift(SI) )
    sireal = si.real
#    siimag = si.imag

    # sometimes (don't know why) the two extreme points blow up:
    if (np.abs(sireal[0] - np.mean(sireal)) > (5 * np.std(sireal))) & (np.abs(sireal[-1] - np.mean(sireal)) > (5 * np.std(sireal))):
        # Spike replacement:       
        indices = np.arange(len(sireal)); not_spk=np.delete(indices,[0,len(sireal)-1])
        sireal=np.interp(indices, not_spk, sireal[not_spk])        
    
    return sireal#,siimag
    
def calc_hs_fw(h):
    # ALTURA SIGNIFICATIVA --------
    # Organiza a serie de H em ordem da menor para a maior:
    sh=np.sort(h)

    #  vai ser a media do ultimo terco da serie de altura em ordem cresente:
    hs=mean(sh[(int(2*len(sh)/3)):len(sh)])
    # ALTURA MAXIMA ---------------
    hmax = max(sh)
    Fw=hmax/hs
    
    return hs, hmax, Fw

def numpy_flat(a):
    return list(np.array(a).flat)
    
def zero_cross(eta):
    eta=eta-mean(eta)
    # Acha os cruzamentos de zero descendente, quando o sinal da serie muda de pos para neg 
    izd1 = ravel(diff(sign(eta))==-2)
    # Acha os cruzamentos de zero ascendente, quando o sinal da serie muda de neg para pos
    iza1 = ravel(diff(sign(eta))== 2)
    izd = numpy_flat(list(np.where(izd1)))
    iza = numpy_flat(list(np.where(iza1)))
    return izd, iza

##
def ind_params(izd, iza, eta):    
    ca_d=[]
    cr_d=[]
    h_d=[]
    ca_a=[]
    cr_a=[]
    cr=[]
    h_a=[]
    t_a=[]
    # Calcula o valor maximo da elevacao entre dois  cruz de zero desc
    for i in range(0, len(izd)-1):
        ca_d.append(min(eta[izd[i]:(izd[i+1]+1)]))
        cr_d.append(max(eta[izd[i]:(izd[i+1]+1)]))
        h_d.append(cr_d[i]-ca_d[i])

        
    # Calcula o valor maximo da elevacao entre dois  cruz de zero asc
    for i in range(0, len(iza)-1):
        eta_ind=eta[iza[i]:(iza[i+1]+1)]
        ca_a.append(min(eta_ind))
        cr_a.append(max(eta_ind))
        icr = ravel(diff(sign(diff(eta_ind)))==-2)+1
        cr.append(eta_ind[icr])
        h_a.append(cr_a[i]-ca_a[i])
        t_a.append(((iza[i+1]+1)-iza[i])/2.)
        
    return ca_d, cr_d, h_d, ca_a, cr_a, cr, h_a, t_a

##    
def espectro(x, dt=0.2, window=512):
    # x = sinal
    # dt = intervalo de amostragem em segundos
    # janseg = janela em segundos
    # n= Numero de pontos da fft, "n", sera igual ao numero de frequencias
    # resultantes sendo que a ultima metade das frequencias e o "reflexo" da
    # primeira metade. Assim, usa-se apenas n/2 frequencias.
    # n = 64; n = janela em segundos * numero de pontos entre cada segundo =
    # frequencia de amostragem (1/dt):
    # n = 128;
    n = 512 #int(window * 1./dt)

    n21 = np.fix(n/2)+1
    lx=len(x)
    x1 = x-mean(x) # retira-se a media do sinal.
    # Janela de Hanning com tamanho dos segmentos, ou seja, "n".
    h = np.hanning(n)

    # Vetor de frequencias:
    f = np.arange(1./(n*dt), (n/2.)/(n*dt), 1./(n*dt)) # Frequencia fundamental (menor freq)=1/ periodo de medicao, nesse caso = "n"
    L = len(x1)-n+1
    y1= np.zeros([n,L],)
    for i in range(0,int(L)):
        #k=n*(i-1.)+1  # L segmentos com "n" pontos, completa as "i"
        k=i
        y1[:,i]=x1[k:k+n]*h  # de y com os L segmentos. e multiplicando os segmentos pela janela "h".
    
    # Fazendo a fft de y:
    z1 = np.fft.fft(y1.conj().transpose())
    z1 = z1.conj().transpose()
    z1 = z1[int(1+1):int(n/2+1),:]
    zx = z1*conj(z1)*(2/(n*.375)) #(2/(n*.375))  #(n*.375)
    S=real(zx.mean(axis=1))
    return S, f

##

# Power Spectral Density
def spec1(x,nfft,fs):
    """
    ======================================================================#
    
     Estimate the Power Spectral Density (v**2/Hz, where v is unit of x)
     using the Welch's average periodogram method
    
     Input:            x    : serie real
                       nfft : number of points for FFT
                       fs   : sampling frequency
     Output: [aa] -    col0 : frequency array
                       col1 : auto-espectrum
                       col2 : lower confidence limit
                       col3 : higher confidence limit
    
     Infos:	detrend - mean
    			window - hanning
    			noverlap (welch) - 50%
    ======================================================================#
    """
    # obtaining spectrum
    sp = mlab.psd(x,NFFT=nfft,Fs=fs,detrend=mlab.detrend_mean,scale_by_freq=True,window=mlab.window_hanning,noverlap=nfft/2)
    
    f, sp = sp[1][1:],sp[0][1:]
    #graus de liberdade
    gl = len(x) / nfft * 2
    #intervalo de confianca 95%
    ici = sp * gl / 26.12
    ics = sp * gl / 5.63
    aa = np.array([f,sp,ici,ics]).T

    return aa