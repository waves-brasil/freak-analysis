# -*- coding: cp1252 -*-
import sys
sys.path.append("C:\\Python25\\Lib\\site-packages,C:\\Python25\\Lib")
##from pylab import *
##from matplotlib import plot
from scipy import  fix, mean, zeros, arange, fft, conj, transpose, arctan2, cos, where #, find
from scipy import  real, ifft, array, cosh, tanh, pi, sum, sqrt, ones, hanning, sinh, nan, isnan
#from scipy.signal import cspline1d, cspline1d_eval
import scipy.interpolate as sp
import scipy.interpolate
##from scipy.misc.common import factorial
##from scipy.interpolate import spline

#from pylab import  fix, mean, zeros, arange, fft, conj, transpose, find, arctan2, cos

#from pylab import  real, ifft, array, cosh, tanh, pi, sum, sqrt, ones, hanning, sinh
########################################
# DECLARACAO DE VARIAVEIS E CONSTANTES #
########################################

# Gravidade
g=9.81

# Lamina dagua = batim lido no arquivo fsiP.txt

########################################
def crospec(x,x2,dt,janseg):
    # x = sinal
    # dt = intervalo de amostragem em segundos
    # janseg = janela em segundos
    # n= Numero de pontos da fft, "n", sera igual ao numero de frequencias
    # resultantes sendo que a ultima metade das frequencias e o "reflexo" da
    # primeira metade. Assim, usa-se apenas n/2 frequencias.
    # n = 64; n = janela em segundos * numero de pontos entre cada segundo =
    # frequencia de amostragem (1/dt):
    # n = 128;
    n = janseg * 1./dt

    n21=fix(n/2)+1
    lx=len(x)
    x1 = x-mean(x) # retira-se a media do sinal.
    x2 = x2-mean(x2)
    # Janela de Hanning com tamanho dos segmentos, ou seja, "n".
    h=hanning(n)

    # Vetor de frequencias:
    f = arange(1./(n*dt), (n/2.)/(n*dt), 1./(n*dt)) # Frequencia fundamental (menor freq)=1/ periodo de medicao, nesse caso = "n"
    L = len(x1)-n+1
    y1=zeros([n,L],)
    y2=zeros([n,L],)
    Z=[]
    for i in range(0,int(L)):
       #k=n*(i-1.)+1  # L segmentos com "n" pontos, completa as "i"
        k=i
        y1[:,i]=x1[k:k+n]*h  # de y com os L segmentos. e multiplicando os segmentos pela janela "h".
        y2[:,i]=x2[k:k+n]*h
        z=abs(y1[:,i]*conj(y2[:,i]))
        Z.append(z.sum())
    Z=array(Z)
    mZ=mean(Z)
    #iZ=find(Z>0.8*mZ)
    iZZ=where(Z>0.8*mZ)
    iZ=iZZ[0]

    # Fazendo a fft de y:
    z1=fft(y1[:,iZ].conj().transpose())
    z2=fft(y2[:,iZ].conj().transpose())
    z1=z1.conj().transpose()
    z2=z2.conj().transpose()
    z1=z1[(1+1):(n/2+1),:]
    z2=z2[(1+1):(n/2+1),:]
    zx=z1*conj(z2)*(2/(n*.375)) #(2/(n*.375))  #(n*.375)
    S=real(zx.mean(axis=1))
    return S, f



def espectro(x,dt,nfft):
    # x = sinal
    # dt = intervalo de amostragem em segundos
    # janseg = janela em segundos
    # n= Numero de pontos da fft, "n", sera igual ao numero de frequencias
    # resultantes sendo que a ultima metade das frequencias e o "reflexo" da
    # primeira metade. Assim, usa-se apenas n/2 frequencias.
    # n = 64; n = janela em segundos * numero de pontos entre cada segundo =
    # frequencia de amostragem (1/dt):
    # n = 128;
    n = nfft

    n21=fix(n/2)+1
    lx=len(x)
    x1 = x-mean(x) # retira-se a media do sinal.
    # Janela de Hanning com tamanho dos segmentos, ou seja, "n".
    h=hanning(n)
    Novlp=nfft/2
    continuo=0
    Overlap=1

    # Vetor de frequencias:
    f = arange(1./(n*dt)+ 1./(n*dt), (n/2.)/(n*dt) + 1./(n*dt), 1./(n*dt)) # Frequencia fundamental (menor freq)=1/ periodo de medicao, nesse caso = "n"
    L = len(x1)-n+1
    if continuo==1:
        y1=zeros([n,L],)
        for i in range(0,int(L)):
            #k=n*(i-1.)+1  # L segmentos com "n" pontos, completa as "i"
            k=i
            y1[:,i]=x1[k:k+n]*h  # de y com os L segmentos. e multiplicando os segmentos pela janela "h".
    elif Overlap==1:
        L=len(x1)
        i=0
        J=0
        M=nfft-Novlp
        y1=zeros([n,(L/M)-1],)
        while J+nfft<=L:
            y1[:,i]=x1[J:(J+nfft)]*h;
            J=J+M
            i=i+1
        
    # Fazendo a fft de y:
    z1=fft(y1.conj().transpose())
    z1=z1.conj().transpose()
    z1=z1[(1+1):(n/2+1),:]
    zx=z1*conj(z1)*(2/(n*.375)) #(2/(n*.375))  #(n*.375)
    S=real(zx.mean(axis=1))
    #contour(zx)
    #plot(f,Sp)
    #Tp = 1/f[ind]
    return S, f


def TranPres(p,u,v,w,nfft,dt,hlocal,cota):
    n21 = nfft/2.
    # Criando vetor de frequencia:
    f = arange(0.,1./(2.*dt),1./(nfft * dt))
    # Funcao de transferencia:
    Kp,Kvel=FuncTran(f,hlocal,cota)
    # Vetor linha
    Kp=Kp.flatten()
    Kvel=Kvel.flatten()

    # Retirando a media:
    P=p[0:nfft]-mean(p[0:nfft])
    U=u[0:nfft]-mean(u[0:nfft])
    V=v[0:nfft]-mean(v[0:nfft])
    W=w[0:nfft]-mean(w[0:nfft])
    # FFT da serie bruta:
    fp=fft(P)
    fu=fft(U)
    fv=fft(V)
    fw=fft(W)
    # Passando para vetor linha
    fp=fp.flatten()
    fu=fu.flatten()
    fv=fv.flatten()
    fw=fw.flatten()
    # Aplica funcao de transferencia na primeira metade do espectro:
    fp[0:n21] = fp[0:n21] * Kp
    fu[0:n21] = fu[0:n21] * Kvel
    fv[0:n21] = fv[0:n21] * Kvel
    fw[0:n21] = fw[0:n21] * Kvel

    # A segunda metade do espectro eh o complexo conjugado da primeira:
    # parte real eh invertida a complexa eh invertida e o sinal eh negativo = complex conj
    fp[:n21:-1]=conj(fp[1:n21])
    fu[:n21:-1]=conj(fu[1:n21])
    fv[:n21:-1]=conj(fv[1:n21])
    fw[:n21:-1]=conj(fw[1:n21])

    Sp =fp*conj(fp)*(2./nfft)

##    Su =fu*conj(fu)*(2/nfft)
##    Sv =fv*conj(fv)*(2/nfft)
    Sp =Sp[0:n21]
    #hf=find(f>1./4)
    hff=where(f>1./4)
    hf=hff[0]
    #lf=find(f<1./20)
    lff=where(f<1./20)
    lf=[0]

    fp[hf]=0; fp[lf]=0; fu[hf]=0; fu[lf]=0; fv[hf]=0; fv[lf]=0;fw[hf]=0; fw[lf]=0;
    ind=Sp.argmax()

    Tp = 1/f[ind]
    # Transf inversa de Fourier
    Pp = real(ifft(fp))
    Uu = real(ifft(fu))
    Vv = real(ifft(fv))
    Ww = real(ifft(fw))

    return Pp,Uu,Vv,Ww,Tp

#-----------------------------------------------------------------------------
# Aplica funcao de transferencia Kp a um vetor dado (pressao, espectro, etc)
#-----------------------------------------------------------------------------
def FuncTran(f, hlocal,cota):
    h=hlocal
    prof=cota
    n21 = len(f)

    # Calculo do Numero de Onda %
    # Retorna o comprimento de onda atraves da equacao de dispersao da
    # teoria de onda linear pela formula de HUNT (1979), conforme
    # descrito no livSp=FuncTran(Ppp, freqs, prof)ro "Water Wave Mechanics for Engineers and
    # Scientists", Dean and Dlrymple (1984), pag. 72
    y = pow(2.*pi*f,2)* (prof+h) / 9.81  # frequencia angular
    d = [0.6666666667, 0.3555555555, 0.1608465608, 0.0632098765, 0.0217540484, 0.0065407983]
    k=zeros([len(y),1],)
    for i in range(0,len(y)):
        aux = pow(y[i],2) + y[i]/(1 + sum(d[0]*y[i]+d[1]*pow(y[i],2)+d[2]*pow(y[i],3)+d[3]*pow(y[i],4)+d[4]* pow(y[i],5)+d[5]*pow(y[i],6)))
        k[i] = sqrt(aux)/(prof+h)

    # CALCULO DAS FUNCOES DE TRANSFERENCIAS
    # Pressao - Funcao de transferencia invertida !!rp=(rho*g)*cosh(k*(h_bat+z))./cosh(k*h_bat);
    Kvel=ones(len(f))
    Kvel=Kvel.flatten()
    for ff in range(0,len(f)):
        Kvel[ff] = (2.*pi*f[ff])*sinh(k[ff]*h)/cosh(k[ff]*(prof+h))
    Kp = cosh(k*(prof+h))/cosh(k*h)

    # Aplicando filtro passa-baixa na funcao de transferencia
    Tlpass = 4.0 # Periodo para zerar a funcao de transferencia (segundos)
    #ki = find(f >= 1./Tlpass)
    kii = where(f >= 1./Tlpass)
    ki=kii[0]

    klpass = ki[0]
    kT = (n21-klpass)*4.
    lowpass = ones([n21,1],)
    lowpass[klpass:n21] = zeros([n21 - klpass,1],) # Zerando a funcao de transferencia a partir de Tlpass
    Tlpass = 6.0  # Periodo para inicio da atenuacao (segundos)
    #iK = find(f >= 1./Tlpass)
    ikk = where(f >= 1./Tlpass)
    iK=ikk[0]

    klpass =iK[0]
    kT = (n21 - klpass) * 4.

    # Lowpass
    Cos=cos(2.*pi*((arange(klpass,n21))-klpass)/kT)
    COS=zeros([len(Cos),1],)
    for i in range(0,len(Cos)):
        COS[i]=pow(Cos[i],32)
    lowpass[klpass:n21] = lowpass[klpass:n21]*COS

    Kp = Kp * lowpass
    for g in range(0,len(Kvel)):
        Kvel[g] = Kvel[g]*lowpass[g]

    return Kp,Kvel

def FuncTran2(f, hlocal,cota):
    prof=cota
    # Numero de onda
    vk=wtok(f,prof)
    vk=array(vk)




    fc1=1./4. # 4 seg
    fc2=1./25.  # 25 segundos
    # Funcao de transferencia
    Kp=cosh(vk*(h+cota))/cosh(vk*h)
    #Ku=vw*sinh(vk*h)/cosh(vk*(h+cota))


    ifc11=where(f>fc1); ifc22=where(f<fc2)
    ifc1=ifc11[0]; ifc2=ifc22[0]

    Kp[ifc1]=0; Kp[ifc2]=0

    return Kp


# Funcao de Transformacao de frequencia angular para Num de onda:
def wtok(vf,h):
    vk=[]
    for f in vf:
        w=2*pi*f
        ki=w*w/g
        # Primeiro chute para k
        kf=3*ki
        while(abs(kf-ki) > kf*.0001):
            ki=kf
            kf  = w*w/(g*tanh(ki*h))
        vk.append(kf)
    return vk




#############################################################
#                   LEITURA DOS ARQUIVOS                    #
#############################################################
# Numeros das colunas com os dados nos arquivos do fsi:
#velx=0; vely=1;
velN=0; velE=1; pres=2;
Demag=0;Aproa=1;Batim=2;

#~~~~~~~~~~~~~~~~~~~~ Abre os arquivos ~~~~~~~~~~~~~~~~~~~#
fP = open('fsiP.txt')
for lineP in fP.readlines():
    try:
        demag=float(lineP.split()[Demag])
        aproa=float(lineP.split()[Aproa])
        batim=float(lineP.split()[Batim])
    except ValueError:
        print 1
        demag=0.
        aproa=0.
        batim=100.


arq = open('fsiE.txt')
PRES=[]; VN=[]; VE=[];

#~~~~ Le as linhas do arquivo aberto e separa as variaveis:
for line in arq.readlines():
    try:
        VN.append(float(line.split()[velN]))
    except ValueError:
        #print 'Nan na serie'
        VN.append(mean(VN)) #VN[-1])
        
    try:
        VE.append(float(line.split()[velE]))
    except ValueError:
        #print 'Nan na serie'
        VE.append(mean(VE)) #VE[-1])
        
    try:
        PRES.append(float(line.split()[pres]))
    except ValueError:
        #print 'Nan na serie'
        PRES.append(mean(PRES)) #nan) #PRES[-1])
            

try:

    #~~~~~~~~~~~  FIM DA LEITURA DOS ARQUIVOS   ~~~~~~~~~~~~~~#
    #############################################################

    if len(PRES) < 1200:
        while len(PRES) < 1200:
            PRES.append(float(mean(PRES)))
            VN.append(float(mean(VN)))
            VE.append(float(mean(VE)))

    PRES=array(PRES)

    PRES=(PRES-1012.5)/100

    # Profundidade do sensor: (positiva para cima do nivel mEdio do mar)
    cota=(mean(PRES))
    ######################################################################
    # Power Spectral Desity
    dt = 1 # 1 segundo
    #nfft=1024

    # Series temporais
    pp = array(PRES)
    uu = array(VN)     #uu=filtra(uu,2)#uu=uu-mean(uu)
    vv = array(VE)     #vv=filtra(vv,2)#vv=vv-mean(uu)
    ww = vv

    if isnan(pp).any():
        pp[isnan(pp)]=mean(pp[~isnan(pp)])
        
    NFFT1=1024
    # NAO MUDAR DE 1024, EH PARA CALCULAR SERIE DE ELEVACAO
    nfft=1024
    dt = 1
    elev,U,V,W,Tp=TranPres(pp,uu,vv,ww,nfft,dt,batim,cota)

    Ppp, freqs = espectro(elev,1,128)

    df=freqs[1]-freqs[0]

    # Indice da Frequencia de Pico:
    ##ind=Ppp.argmax()
    ##Tp2 = 1/freqs[ind]

    ## ALTURA SIGNIFICATIVA ##
    Hs=4*sqrt(Ppp.sum()*df)

    # Indice da Frequencia de Pico:
    Ppp2, freqs2 = espectro(elev,1,128)
    x0=freqs2
    y0=Ppp2
    dx=x0[1]-x0[0]
    dx1=dx/30.
    x1=arange(x0[0],x0[-1]+dx1,dx1)
    ##y1=spline(x0,y0,x1)

    spl = scipy.interpolate.splrep(x0, y0)
    y1 = scipy.interpolate.splev( x1, spl)

    ##f = sp.interpolate.interp1d(x0, y0,kind='cubic')
    ##y1 = f(x1)

    ##dx=x0[1]-x0[0]
    ##cj=cspline1d(y0)
    ##dx1=dx/30.
    ##x1=arange(x0[0],x0[-1]+dx1,dx1)
    ##y1=cspline1d_eval(cj,x1,dx=dx,x0=x0[0])

    indT=y1.argmax()
    Tp2 = 1/x1[indT]
    Per = Tp2  #1/freqs[indT]

    # fim do calculo do periodo.

    ind=Ppp.argmax()

    NFFT2=128
    UP=crospec(U,elev,1,NFFT2)
    VP=crospec(V,elev,1,NFFT2)

    A1=array(UP[0])
    A2=array(VP[0])

    Dir =180-arctan2(A1[ind],A2[ind])*180/pi+ demag
    if Dir<0:
        Dir=Dir+360

    ##Ppp, freqs = espectro(elev,1,NFFT2) #pp

    arq.close()
    fP.close()

    escrita = open('fsiS.txt', 'w')
    if isnan(Hs):
        escrita.write('%s%s%s%s%s%s%s\n' % ("NaN;","NaN;","NaN;","NaN;","NaN;","NaN;","NaN;"))
    else:
        escrita.write('%s%2.4f%s%3.2f%s%2.2f%s%2.2f%s\n' % ("NaN;",Per,";NaN;",Dir,";NaN;",Hs,";",Hs,";"))

            
    escrita.close()

except ValueError:

    arq.close()
    fP.close()

    escrita = open('fsiS.txt', 'w')
    escrita.write('%s%s%s%s%s%s%s\n' % ("NaN;","NaN;","NaN;","NaN;","NaN;","NaN;","NaN;"))
    escrita.close()


