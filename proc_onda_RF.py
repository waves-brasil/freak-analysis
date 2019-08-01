from pylab import *
import time
import os
import datetime

########################################
# DECLARACAO DE VARIAVEIS E CONSTANTES #
########################################

# Numeros das colunas com os dados nos arquivos do arquivo DF037:
data=0; Hora=1; elev=3;

fileNames = []
mare=[]
# Cria/Limpa as variaveis:   #cr=[];ca=[];h=[]

 ######################################################################
 #                                                                    #
 #     Funcao para calculo dos parametros no dominio do tempo:        #
 #                                                                    #
 ######################################################################

def espectro(x,dt,janseg):
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
    # Janela de Hanning com tamanho dos segmentos, ou seja, "n".
    h=hanning(n)

    # Vetor de frequencias:
    f = arange(1./(n*dt), (n/2.)/(n*dt), 1./(n*dt)) # Frequencia fundamental (menor freq)=1/ periodo de medicao, nesse caso = "n"
    L = len(x1)-n+1
    y1=zeros([n,L],)
    for i in range(0,int(L)):
       #k=n*(i-1.)+1  # L segmentos com "n" pontos, completa as "i"
        k=i
        y1[:,i]=x1[k:k+n]*h  # de y com os L segmentos. e multiplicando os segmentos pela janela "h".
    
    # Fazendo a fft de y:
    z1=fft(y1.conj().transpose())
    z1=z1.conj().transpose()
    z1=z1[(1+1):(n/2+1),:]
    zx=z1*conj(z1)*(2/(n*.375)) #(2/(n*.375))  #(n*.375)
    S=real(zx.mean(axis=1))
    return S, f


def paramts(x,segm,nfft,data):
    # Dividindo a serie temporal em segmentos de 1024 segundos (no caso de 2Hz, 2048pts)
    Data_Seg=[]
    
    cont=1
    for k in range(segm,len(x),segm):
        eta=x[k-segm:k]
        mare.append(mean(eta))  
        # Retira a media
        eta=eta-mean(eta)
        
        Data_Seg.append(data[k-segm])
        
	# Acha os cruzamentos de zero descendente, quando o sinal da serie muda de pos para neg 
        izd=find(diff(sign(eta))==-2)
        # Acha os cruzamentos de zero ascendente, quando o sinal da serie muda de neg para pos
        iza=find(diff(sign(eta))== 2)
        
	# Limpa as cristas e cavas obtidos no loop anterio
        ca=[];cr=[];h_a=[];h_d=[]
        ca_a=[];ca_d=[];
        cr_a=[];cr_d=[];
        t_a=[];sh=[];hs=[];
	
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
            icr=find(diff(sign(diff(eta_ind)))==-2)+1
            cr.append(eta_ind[icr])
            h_a.append(cr_a[i]-ca_a[i])
            t_a.append(((iza[i+1]+1)-iza[i])/2.)

        # ALTURA SIGNIFICATIVA --------
        # Organiza a serie de H em ordem da menor para a maior:
        sh=sort(h_a)
        vh.append(h_a)
        vt.append(t_a)
        CR.append(cr)

	#  vai ser a media do ultimo terco da serie de altura em ordem cresente:
        hs=mean(sh[((2*len(sh)/3)):len(sh)])

        # ALTURA MAXIMA ---------------
        hmax = max(sh)
        Fw=hmax/hs

##        print fw.argmax()
        crista=[]
##        print cont,Fw
        #
        # NUMERO DO SEGMENTO QUE CONTEM A FREAKWAVE
        # 76
        if cont == 76:
#-grafico freakwave            figure(3),plot(eta,'-o'),title('elev do reg selec'),  grid(True)
            for i in flatten(cr): crista.append(i)
#-grafico freakwave            figure(4),plot(crista,'-o'),title('cristas do registro selec'),  grid(True)
            print hmax/hs, hs
        cont=cont+1
##        print cr
        for i in flatten(cr): crista.append(i)
        fw.append(max(crista))

	# As 40 maiores alturas (pedido pela US-SAE)
        hmax40=sh[(len(sh)-40):len(sh)]

        # Razao entre Hmax e Hs
        hmax_hs=hmax/hs
        
        # PERIODOS
        Tind=diff(izd)/2. # 2Hz

        # Acrescenta o valor de Hs e Hmax calculado para a ultima serie de 2048 pts:
        Hs.append(hs)
        
        Hmax.append(hmax)

        S,f=espectro(eta,.5,nfft)

        # Periodo de pico
        ind=S.argmax()
        tp.append(1/f[ind])

        # Altura Significativa Hs, Hmax ,  CR,f,S, vh,vt,hm0,tp, 
        df=f[1]-f[0]
        hm0.append(4*sqrt(S.sum()*df))
    
    return Hs, Data_Seg


#############################################################
#                   LEITURA DOS ARQUIVOS                    #
#############################################################
Hs=[];  Hmax=[]; Hm0=[]; hs=[];hm0=[];h=[];
HM0=[];HMAX=[];vh=[]; HS=[]
Tp=[];tp=[];vt=[]; t=[];T=[];
CRI=[];Cr=[] ; Crist=[];
#Data_Seg=[];
data2=[];DaTa=[];DaTa2=[]
CR=[];fw=[]


# MAIN LOOP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Lista os arquivos com extencao definida abaixo:
for f in os.listdir(os.getcwd()):
    # Define qual extensao dos arquivos de dados:
    if f.endswith('.DF037'):
        fileNames.append(f)
# Organiza os nomes dos arquivos pelo nome (se o nome contiver a data...)
fileNames=sorted(fileNames)

#~~~~~~~~~ Abre os arquivos ~~~~~~~~~~~~~~~~~~~#

for fileName in fileNames:
    # Mostra na tela qual arquivo esta sendo processado:
    print fileName
    # Abre o arquivo
    f = open(fileName)
    # Cria/Limpa as variaveis:
    DATA=[]; ELEV=[];HORA=[];
    data1=[]

    #~~~~~~ Le as linhas do arquivo aberto e separa as variaveis:
    for line in f.readlines():
        # Verifica se o primeiro caracter eh numero
        # se for, eh onde termina o cabecalho e comecam os dados.
        if line[0].isdigit():
            DATA.append(line.split()[data])
            HORA.append(line.split()[Hora])
            ELEV.append(float(line.split()[elev]))

    # DATA:-------------------------------------------------------
    for i in range(0, len(DATA)):
        ano=DATA[i][0:4]
        mes=DATA[i][5:7]
        dia=DATA[i][8:10]
        hora=HORA[i][0:2]
        minu=HORA[i][3:5]
        segu=HORA[i][6:8]
        # passar para inteiro!
        data1.append(int((24*3600*10.)*date2num(datetime.datetime(int(ano),int(mes),int(dia),int(hora),int(minu),int(segu)))))

    # Series temporal de elevacao
    el = array(ELEV)

    # FECHA O ARQUIVO DE DADOS
    f.close()

    iel=[]
##    el=el-mean(el)   
    iel=find(el>0)
    if len(iel)>0:
        el=el[iel]
        for k in iel: data2.append(data1[k])
        
##    print data2[-1]
##--    el=el-mean(el)
    
    Hs, Data_Seg = paramts(el,2048,256,data2)
    
##    hs.append(HS) HS, Hmax ,  CR, f,S, vh,vt,hm0, tp,
    DaTa.append(Data_Seg)


##### MARE TSUNAMI
####mare_filt=[]
####for i in range(1,len(mare)-1):
####    mare_filt.append((mare[i-1]+2*mare[i]+mare[i+1])/4.)
####mare_filt.insert(0, mare[0])
####mare_filt.append(mare[-1])
####alta_freq=array(mare)-array(mare_filt)
####figure(55),plot(alta_freq)

##    HS.append(Hs)

##    Hmax.append(hmax)
##    Cr.append(CR)
##    Hm0.append(hm0)
##    Tp.append(tp)
    
##    figure(1),  plot(f,S)
##    grid(True)


for i in flatten(hm0): HM0.append(i),
for i in flatten(Hmax): HMAX.append(i),
for i in flatten(CR): CRI.append(i),
for i in flatten(tp): T.append(i),
for i in flatten(Hs): HS.append(i),
for i in flatten(DaTa): DaTa2.append(i),
##for i in flatten(vh): h.append(i),
##for i in flatten(vt): t.append(i),

Hs=array(Hs)
T=array(T)
HM0=array(HM0)
HMAX=array(HMAX)


escrita = open('DF037.txt', 'w')
escrita.write("###################" + "\n")
escrita.write("" + "\n")
for k in range(0,len(Hs)):
    escrita.write(str(DaTa[0][k]) + " " + str(Hs[k]) + " " + str(HMAX[k]) + " " + str(HM0[k])  + " " + str(T[k]) + " " + str(T[k]) +'\n')
escrita.close()

escrita = open('DF037_CR.txt', 'w')
escrita.write("###################" + "\n")
for k in range(0,len(CRI)):
    escrita.write(str(CRI[k])+'\n')
escrita.close()

#-figure(1),plot(CRI,'-o'),grid(True),title('Cristas')

#-figure(2),plot(HMAX/Hs,'-or'),grid(True),title('Hmax/Hs')
##figure(2),plot(Hs,'-ob')#,grid(True),title('fw')
#-figure(2),plot(fw,'--og'),grid(True),title('fw')

##h=array(Hs)
##hm=array(Hmax)
#hi=array(Hind)




