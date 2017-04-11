from numpy import *
import numpy as np
import math
import NLFG
import matplotlib.pyplot as plt

def loadDataset(filename):
    print("Loading data...")
    dataMat=[]
    nameMat=[]
    fr=open(filename,'r',encoding='utf-8')
    for line in fr.readlines():
        currLine=line.strip('\t').split()
        lineArr=[]
        nameMat.append([currLine[0],currLine[1]])
        for i in range(5):
            lineArr.append(float(currLine[i+2]))
        dataMat.append(lineArr)
    fr.close()
    print("Loading completely!")
    dataMat=mat(dataMat)
    nameMat=mat(nameMat)
    return dataMat, nameMat

def calVe(dataMat):
    e_mat=np.multiply(dataMat[:,0],dataMat[:,1])
    return e_mat

def calDp(dataMat):
    dp_mat=dataMat[:,2]+0.5*dataMat[:,3]
    return dp_mat

def getThetaE(dataMat):
    return dataMat[:,4]/100

def cal_Va_ThetaA(Ve_mat,dp_mat,r,t,ThetaE):  #core function--->V(A)  and Theta(A)
    VaArr=[]
    ThetaArr=[]
    Ve=Ve_mat.flatten().A[0]
    dp=dp_mat.flatten().A[0]
    ThetaE=ThetaE.flatten().A[0]
    i=0
    while i<len(Ve):
        Va,ThetaA=NLFG.f(Ve[i],ThetaE[i],dp[i],r,t)
        VaArr.append(Va)
        ThetaArr.append(ThetaA)
        i += 1
    num=len(VaArr)
    VaArr=array(VaArr).reshape(num,1)
    ThetaArr=array(ThetaArr).reshape(num,1)
    return mat(VaArr),mat(ThetaArr)

def calDD(Va_mat,ThetaA_mat,dp_mat):
    DD_mat=np.multiply((Va_mat-dp_mat),np.power(np.multiply(Va_mat,ThetaA_mat),-1))
    for i in range(len(DD_mat)):
        if DD_mat[i]<-1.67:
            DD_mat[i]=-1.67
        elif DD_mat[i]>8:
            DD_mat[i]=1.67
    return DD_mat

def calEDF(DD_mat):
    m=len(DD_mat)
    EDF=zeros((m,1))
    EDF=mat(EDF)
    i=0
    while i<m:
        EDF[i]=NLFG.N(-DD_mat[i])
        i += 1
    return EDF

def plz_DD_EDF(X,Y,index):
    fig1=plt.figure('DD_EDF')
    ax=fig1.add_subplot(111)
    ax.set_title('DD-EDF Graph')
    plt.xlabel('DD')
    plt.ylabel('EDF/%')
    if index==1:
        cc='b'
    elif index==2:
        cc='k'
    else:
        cc='y'
    ax.scatter(X.flatten().A[0],Y.flatten().A[0]*100,c=cc,marker='o',s=10)
    return

def plz_DD_Smp(Y,index,num):
    fig2=plt.figure('DD_Samples')
    ax=fig2.add_subplot(111)
    ax.set_title('DD-Samples Graph')
    plt.xlabel('Samples Index')
    plt.ylabel('DD')
    if index == 1:
        mk = '-ob'
    elif index == 2:
        mk = '-^r'
    else:
        mk = '-sk'
    X=range(num)
    plt.plot(X,Y.flatten().A[0],mk,markersize=1)
    return

def calall(inputfilename,outputfilename):
    r=0.015  # risk-free interest rate
    t=1 # duration
    dm,nm=loadDataset(inputfilename)
    Ve_mat=calVe(dm)
    dp_mat=calDp(dm)
    ThetaE_mat=getThetaE(dm)
    Va_mat,ThetaA_mat=cal_Va_ThetaA(Ve_mat,dp_mat,r,t,ThetaE_mat)
    dd_mat=calDD(Va_mat,ThetaA_mat,dp_mat)
    edf_mat=calEDF(dd_mat)
    fw=open(outputfilename,'w')
    fw.write('EquityValue'+'\t\t\t'+'ThetaEquity'+'\t\t'+'DefaultPoint'+'\t\t'+'AssetValue'+'\t\t'+'ThetaAsset'
             +'\t\t'+'DefaultDistance'+'\t\t\t'+'Edf'+'\n')
    for i in range(0,len(dm)):
        Ve=Ve_mat[i].flatten().A[0][0]
        Te=ThetaE_mat[i].flatten().A[0][0]
        dp=dp_mat[i].flatten().A[0][0]
        Va=Va_mat[i].flatten().A[0][0]
        Ta=ThetaA_mat[i].flatten().A[0][0]
        dd=dd_mat[i].flatten().A[0][0]
        edf=edf_mat[i].flatten().A[0][0]
        fw.write(str(Ve) + '\t\t' + str(Te) + '\t\t' + str(dp) + '\t\t' + str(Va) + '\t\t'
                 + str(Ta) + '\t\t' + str(dd) + '\t\t' + str(edf)+'\n')
    fw.close()
    return dd_mat,edf_mat

if __name__ == '__main__':
    #dd_mat,edf_mat=calall('2016Data.txt','2016Data_output.txt')
    #dd_mat,edf_mat=calall('2015Data.txt','2015Data_output.txt')
    #dd_mat,edf_mat=calall('2014Data.txt','2014Data_output.txt')
    #dd_mat,edf_mat=calall('2016AllStocksData.txt','2016AllStocksData_output.txt')
    dd_mat,edf_mat=calall('2015AllStocksData.txt','2015AllStocksData_output.txt')
    #dd_mat,edf_mat=calall('2014AllStocksData.txt','2014AllStocksData_output.txt')
    '''plz_DD_EDF(dd_mat,edf_mat,1)
    x=linspace(-2,8,10000)
    xm=mat(x).reshape(len(x),1)
    ym=calEDF(xm)
    y=ym.flatten().A[0]*100
    plt.plot(x,y,'r-')
    plz_DD_Smp(dd_mat,1,2724) # The forth parameter is the number of samples
    plt.show()'''