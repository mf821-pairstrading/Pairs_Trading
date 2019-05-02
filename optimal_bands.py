#Creator: Yuze Bai
#MF 821 Project
#Optimal bands

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

def c(n):
    df = pd.read_csv('data.csv', usecols=range(8*n - 6, 8*n)).dropna()
    M1 = sum(df.iloc[:,1]-df.iloc[:,0])/len(df)
    M2 = sum(df.iloc[:,3]-df.iloc[:,2])/len(df)
    c = (M1 + M2)/4
    return c

def HPDE(T,sigma,keppa,theta,smin,smax,Nt,Ns,c):
    ht = T/Nt
    hs = (smax-smin)/Ns
    
    s = np.arange(smin, smax+hs, hs)
    a = np.full(len(s),1 - sigma**2*ht/hs**2)
    u = 0.5 * (keppa*(theta-s)*ht/hs + sigma**2*ht/hs**2)
    l = 0.5 * (sigma**2*ht/hs**2- keppa*(theta-s)*ht/hs)
    A = np.diag(a[1:Ns])
    
    l1 = l[2:Ns]
    u1 = u[1:Ns-1]
    
    for i in range(len(l1)):
        A[i+1][i] = l1[i]
        A[i][i+1] = u1[i]
    
    
    hend = s - c
    h1 = hend[1:Ns]
    h2 = h1
    H = np.zeros((Ns-1, Nt))
    band = []

    for i in range(Nt):
        h2 = A.dot(h2)
        h2[-1] = h2[-1] + u[-2]*(smax - c)
        h2[0] = h2[0] + l[1] *(smin-c)
        h3 = []
        for j in range(len(h2)):
            if h2[j] < h1[j]:
                h3.append(h1[j])
                h2[j] = h1[j]
            H[j][i] = h2[j]
        band.append(min(h3) + c)
    list.reverse(band)
    r = [H,band]        
    return r

def GPDE(T,sigma,keppa,theta,smin, smax, Nt,Ns,c,H):
    ht = T/Nt
    hs = (smax-smin)/Ns
    
    s = np.arange(smin, smax+hs, hs)
    a = np.full(len(s),1 - sigma**2*ht/hs**2)
    u = 0.5 * (keppa*(theta-s)*ht/hs + sigma**2*ht/hs**2)
    l = 0.5 * (sigma**2*ht/hs**2- keppa*(theta-s)*ht/hs)
    A = np.diag(a[1:Ns])
    
    l1 = l[2:Ns]
    u1 = u[1:Ns-1]
    
    for i in range(len(l1)):
        A[i+1][i] = l1[i]
        A[i][i+1] = u1[i]
    
    
    hend = [-2 * c for x in range(Ns)]
    h2 = hend[1:Ns]
    G = np.zeros((Ns-1, Nt))
    band = []
    hr = s - c
    hr = hr[1:Ns]
    for i in range(Nt):
        h2 = A.dot(h2)
        h2[-1] = h2[-1] + u[-2]*(smax - c)
        h2[0] = h2[0] + l[1] *(smin-c)
        h3 = []
        for j in range(Ns-1):
            if h2[j] < H[j][i]-hr[j]:
                h3.append([H[j][i]-hr[j],H[j][i]])
                h2[j] = H[j][i]-hr[j]
            G[j][i] = h2[j]
        a = min(h3)
        band.append(a[1]-a[0]-c)
    list.reverse(band)
    r = [G,band] 
    return r 


for n in range(1):
    
    sigma = [3.541689226183422, 3.311109421758028, 2.3639786461316663, 2.6766354554911476, 3.4596600694045962, 7.293134132042836, 3.261049714507078]
    
    kappa = pd.read_csv('coi_drifts.csv').iloc[n, 1]
    theta = pd.read_csv('coi_thetas_solveu_3.csv').iloc[n, 1] 
    
    r = HPDE(1,sigma[n],kappa,theta,theta-10*sigma[n],theta+10*sigma[n],39000,3000,c(n+1))
    H = r[0]
    band1 = r[1]
    r2 = GPDE(1,sigma[n],kappa,theta,theta-10*sigma[n],theta+10*sigma[n],39000,3000,c(n+1),H)
    G = r2[0]
    band2 = r2[1]

    plt.plot(band1)
    plt.plot(band2)
    plt.legend(['Exiting Band', 'Entering Band'])
    plt.show()    

    exitband = []
    enterband = []
    for i in [x*100 for x in range(390)]:
        exitband.append(band1[i])
        enterband.append(band2[i])
    

    csvData = [[x, y] for x, y in zip(enterband, exitband)]
    with open('optimal_bands_'+str(n+1)+'.csv', 'w') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(csvData)
    
    csvFile.close()


#print(H[:,-1])
#print(c(1))
#print(band)




  
    
    
        
        
        