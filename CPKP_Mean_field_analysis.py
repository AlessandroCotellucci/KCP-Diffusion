import matplotlib
from matplotlib import pyplot as plt
import vegas
import math
import random
from numpy import *

def Deriv(x,y,p):
    #Parameters of the model
    mu_c=0.05     #Colonized discharge rate
    mu_u=0.1      #Uncolonized discharge rate
    alpha=0.21    #Probability of colonization
    phi=0.05      #Probability that the patient added is colonized
    lambd=1       #Probability to add a patient
    mu_HC=0.5     #Probability of a colonized HCW to become uncolonized
    #p=0
    b=0.1
    c=0.5
    B=80
    dy=zeros(n,'double')

    dy[0]=(1-phi)*lambd*(B-y[0]-y[1])-mu_u*y[0]-alpha*b*y[0]*y[3]
    dy[1]=phi*lambd*(B-y[0]-y[1])-mu_c*y[1]+alpha*b*y[0]*y[3]
    dy[2]=mu_HC*y[3]-alpha*b*y[1]*y[2]-alpha*c*y[2]*y[3]
    dy[3]=-mu_HC*y[3]+alpha*b*y[1]*y[2]+alpha*c*y[2]*y[3]
    return dy


def Derivmod(x,y,p):
    #Parameters of the model
    mu_c=0.05     #Colonized discharge rate
    mu_u=0.1      #Uncolonized discharge rate
    alpha=0.21    #Probability of colonization
    phi=0      #Probability that the patient added is colonized
    lambd=1       #Probability to add a patient
    mu_HC=0.5     #Probability of a colonized HCW to become uncolonized
    #p=0.7
    b=0.1
    c=0.5
    B=80
    dy=zeros(n,'double')

    dy[0]=(1-phi)*lambd*(B-y[0]-y[1])-mu_u*y[0]-alpha*b*y[0]*y[3]
    dy[1]=phi*lambd*(B-y[0]-y[1])-mu_c*y[1]+alpha*b*y[0]*y[3]
    dy[2]=mu_HC*y[3]-(1-p)*alpha*b*y[1]*y[2]-(1-p)*alpha*c*y[2]*y[3]
    dy[3]=-mu_HC*y[3]+(1-p)*alpha*b*y[1]*y[2]+(1-p)*alpha*c*y[2]*y[3]
    return dy



def RK4(x,y,n,h,xf,p):
    k1=zeros(n, 'double')
    k2=zeros(n, 'double')
    k3=zeros(n, 'double')
    k4=zeros(n, 'double')
    ym=zeros(n, 'double')
    ye=zeros(n, 'double')
    slope=zeros(n, 'double')

    if x<xf/2:
        k1=Deriv(x,y,p)
    else:
        k1=Derivmod(x,y,p)

    for i in range(n):
        ym[i]=y[i]+k1[i]*h/2

    if x<xf/2:
        k2=Deriv(x+h/2,ym,p)
    else:
        k2=Derivmod(x+h/2,ym,p)

    for i in range(n):
        ym[i]=y[i]+k2[i]*h/2

    if x<xf/2:
        k3=Deriv(x+h/2,ym,p)
    else:
        k3=Derivmod(x+h/2,ym,p)

    for i in range(n):
        ye[i]=y[i]+k3[i]*h
    if x<xf/2:
        k4=Deriv(x+h,ye,p)
    else:
        k4=Derivmod(x+h,ye,p)

    for i in range(n):
        slope[i]=(k1[i]+2*(k2[i]+k3[i])+k4[i])/6
        y[i]=y[i]+slope[i]*h
    x=x+h
    return x, y

def integrator(x,y,n,h,xend,xf,p):
    while x<xend:
        if xend-x<h:
            h=xend-x
        x,y=RK4(x,y,n,h,xf,p)
    return x ,y



def Simulation(p):
    n=4
    y=zeros(n,'double')
    yi=zeros(n,'double')
    yi[0]=80
    yi[1]=0
    yi[2]=10
    yi[3]=0
    xi=0
    dx=0.01
    xf=1000
    xout=1

    j=int(dx*xf)
    xp=zeros(1001,'double')
    yp=zeros((n,1001), 'double')


    x=xi
    m=0
    xp[m]=x
    for i in range(n):
        yp[i,m]=yi[i]
        y[i]=yi[i]

    while x<xf:
        xend=x+xout
        if xend>xf:
            xend=xf
        h=dx
        x,y=integrator(x,y,n,h,xend,xf,p)
        m=m+1
        xp[m]=x
        for i in range(n):
            yp[i,m]=y[i]

    Num_col=(abs(yp[1,1000])+abs(yp[3,1000]))/(abs(yp[0,1000])+abs(yp[1,1000])+abs(yp[2,1000])+abs(yp[3,1000]))
    return Num_col


n=21
p=zeros(n,'double')
Perc=zeros(n,'double')
for j in range(n):
    p[j]=j*1/(n-1)
    Perc[j]=Simulation(p[j])



plt.plot(p,Perc,'r')
plt.xlabel(r'$p$')
#plt.xlabel(r'$\alpha$')
plt.ylabel('Fraction of colonized nodes')
plt.title('KCP Diffusion in mean field approximation')
plt.show()
