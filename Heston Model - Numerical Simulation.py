# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 14:37:26 2018

@author: volle
"""

import random
from scipy.stats import norm,t
from math import exp,log
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time
import textwrap

days = 360
gamma1 = 0.5 #averaging factors for the discretivazion of Xt
gamma2 = 0.5 #averaging factors for the discretivazion of Xt
phiC = 1.5 #Threshold for the initiation of the two aproximate distribution of V(t+1 | Vt)
InitialP = 100  #DELETE TILL pricepath = []
Strike = 100
Vol = (0.18)**2   
T = 1
Q = 0
Riskfree = 0.02
mu = 0.05
step = 1/360
sim = 100
accuracy = 0.0001
barriertype = 0
barrierlevel = 0
rho = -0.15
theta = (0.25)**2
k = 2
eta= 0.5

def validinput(prompt, type_=None, min_=None, max_=None, range_=None,base_=None):
    if min_ is not None and max_ is not None and max_ < min_:
        raise ValueError("min_ must be less than or equal to max_.\r")
    while True:
        ui = input(prompt)
        if type_ is not None and base_ is not None:
            try:
                ui = type_(ui)
            except ValueError:
                ui = base_
                print('Input not valid. Variable set to default value: %i' %ui)
                        
        if type_ is not None:
            try:
                ui = type_(ui)
            except ValueError:
                print("Input type must be {0}.\r".format(type_.__name__))
                continue
        if max_ is not None and ui > max_:
            print("Input must be less than or equal to {0}.\r".format(max_))
        elif min_ is not None and ui < min_:
            print("Input must be greater than or equal to {0}.\r".format(min_))
        elif range_ is not None and ui not in range_:
            if isinstance(range_, range):
                template = "Input must be between {0.start} and {0.stop}.\r"
                print(template.format(range_))
            else:
                template = "Input must be {0}.\r"
                if len(range_) == 1:
                    print(template.format(*range_))
                else:
                    print(template.format(" or ".join((", ".join(map(str,
                                                                     range_[:-1])),
                                                       str(range_[-1])))))
        else:
          return ui
      
def payoff(strikeprice,price,maxprice,minprice,callput,barriertype,barrierlevel):
    'Calculate the payoff based on input for Plain Vanilla and Barrier Option'
    k = strikeprice
    x = price
    if barriertype == 0: #Plain Vanilla Option
            if callput==1: #PutPayoff
                payoff = float(max(0, k - x ))
            elif callput==0: #CallPayoff     
                payoff = float(max(0, x - k ))
            return payoff
    elif barriertype == 1: #Knock-In Option
         if barrierlevel > strikeprice and maxprice > barrierlevel:  #Up-and-in
             if callput==1: #PutPayoff
                payoff = float(max(0, k - x ))
             elif callput==0: #CallPayoff     
                payoff = float(max(0, x - k ))
             return payoff
         if barrierlevel < strikeprice and minprice < barrierlevel: #Down-and-in
             if callput==1: #PutPayoff
                payoff = float(max(0, k - x ))
             elif callput==0: #CallPayoff     
                payoff = float(max(0, x - k ))
             return payoff
         else: 
            payoff = 0
            return payoff
    
    elif barriertype == 2: #Knock-Out Option
         if barrierlevel > strikeprice and maxprice < barrierlevel:  #Up-and-out
             if callput==1: #PutPayoff
                payoff = float(max(0, k - x ))
             elif callput==0: #CallPayoff     
                payoff = float(max(0, x - k ))
             return payoff
         if barrierlevel < strikeprice and minprice > barrierlevel: #Down-and-out
             if callput==1: #PutPayoff
                payoff = float(max(0, k - x ))
             elif callput==0: #CallPayoff     
                payoff = float(max(0, x - k ))
             return payoff
         else: 
            payoff = 0
            return payoff


def variancenew (Vt) :      #QE discretization algorithm for the Variance Process {Vt+1}, yield one step ahead variance conditional on Vt
    m = theta+( Vt - theta) * exp(-k*step)
    s2 = ((Vt*eta**2) * exp(-k*step)/k)*(1-exp(-k*step)) + ( theta*eta**2)*((1-exp(-k*step))**2)/(2*k)
    phi = s2/m**2
    
    if phi <= phiC:
        b2 = (2*phi**(-1))-1+2*((2*phi**(-1))**0.5)*((-1+2*phi**(-1))**0.5) 
        a = m/(1+b2)
        Zv = norm.ppf(random.random())
        Vnew = a*(Zv + b2**0.5)**2    #Non central Chi square variable aproximate sufficiently big value of Vt
        return Vnew
    elif phi > phiC:
        p = (phi-1)/(phi+1)
        beta = 2/ (m+m*phi)           #Function of Delta Dirac variable for sufficiently small value of Vt
        Uv = random.random()
        
        if Uv <=p:
            Vnew=0
            return Vnew
        elif Uv> p:
            Vnew= (beta**(-1)) * log((1-p)/(1-Uv))
            return Vnew
        
def returnnew (Vt,Xt):    #QE discretization Agorithm for the Return Process {Xt+1}, yield one step ahead return conditional on Vt and Xt
    Vnew = variancenew(Vt)   
    K0 = -1*(rho*k*theta)*step/eta
    K1 = gamma1*step*(-0.5+(k*rho/eta))-(rho/eta)
    K2 = gamma2*step*(-0.5+(k*rho/eta))+(rho/eta)
    K3 = gamma1*step*(1-rho**2)
    K4 = gamma2*step*(1-rho**2)
    epsilon = norm.ppf(random.random())
    
    Xnew = Xt + Riskfree*step + K0 + K1*Vt + K2*Vnew + ((K3*Vt+K4*Vnew)**0.5)*epsilon
    VarRet = [Vnew,Xnew]
    return VarRet


def path():
    path = [[],[],[InitialP]]
    i=0
    x = Vol**2     #Initiation Value for the Dynamic Variance Process
    y = Riskfree   #Initiation Value for the Dynamic Return Process
    Snew= InitialP
    while i < (T/step): 
        new = returnnew(x,y)
        Snew = InitialP * exp(new[1])
        path[0].append(new[0]**0.5)
        path[1].append(new[1])
        path[2].append(Snew)
        x,y = new[0],new[1]   #Update Vt and Xt for recursive computation of Vt+1 and Xt+1
        i = i+1
    VT = Snew
    minPrice= min(path[2])
    maxPrice= max(path[2])
    
    return [VT,maxPrice, minPrice, path]  #VT:Final Price  minPrice, maxPrice: max and min of the path, path: list of all Price value along the path

def simulation(): #Run the simulation
    start = time.time()
    underprice = 0
    optionpayoff = 0
    meanoption=0
    disprice = 0
    x = 0
    timer = 0
    check = 0
    ask = 0
    vector = [[],[],[],[],[],[],[],[]]
    timeseries = [[],[],[],[]]
     
    kind = int(validinput('Type 0: Call    1: put    ---> ',int,range_=(0,1)))
    
    portion = 0.1
    
    while check < 10:  #timer < sim:   #Loop for the simulation
         newprice = path()
         underprice = float(newprice[0])   #Calculate new price
        
         optionpayoff = float(payoff(Strike,underprice,newprice[1],newprice[2],kind,barriertype,barrierlevel )) #calculate option payoff
        
         meanoption = (meanoption*timer + optionpayoff)/(timer+1) #Mean optionpayoff, updates at each simulation
         disprice = float((meanoption * exp(- Riskfree * T)))  #Discount the average payoff
         vector[0].append(timer)           #Create a matrix with data of each simulation
         vector[1].append(underprice)
         vector[2].append(optionpayoff)
         vector[3].append(meanoption)
         vector[4].append(disprice)
         if timer > 105:
             vector[5].append((1-np.std(vector[4][-100:])/disprice)*100)
             vector[6].append((1-np.std(vector[4][-round(len(vector[4])*0.1):])/np.average(vector[4][-round(len(vector[4])*0.1):]))*100)
             #vector[7].append(abs(disprice-np.average(vector[4][-round(len(vector[4])*0.1):]))/np.std(vector[4][-round(len(vector[4])*0.1):])) 
             vector[7].append(1-np.std(vector[4][-round(len(vector[4])*0.05):])/np.std(vector[4][-round(len(vector[4])*0.2):])) 
                         #Need to convert the value collected in vector[7] by a function : f(0) = 1 and f'(x) < 0 for for x > 0 ???
                         # the information given arise from the slope of the plot. If the quantity above increase then accuracy decrease, and viceversa
                         # only give information about how the simulation is doing, if improving accuracy or not
         timeseries[0].append(newprice[3][2])
         timeseries[1].append(newprice[3][0])
         timeseries[2].append(newprice[3][1])
         timer = timer + 1
         
         if timer % (sim*0.1) == 0 and timer <= sim: #Give an idea on the current state of the simulation
             percent = (timer/sim)*100
             print('# of Simulations: %i \r' %timer)
     
         
         if timer % (sim*portion) == 0 and timer > 100: #Check if after portion% of the # of simulation, accuracy level is reached
            if np.std(vector[4][-100:])/disprice < accuracy:
                ask = validinput('Accuracy Level reached. Do you want to continue? Y=1 /N=0   \r',range_=('0','1','Y','N','n','y' ))
                if  ask == '1' or ask == 'Y' or ask == 'y':
                    check = check + 1
                    portion = portion + 0.1
                    print('Resuming Simulation')
                    
                elif  ask == '0' or ask == 'n' or ask == 'N':
                    check = 15
                    
                else:
                    check = check + 1
                    
            elif timer>sim:  #ask to stop S. if reaches the requested # of S. even if accuracy is not Reached
                 ask1 = validinput('Number of Simulation reached BUT Accuracy Level NOT reached. Do you want to continue? Y=1 /N=0   \r',range_=('0','1','Y','N','n','y' ))
                 if  ask1 == '1' or ask1 == 'Y' or ask1 == 'y':
                    check = check + 1
                    portion = portion + 0.1
                    
                 elif ask1 == '0' or ask1 == 'n' or ask1 == 'N':
                    check = 15
                    print('Simulation Finished. Elaborating results')
                    
    end = time.time()
        
    probability = (1 - vector[2].count(0)/timer)*100
    print('+' * 50)
    print('        OptionPrice: ' + str(round(disprice,3)) + ' \n        Probability of Exercise: ' + str(round(probability,3)) +' %  ' + ' \n        Accuracy(on last 100 obs): ' + str(round(vector[5][-1],3)) + ' \n        Accuracy(on last 10%): ' + str(round(vector[6][-1],3)) + '  \n        Time: ' + str(round((end - start),3)) + ' seconds \n        Speed: ' + str(round((end-start)/sim*1000,3)) + ' ms per simulation')
    print('+' * 50 + '\n')
    
    
    plt.figure(0)
    #Plot all the simulated Path
    for i in range(0,int(timer),int(timer/100)):
        matplotlib.pyplot.plot (range(len(timeseries[0][i])),timeseries[0][i],'g-',linewidth=0.1) 
    
    if barrierlevel != 0:
        matplotlib.pyplot.plot([0, len(timeseries[0][1])], [barrierlevel, barrierlevel], color='r', linestyle='-', linewidth=1,label='Barrier Level')
    
    
    #Plot Barrier and Strike level on plot above
    matplotlib.pyplot.plot([0, len(timeseries[0][1])], [Strike, Strike], color='b', linestyle='-', linewidth=1,label='Strike Price')
    plt.title('Path Simulations:' + str(timer) + 'Paths Shown: 100')
    plt.legend()
    plt.ylabel('Underlying Price')
    plt.xlabel('Time')
    
    #Plot all the Underlying prices reached during simulation
    plt.figure(1)     
    plt.subplot(2, 1, 1)     
    matplotlib.pyplot.plot (vector[0],vector[1],'g-',label='Underlying Price',linewidth=0.5)
    plt.title('Simulations:' + str(timer))
    plt.legend()
    plt.ylabel('Underlying Price in T')
    
    #Plot all the Payoff generated during the simulation
    plt.subplot(2, 1, 2)
    matplotlib.pyplot.plot (vector[0],vector[2],'r-',label='Option Payoff',linewidth=0.5)
    plt.ylabel('Option Payoff in T')
    plt.legend()
    plt.xlabel('Number of Simulations')
    
    #Plot Option Price and Accuracy Level
    plt.figure(2)
    plt.subplot(2, 1, 1)
    matplotlib.pyplot.plot (vector[0][-round(0.1*len(vector[0])):],vector[4][-round(0.1*len(vector[4])):],'b-',label='Option Price %f' %vector[4][-1],linewidth=0.5)
    plt.ylabel('Option Price')
    plt.title('Price and accuracy for the last 10% of Simulations')
    plt.legend()
    plt.xlabel('Number of Simulations')
   
    plt.subplot(2, 1, 2)
    matplotlib.pyplot.plot (vector[5][-round(0.1*len(vector[5])):],'r-',label='Accuracy %f' %vector[5][-1],linewidth=0.5)
    plt.ylabel('Accuracy Level')
    plt.legend()
    plt.xlabel('Number of Simulations')
    
    plt.figure(3)
    plt.subplot(2, 1, 1)
    matplotlib.pyplot.plot (vector[6][-round(0.1*len(vector[6])):],'r-',label='Accuracy %f' %vector[6][-1],linewidth=0.5)
    plt.title('Accuracy Level based on last 10% of simulations')
    plt.ylabel('Accuracy Level')
    plt.legend()
    plt.xlabel('Number of Simulations')
    
    """
    acc1 = np.array(vector[7][-round(0.1*len(vector[7])):])     
    acc_1 = np.array(vector[7][-1-round(0.1*len(vector[7])):-1])
    accvec = -(acc1 - acc_1)/acc_1                               #Accuracy indicators tests, may be eliminated
    """
    
    plt.subplot(2, 1, 2)
    matplotlib.pyplot.plot (vector[7],'b-',label='Improvement: %f' %vector[7][-1],linewidth=0.8)
    
    plt.ylabel('Accuracy Performance %')
    plt.legend()
    plt.xlabel('Number of Simulations')
   
    plt.figure(4)
    plt.subplot(2, 1, 1)
    #Plot the Variance Process
    for i in range(0,int(timer),int(timer/100)):
        matplotlib.pyplot.plot (range(len(timeseries[1][i])),timeseries[1][i],'g-',linewidth=0.1) 
    
    plt.subplot(2, 1, 2)
    #Plot all the Return Process
    for i in range(0,int(timer),int(timer/100)):
        matplotlib.pyplot.plot (range(len(timeseries[2][i])),timeseries[2][i],'g-',linewidth=0.1) 
    
def disclamer():
    print('+'*60 +'\n')
    message = (textwrap.wrap('NOTE: Accuracy Performance is an indicator of the improvement in the quality of the estimate of the option price. The indicator is built as the % difference between std(n=20,option price) and std(n=5,option price). So any positive value indicates a corrispondent decrease in the standard deviation of the estimated price by the simulations',width=58))
    for i in range(len(message)):
        print(' ' + message[i])
    print('\n  FORMULA :  1-  \n             std(last 20% of obs,\'OptionPrice\') / \n             std(last  5% of obs,\'OptionPrice\')  \n')
    
    message = (textwrap.wrap('Simulations shall stop when this indicator remain positive and close to zero. If the indicator is negative, simulations are losing accuracy. If positive, simulations are improving accuracy. If null, simulation are not improving accuracy anymore.',width=58))
    for i in range(len(message)):
        print(' ' + message[i])
    print('\n' +'+'*60 +'\n')
    

simulation()
time.sleep(2)
disclamer()





