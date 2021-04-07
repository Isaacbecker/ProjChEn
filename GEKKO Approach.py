#%% - GEKKO Approach
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from scipy.integrate import odeint

for Section in range(1): #Input parameters
    
    T0=200+273.15 #[K] Initial Temperature
    P0=51 #[bar] Initial pressure
    NI=1000 #Number of intervals in the model
    FDME0=0 #[mol/s]
    FMetOH0=0 #[mol/s]
    FH2O0=0 #[mol/s]
    FH20=3 #[mol/s]
    FCO20=1 #[mol/s]
    FCO0=1 #[mol/s]

for Section in range(1): #Kinetic Parameters
    """
    RKinetics=8.3145 #[J/molK]
    k1pre=1.50E9
    k2pre=6.08E9
    k3pre=7.04E5
    k4pre=8.28E11
    KCOpre=7.99E-7
    KCO2pre=1.02E-7
    KH2OH2pre=4.13E-11
    KMetOHpre=7.9E-7
    KH2Opre=0.84E-4
    
    k1B=-125757
    k2B=-116370
    k3B=-88777
    k4B=-122427
    KCOB=58100
    KCO2B=67400
    KH2OH2B=104500
    KMetOHB=70500
    KH2OB=41100
    """
    
    RKinetics=8.3145 #[J/molK]
    k1pre=2.69E7
    k2pre=7.31E8
    k3pre=4.36E2
    k4pre=1.0278E10
    KCOpre=7.99E-7
    KCO2pre=1.02E-7
    KH2OH2pre=4.13E-11
    KMetOHpre=7.9E-7
    KH2Opre=0.84E-4
    
    k1B=-109900
    k2B=-123400
    k3B=-65200
    k4B=-105000
    KCOB=58100
    KCO2B=67400
    KH2OH2B=104500
    KMetOHB=70500
    KH2OB=41100
    
for Section in range(1): #Catalyst parameters
    W0=0 #[Kg]
    Wf=1000 #[Kg]
    W=np.linspace(W0,Wf,NI)

for Section in range(1): #Prelimiary calculations
    
    R=8.3145E-5 #[bar*m³/molK]
    FT0=FDME0+FMetOH0+FH2O0+FH20+FCO20+FCO0 #Initial total flowrate
    v0=FT0*R*T0/P0
    
#Initialize model
m=GEKKO()
m.time=W

for Section in range(1): #Set variables
    FDME=m.Var(value=FDME0,lb=0)
    FMetOH=m.Var(value=FMetOH0,lb=0)
    FH2O=m.Var(value=FH2O0,lb=0)
    FH2=m.Var(value=FH20,lb=0)
    FCO2=m.Var(value=FCO20,lb=0)
    FCO=m.Var(value=FCO0,lb=0)
    
    T=m.Var(value=T0,lb=0) #Optimize - Isothermal
    P=m.Var(value=P0,lb=1) #Optimize

    FT=m.Var(value=FT0)
    v=m.Var(value=v0)

#for Section in range(1): #Set intermediates

p=m.Intermediate(P/P0)
k1=m.Intermediate(k1pre*m.exp(k1B/(RKinetics*T)))
k2=m.Intermediate(k2pre*m.exp(k2B/(RKinetics*T)))
k3=m.Intermediate(k3pre*m.exp(k3B/(RKinetics*T)))
k4=m.Intermediate(k4pre*m.exp(k4B/(RKinetics*T)))
KCO=m.Intermediate(KCOpre*m.exp(KCOB/(RKinetics*T)))
KCO2=m.Intermediate(KCO2pre*m.exp(KCO2B/(RKinetics*T)))
KH2OH2=m.Intermediate(KH2OH2pre*m.exp(KH2OH2B/(RKinetics*T)))
KMetOH=m.Intermediate(KMetOHpre*m.exp(KMetOHB/(RKinetics*T)))
KH2O=m.Intermediate(KH2Opre*m.exp(KH2OB/(RKinetics*T)))

Keq1=m.Intermediate(10**(5139/T-12.621))
Keq2=m.Intermediate(10**(-2073/T+2.029))
Keq3=m.Intermediate(Keq1*Keq2)
Keq4=m.Intermediate(m.exp(4019/T+3.707*m.log(T)-2.783E-3*T+3.8E-7*T**2-6.561E4/T**3-26.64))

PDME=m.Intermediate(P*(FDME/FT)*p)
PMetOH=m.Intermediate(P*(FMetOH/FT)*p)
PH2O=m.Intermediate(P*(FH2O/FT)*p)
PH2=m.Intermediate(P*(FH2/FT)*p)
PCO2=m.Intermediate(P*(FCO2/FT)*p)
PCO=m.Intermediate(P*(FCO/FT)*p)

fDME=m.Intermediate(PDME)
fMetOH=m.Intermediate(PMetOH)
fH2O=m.Intermediate(PH2O)
fH2=m.Intermediate(PH2)
fCO2=m.Intermediate(PCO2)
fCO=m.Intermediate(PCO)

CDME=m.Intermediate(FDME/v)
CMetOH=m.Intermediate(FMetOH/v)
CH2O=m.Intermediate(FH2O/v)
CH2=m.Intermediate(FH2/v)
CCO2=m.Intermediate(FCO2/v)
CCO=m.Intermediate(FCO/v)

#R1=m.Intermediate(k1*KCO*(fCO*fH2**(3/2)-fMetOH/(fH2**(1/2)*Keq1))/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O)))
R1=m.Intermediate(k1*KCO*(fCO*fH2**(2)-fMetOH/(Keq1))/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O)*(fH2**(1/2))))
R2=m.Intermediate(k2*KCO2*( fCO2*fH2-fH2O*fCO/Keq2 )/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O)))
#R3=m.Intermediate(k3*KCO2*( fCO2*fH2**(3/2)-fMetOH*fH2O/(fH2**(3/2)*Keq3) )/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O)))
R3=m.Intermediate(k3*KCO2*( fCO2*fH2**(3)-fMetOH*fH2O/(Keq3) )/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O)*(fH2**(3/2))))
R4=m.Intermediate(k4*KMetOH**2*( CMetOH**2-CH2O*CDME/(Keq4) )/( (1+2*(KMetOH*CMetOH)**(1/2)+KH2O*CH2O)**4 ))

#for Section in range(1): #Set Equations
    
m.Equation(FT==FDME+FMetOH+FH2O+FH2+FCO2+FCO)
m.Equation(v==v0*(FT/FT0)*(P0/P)*(T/T0))

m.Equation(FDME.dt()==R4/2)# →Obj: Max(FDME)
m.Equation(FMetOH.dt()==R1+R3-R4)
m.Equation(FH2O.dt()==R2+R3+R4/2)
m.Equation(FH2.dt()==-2*R1-R2-3*R)
m.Equation(FCO2.dt()==-R2-R3)
m.Equation(FCO.dt()==-R1+R2)

m.Equation(T.dt()==0) #Ideal
m.Equation(P.dt()==0) #Ideal    

#Set solver options
m.options.IMODE=8

#Solve model
m.solve()

#Plot resutls

plt.plot(W,FDME)
plt.plot(W,FMetOH)
plt.plot(W,FH2O)
plt.plot(W,FH2)
plt.plot(W,FCO2)
plt.plot(W,FCO)
plt.legend(['$F_{DME}','$F_{MetOH}','$F_{H2O}$','$F_{H2}','$F_{CO2}$','$F_{CO}$'])
plt.ylabel('$F_i$ [mol/s]')
plt.xlabel('Catalyst mass [Kg]')
plt.title('Molar flowrate vs. Catalyst Mass')
plt.show()

# %%
