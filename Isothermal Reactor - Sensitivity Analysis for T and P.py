#%% Isothermal Reactor - Temperature and Pressure Sensitivity Analysis
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import cm

for Section in range(1): #Set Sensitivity Analysis parameters
    
    NTSA=30 #Number of points for Temperatures Sensitivity Analysis
    T0SA=10+273.15 #Initial Temperature for Sensitivity Analysis
    TfSA=200+273.15 #Final Temperature for Sensitivity Analysis
    TSA=np.linspace(T0SA,TfSA,NTSA)


    NPSA=30 #NPSA=100 #Number of point in Pressure for Sensitivity Analysis
    P0SA=1 #Initial Pressure for Sensitivity Analysis
    PfSA=80 #Final Pressure for Sensitivity Analysis
    PSA=np.linspace(P0SA,PfSA,NPSA)

    ANS=np.zeros((NTSA*NPSA,3))

    Count=0

for j in range(NPSA):
        
    for i in range(NTSA):

        for Section in range(1): #Input parameters
            
            T0=TSA[i]+273.15 #[K] Initial Temperature
            P0=PSA[j] #[bar] Initial pressure
            NI=50000 #Number of intervals in the model
            FDME0=0 #[mol/s]
            FMetOH0=0 #[mol/s]
            FH2O0=0 #[mol/s]
            FH20=300 #[mol/s]
            FCO20=100 #[mol/s]
            FCO0=100 #[mol/s]

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
            Wf=40000 #[Kg]
            W=np.linspace(W0,Wf,NI)

        for Section in range(1): #Prelimiary calculations
            
            R=8.3145E-5 #[bar*m³/molK]
            FT0=FDME0+FMetOH0+FH2O0+FH20+FCO20+FCO0 #Initial total flowrate
            v0=FT0*R*T0/P0
            

        # function that returns dz/dt
        def model(z,W):
            
            for Section in range(1): #Assign Fi to z[i]
                FDME=z[0]
                FMetOH=z[1]
                FH2O=z[2]
                FH2=z[3]
                FCO2=z[4]
                FCO=z[5]
                T=z[6]
                P=z[7]
            
            for Section in range(1): #Other calculations
                
                p=P/P0
                FT=FDME+FMetOH+FH2O+FH2+FCO2+FCO
                v=v0*(FT/FT0)*(P0/P)*(T/T0)
                
            for Section in range(1): #Calculate kinetic and equilibrium expressions
                
                k1=k1pre*np.exp(k1B/(RKinetics*T))
                k2=k2pre*np.exp(k2B/(RKinetics*T))
                k3=k3pre*np.exp(k3B/(RKinetics*T))
                k4=k4pre*np.exp(k4B/(RKinetics*T))
                KCO=KCOpre*np.exp(KCOB/(RKinetics*T))
                KCO2=KCO2pre*np.exp(KCO2B/(RKinetics*T))
                KH2OH2=KH2OH2pre*np.exp(KH2OH2B/(RKinetics*T))
                KMetOH=KMetOHpre*np.exp(KMetOHB/(RKinetics*T))
                KH2O=KH2Opre*np.exp(KH2OB/(RKinetics*T))
                
                Keq1=10**(5139/T-12.621)
                Keq2=10**(-2073/T+2.029)
                Keq3=Keq1*Keq2
                Keq4=np.exp(4019/T+3.707*np.log(T)-2.783E-3*T+3.8E-7*T**2-6.561E4/T**3-26.64)

            for Section in range(1): #Calculate Pi (Partial pressure of component i)
                
                PDME=P*(FDME/FT)*p
                PMetOH=P*(FMetOH/FT)*p
                PH2O=P*(FH2O/FT)*p
                PH2=P*(FH2/FT)*p
                PCO2=P*(FCO2/FT)*p
                PCO=P*(FCO/FT)*p
            
            for Section in range(1): #Calculate fi (fugacity)
                
                fDME=PDME
                fMetOH=PMetOH
                fH2O=PH2O
                fH2=PH2
                fCO2=PCO2
                fCO=PCO
            
            for Section in range(1): #Calculate Ci (Concentration)
                
                CDME=FDME/v
                CMetOH=FMetOH/v
                CH2O=FH2O/v
                CH2=FH2/v
                CCO2=FCO2/v
                CCO=FCO/v
            
            for Section in range(1): #Set Ri (Reaction rate)
                
                R1=k1*KCO*(fCO*fH2**(3/2)-fMetOH/(fH2**(1/2)*Keq1))/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O))
                #R1=k1*KCO*(fCO*fH2**(2)-fMetOH/(Keq1))/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O)*(fH2**(1/2)))
                R2=k2*KCO2*( fCO2*fH2-fH2O*fCO/Keq2 )/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O))
                #R3=k3*KCO2*( fCO2*fH2**(3/2)-fMetOH*fH2O/(fH2**(3/2)*Keq3) )/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O))
                R3=k3*KCO2*( fCO2*fH2**(3)-fMetOH*fH2O/(Keq3) )/((1+KCO*fCO+KCO2*fCO2)*(fH2**(1/2)+KH2OH2*fH2O)*(fH2**(3/2)))
                R4=k4*KMetOH**2*( CMetOH**2-CH2O*CDME/(Keq4) )/( (1+2*(KMetOH*CMetOH)**(1/2)+KH2O*CH2O)**4 )
            
            dFDMEdW = R4/2
            dFMetOHdW = R1+R3-R4
            dFH2OdW = R2+R3+R4/2
            dFH2dW = -2*R1-R2-3*R
            dFCO2dW = -R2-R3
            dFCOdW = -R1+R2
            dTdW = 0
            dPdW = 0
            
            dzdW = [dFDMEdW,dFMetOHdW,dFH2OdW,dFH2dW,dFCO2dW,dFCOdW,dTdW,dPdW]
            return dzdW

        # initial condition
        z0 = [FDME0,FMetOH0,FH2O0,FH20,FCO20,FCO0,T0,P0]

        # solve ODE
        z = odeint(model,z0,W)
        
        print('Loopi=',i,'Loopj=',j,'     |     T=',TSA[i],'     |     ','P=',PSA[j])
        ANS[Count,0]=TSA[i]
        ANS[Count,1]=PSA[j]
        ANS[Count,2]=max(z[:,0])
        Count=Count+1
        

x = np.reshape(ANS[:,0], (NPSA, NTSA))
y = np.reshape(ANS[:,1], (NPSA, NTSA))
z = np.reshape(ANS[:,2], (NPSA, NTSA))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z,cmap=cm.coolwarm)
ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Pressure [bar]')
ax.set_zlabel('FDME [mol/s]')
plt.title('Sensitivity Analysis for FDME')

plt.show()


# %%
