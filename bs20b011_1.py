#imported libraried 
# numpy for making arrays and using some functions 
# scipy.interpolate to convinently interpolate the data of fixed points anywhere in the code 
# matplotlib for plotting the graphs of the data and using some functions of it

import matplotlib.pyplot as plt    
from scipy.interpolate import interp1d
import numpy as np  

#At constant Temperature = 313.14K
pData=[25.99, 26.751, 27.298, 27.48, 27.355, 26.895, 26.003] #ddbst data of pressure values
x1=[0.1282, 0.2354, 0.3685, 0.4932, 0.6143, 0.7428, 0.8656]  #liquid mol fraction of Benzene
y1=[0.1657, 0.2766, 0.3912, 0.495,0.5909, 0.6979, 0.8205] #vapor mol fraction of Benzene
x2=[0.8718, 0.7646, 0.6315, 0.5068, 0.3857, 0.2572, 0.1344] #liquid mol fraction of Cyclohexane
r1=3.19#Relative volume constant of Benzene
r2=3.97#Relative volume constant of Cyclohexane
q1=2.4#Relative surface area constant of Benzene
q2=3.01#Relative surface area constant of Cyclohexane
q_2=q2
q_1=q1
z=10
l1=z/2*(r1-q1)-(r1-1)
l2=z/2*(r2-q2)-(r2-1)

R=8.314 # in J⋅mol−1⋅K−1.
Tau21=0.825001001
Tau12=0.999000001
# antoine Coefficients A,B,C of both Benzene and Cyclohexne
A_b=4.01814	
B_b=1203.835
C_b=-53.226
A_c=4.017914
B_c=1203.79935
C_c=-53.19826
#empty lists to store the activity coefficient values and vapour pressure values
G1=[]
G2=[] 
newY=[]
P_313=[]

T=313.14

#calculation of saturated pressure values of the componets at the fixed temperature
pSat1=100*(10**((A_b- B_b/(C_b + T))))
pSat2=100*(10**((A_c- B_c/(C_c + T))))

#used for loop for find the P value and the activity coefficent at different molefraction
for i in range(len(x2)):
#calculation of activity coefficent models using UNIQUAC model formulae 
    G1.append(np.exp(np.log(r1/(x1[i]*r1 + x2[i]*r2)) + z*q1/2*np.log(q1*((x1[i]*r1 + x2[i]*r2))/(((x1[i]*q1 + x2[i]*q2))* r1)) + (x2[i]*r2/(x1[i]*r1 + x2[i]*r2))  *  (l1- (r1/r2*l2)) - q_1* (np.log((x1[i]*q_1/((x1[i]*q_1) + (x2[i]*q_2)))  +  (x2[i]*q_2/((x1[i]*q_1) + (x2[i]*q_2)))* Tau21)) +  ((x2[i]*q_2)/(x1[i]*q_1 + x2[i]*q_2))  * q_1*(Tau21/( (x1[i]*q_1/((x1[i]*q_1) + (x2[i]*q_2)))  +  (x2[i]*q_2/((x1[i]*q_1) + (x2[i]*q_2)))*Tau21)  - Tau12/((x2[i]*q_2/((x1[i]*q_1) + (x2[i]*q_2))) + (x1[i]*q_1/((x1[i]*q_1) + (x2[i]*q_2)))*Tau12 ))))
    G2.append(np.exp(np.log(r2/(x1[i]*r1 + x2[i]*r2)) + z*q2/2*np.log(q2*((x1[i]*r1 + x2[i]*r2))/(((x1[i]*q1 + x2[i]*q2))* r2)) + (x1[i]*r1/(x1[i]*r1 + x2[i]*r2)) * (l2- (r2/r1*l1)) - q_2* (np.log((x2[i]*q_2/((x1[i]*q_1) + (x2[i]*q_2)))  +  (x1[i]*q_1/((x1[i]*q_1) + (x2[i]*q_2)))* Tau12)) + ((x1[i]*q_1)/(x1[i]*q_1 + x2[i]*q_2))  *  q_2*(Tau12/( (x2[i]*q_2/((x1[i]*q_1) + (x2[i]*q_2)))  +  (x1[i]*q_1/((x1[i]*q_1) + (x2[i]*q_2)))*Tau12)  - Tau21/((x1[i]*q_1/((x1[i]*q_1) + (x2[i]*q_2))) + (x2[i]*q_2/((x1[i]*q_1) + (x2[i]*q_2)))*Tau21 ))))
    P_313.append((x1[i]*G1[i]*pSat1 + x2[i]*G2[i]*pSat2))
    newY.append((x1[i]*G1[i]*pSat1)/(x1[i]*G1[i]*pSat1 + x2[i]*G2[i]*pSat2))
   
#code for plotting graphs
figure, axis = plt.subplots(1, 2)

axis[0].set_xlabel("x,y(Benzene)[mol/mol]")
axis[0].set_ylabel("P[Kpa]")
axis[0].set_title("VLE P/xy-chart Benzene+ cyclohexane (BS20B017)")
axis[0].scatter(x1, pData, color="blue")
axis[0].plot(x1, pData, color="blue")
axis[0].scatter(y1, pData, color="blue")
axis[0].plot(y1, pData, color="blue")



axis[0].scatter(x1, P_313, color="black")
axis[0].plot(x1, P_313, color="black")
axis[0].scatter(newY, P_313, color="black")
axis[0].plot(newY, P_313, color="black")




axis[1].set_xlabel("x(b)[mol/mol]")
axis[1].set_ylabel("y(b)[mol/mol]")
axis[1].set_title("VLE y/x-chart Benxene+ cyclohexane (BS20B017)")
axis[1].scatter(x1, y1, color="green")
axis[1].plot(x1, y1, color="green")

plt.show()