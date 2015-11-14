# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 01:07:54 2015

@author: amine


=========================================================================
Simulation du modèle du Chemostat
=========================================================================
HADJ ABDELKADER Oussama
Le 07 Juin 2015 à Montpellier
-------------------------------------------------------------------------
F. Campillo, M. Joannides, I. Larramendy-Valverde, Stochastic modeling of
the Chemostat Ecological Modelling
B. Benyahia, F. Campillo, B. Cherki, J. Harmand, Particle filtering for
the Chemostat.
-------------------------------------------------------------------------
function y=Chemostat()
Cette fonction est une simulation des EDS du Chemostat par la méthode 
d'Euler Maruyama. La Sortie de cette fonction est un vecteur
d'observations Y celculées à partir de l'equation de sortie de ce système
=========================================================================
"""
from __future__ import print_function
import numpy as np
import math 
import matplotlib.pyplot as plt

T = 1000                 # T: temps total de simulation
Nob = 1000               # nombre d'observations
Ntk = 10                 # Ntk: nombre d'itérations de simulation
dt = float(T)/(Ntk*Nob)         # dt: Constante de temps
mu_m=0.3                 # mu_max: taux de croissance maximum
ks=10                    # ks: constante de demi saturation
D=0.01                   # D: taux de dilution
cs=10                    # k: coefficicent stoechiométrique
Sin=100                  # Sin: concentration de la biomasse à l'entrée
c1=0.03                  # ci: écarts-types des Bruits
c2=0.03                  # ci: écarts-types des Bruits
sigma=0.2                # sigma: écart-type du bruit d'observation
"""
code Matlab : à convertir en Python !!
% tirage des nombres aléatoires:
% dw1=randn(N,1);
% dw2=randn(N,1);
% v=randn(N,1);
"""

# allocation espace memoire
B = np.zeros(Ntk)
S = np.zeros(Ntk)
y = np.zeros(Nob)


# initialisation
S[0] = np.random.normal(2,2)
B[0] = np.random.normal(2,2)
S[0] = max(0, S[0])
B[0] = max(0, B[0])

Bt = []
St = []
# itérations avec vectorisation dans la premiere boucle
for k in range(0,Nob): 
    mu=mu_m * S[0:Ntk-1]/(ks+S[0:Ntk-1]) # mu: taux de croissance
    # Calcul des états
    B[1:] = np.maximum(np.zeros(Ntk-1), B[0:Ntk-1]+(mu[Ntk-2]-D)*B[0:Ntk-1]*dt+c1*np.sqrt(B[0:Ntk-1])*np.sqrt(dt)* np.random.normal(0,1,Ntk-1))
    S[1:] = np.maximum(np.zeros(Ntk-1) , S[0:Ntk-1]+((-cs*mu[Ntk-2]*B[0:Ntk-1])+(D*(Sin-S[0:Ntk-1])))*dt+c2*np.sqrt(S[0:Ntk-1])*np.sqrt(dt)*np.random.normal(0,1,Ntk-1) )
    
    Bt = np.concatenate([Bt, B])
    St = np.concatenate([St, S])

    B[0] = max(0 , B[Ntk-1]+(mu[Ntk-2]-D)*B[Ntk-1]*dt+c1*np.sqrt(B[Ntk-1])*np.sqrt(dt)*np.random.normal())
    S[0] = max(0 , S[Ntk-1]+((-cs*mu[Ntk-2]*B[Ntk-1])+(D*(Sin-S[Ntk-1])))*dt+c2*np.sqrt(S[Ntk-1])*np.sqrt(dt)*np.random.normal())
    
    # l'équation de sortie
    y[k] = S[Ntk-1]+sigma*S[Ntk-1]*np.random.normal();


# tracer les resultats
nbElements = (float(T-dt))/(10*Ntk*dt)
instants = np.linspace(0,T,nbElements+1)

bPlot = Bt[0:(Ntk*Nob):(10*Ntk)]
sPlot = St[0:(Ntk*Nob):(10*Ntk)]


fig1 = plt.figure()
plt.plot(instants,bPlot,'b')
#plt.xlim(0, 1)
#plt.ylim(0, 1)
plt.xlabel('t')
plt.ylabel('b(t)')
plt.title('Etat s du systeme')
plt.legend('b(t)')
plt.show()

fig2 = plt.figure()
plt.plot(instants,sPlot,'r')
#plt.xlim(0, 1)
#plt.ylim(0, 1)
plt.xlabel('t')
plt.ylabel('b(t)')
plt.title('Etat s du systeme')
plt.legend('s(t)')
plt.show()
