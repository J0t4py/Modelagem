import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy import integrate
from scipy.spatial.transform import Rotation
import control as ct
#constantes
I_m = 3*(10**-4) #Inércia do motor em kg.m^2
I_1 = 3*(10**-4) #Inércia do motor em kg.m^2
I_2 = 1.5 #Inércia da engrenagem em kg.m^2
I_L = 50 #Inércia da carga em kg.m^2
k_L = 5*(10**7) #Rigidez torsional da carga em Nm/rad
k_m = 1.2*(10**3) #Rigidez torsional do motor em Nm/rad
c_L = 10 #coeficiente de amortecimento do rolamento do munhão em Nms/rad
c_m = 0.3 #coeficiente de amortecimento do motor
K = 60 #constante de torque do motor em Nm/Amp  #TESTAR COM 60
R_a = 0.8 #Resistência elétrica do motor em Ohms
L_a = 1*10**-3 #Indutância do motor em henry
Np = 14 #número de dentes do pinhão
Ne = 408 #número de dentes do munhão
n = 29.14#Ne/Np # Relação de transmissão

# Condições iniciais
# Ângulo inicial do motor
theta_m_0 = 0.02
# Ângulo inicial do pinhão
theta_1_0 = 0.03
# Ângulo inicial da carga
theta_L_0 = 0.04
# Corrente no instante inicial
i_a_0 = 0.00
# Velocidade angular inicial do motor
theta_m_dot_0 = 0.00
# Velocidade angular inicial do pinhão
theta_1_dot_0 = 0.00
# Velocidade angular inicial da carga
theta_L_dot = 0.00
# vetor com as condições iniciais
y_0 = np.array([theta_m_0, theta_1_0, theta_L_0, theta_m_dot_0, theta_1_dot_0, theta_L_dot, i_a_0])
#y_0 = y_0.reshape((8,1))
Q_deg_lista = ([]) # Esforço aplicado à torreta - entrada em degrau
e_a_lista = ([])
#calcula as entradas para cada instante de tempo
for t in sol.t:
  e_a_lista.append(32*np.sin(2*np.pi*60*t))
  if 2 < t < 5:
    Q_deg_lista.append(3.8)
  elif 6.5 < t < 9:
    Q_deg_lista.append(-3.8)
  else:
    Q_deg_lista.append(0)

t_eval=np.linspace(0, 10, 100)

#Plot das entradas
plt.figure(3)
tempo=np.linspace(0, 10, 100)
plt.plot(tempo, Q_deg_lista)
plt.title("Entrada - Carga em Degrau (N)")
plt.legend(['Carga Q (N)'])
plt.grid()

plt.figure(4)
plt.plot(tempo, e_a_lista, 'r')
plt.title("Entrada - FEM Senoidal do Motor (V)")
plt.legend(['Tensão (V)'])
plt.grid()

#Análises
# Estabilidade - Polos do Sistema 
A = np.array([[0, 0, 0, 1, 0, 0, 0],
              [0, 0, 0, 0, 1, 0, 0],
              [0, 0, 0, 0, 0, 1, 0],
              [-k_m/I_m, k_m/I_m, 0, -c_m/I_m, c_m/I_m, 0, K/I_m],
              [(k_m*n**2)/(I_1*n**2 + I_2), -(k_m*n**2 + k_L)/(I_1*n**2 + I_2), (k_L*n)/(I_1*n**2 + I_2), (c_m*n**2)/(I_1*n**2 + I_2), -(c_m*n**2 + c_L)/(I_1*n**2 + I_2), (c_L*n)/(I_1*n**2 + I_2), 0],
              [0, k_L/(n*I_L), -k_L/I_L, 0, c_L/(n*I_L), -c_L/I_L, 0],
              [0, 0, 0, -K/L_a, 0, 0, -R_a/L_a]])

# plot dos autovalores no plano complexo
plt.figure(5, figsize=(10, 6))
plt.plot(np.real(np.linalg.eig(A)[0]), np.imag(np.linalg.eig(A)[0]), 'x')
plt.title('Posição dos polos no plano complexo')
plt.xscale('linear')
plt.xlabel('Eixo Real')
plt.ylabel('Eixo Imaginário')
plt.gca().set_facecolor('xkcd:light grey')
plt.grid()
plt.show()

B = np.array([
      [0],
      [0],
      [0],
      [0],
      [0],
      [1/I_L],
      [0]
  ])

C = np.array([0, 0, 0, 0, 0, 1, 0])  # Só funciona para sistemas SISO

D = np.array([0])


# Resposta do sistema linear no domínio da frequência
# Obtenção dos diagramas de Bode a partir do Modelo Linearizado
#OBS.: ALTERAR A MATRIZ C DEPENDENDO DO PARÂMETRO ANALISADO - SISTEMA DEVE SER SISO
sys = scp.signal.lti(A, B, C, D)
w, mag, phase = scp.signal.bode(sys, w=np.logspace(1e-2,1e4,10000))


# calculate the transfer function
print(scp.signal.ss2tf(A, B, C, D))
