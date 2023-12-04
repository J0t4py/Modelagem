import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy import integrate
from scipy.spatial.transform import Rotation
!pip install control
import control as ct
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

def fun(t, y):
  # Definição do sinal de entrada
  #Entradas
  Q = 0 # Esforço aplicado à torreta - entrada em degrau unitário
  e_a = 0 #Força eletromotriz do motor


  y_dot = np.array([y[3], y[4], y[5], (K*y[6] - k_m*(y[0]-y[1]) - c_m*(y[3] - y[4]))/I_m,
                    (n**2/(n**2*I_1 + I_2))*(k_m*y[0] + c_m*y[3] + (k_L/n)*y[2] + (c_L/n)*y[5] - y[4]*(c_m + c_L/n**2) - y[1]*(k_m + k_L/n**2)),
                    (1/I_L)*(Q - c_L*y[5] - k_L*y[2] + (k_L/n)*y[1] + (c_L/n)*y[4]), (1/L_a)*(e_a - K*y[3] - R_a*y[6])])

  # velocidades angulares do motor, pinhão e carga
  theta_dot = np.array([y_dot[0], y_dot[1], y_dot[2]])
  # acelerações angulares do motor, pinhão e carga
  theta2_dot = np.array([y_dot[3], y_dot[4], y_dot[5]])
  # derivada temporal da corrente
  i_p = np.array(y_dot[6])
  #return theta_dot, theta2_dot, i

  return y_dot

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

# resolução do sistema de EDO's pelo método de Adams (LSODA)
sol = solve_ivp(fun, (0,10), y_0, method='LSODA', t_eval=np.linspace(0, 10, 100))

Q_deg_lista = ([]) # Esforço aplicado à torreta - entrada em degrau
e_a_lista = ([])
#calcula as entradas para cada instante de tempo
for t in sol.t:
  e_a_lista.append(75*np.sin(0.03*np.pi*t))
  if 20 < t < 50:
    Q_deg_lista.append(4)
  elif 75 < t < 85:
    Q_deg_lista.append(-4)
  else:
    Q_deg_lista.append(0)

t_eval=np.linspace(0, 100, 200)

plt.figure(0)
plt.plot(sol.t, sol.y[0, :])
plt.plot(sol.t, sol.y[1, :])
plt.plot(sol.t, sol.y[2, :])

plt.legend(['theta_m', 'theta_1','theta_L'])

plt.grid()

plt.figure(1)

plt.plot(sol.t, sol.y[3, :])
plt.plot(sol.t, sol.y[4, :])
plt.plot(sol.t, sol.y[5, :])
plt.plot(sol.t, sol.y[6, :])

plt.legend(['theta_p_m', 'theta_p_1', 'theta_p_L', 'i_a'])

plt.grid()
plt.show
