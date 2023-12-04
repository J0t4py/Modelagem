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
def fun2(t, y,u1,u2):
  # Definição do sinal de entrada
  #Entradas
  Q = 3.8 
  e_a = 32/(np.sqrt(2))
  # Q = 0
  # e_a=0

  #e_a = 75*np.sin(0.03*np.pi*time_vector)
  #e_a =  #fem do motor
  #Q =  # Esforço aplicado à torreta - entrada em degrau unitário

  y_dot = np.array([y[3], y[4], y[5], (K*y[6] - k_m*(y[0]-y[1]) - c_m*(y[3] - y[4]))/I_m,
                    (n**2/(n**2*I_1 + I_2))*(k_m*y[0] + c_m*y[3] + (k_L/n)*y[2] + (c_L/n)*y[5] - y[4]*(c_m + c_L/n**2) - y[1]*(k_m + k_L/n**2)),
                    (1/I_L)*(Q - c_L*y[5] - k_L*y[2] + (k_L/n)*y[1] + (c_L/n)*y[4]), (1/L_a)*(e_a - K*y[3] - R_a*y[6])])

  #y_dot = y_dot.reshape((8,1))
  # velocidades angulares do motor, pinhão e carga
  theta_dot = np.array([y_dot[0], y_dot[1], y_dot[2]])
  # acelerações angulares do motor, pinhão e carga
  theta2_dot = np.array([y_dot[3], y_dot[4], y_dot[5]])
  # derivada temporal da corrente
  i_p = np.array(y_dot[6])
  #return theta_dot, theta2_dot, i

  return y_dot
time_vector = np.linspace(start=0, stop=10, num=100)
U = np.zeros((1, 100))
def output(t, y,u1,params):
  return y[0],y[1],y[2],y[3],y[4],y[5],y[6]
model = ct.NonlinearIOSystem(fun2, outfcn=output, inputs=['dy'], output=('y','y2','y3','y4','y5','y6','y7'))
ss = ct.linearize(sys=model, xeq=[0, 0, 0,0,0,0,0], ueq=[0])
X0 = [0.02, 0.03, 0.04, 0, 0, 0, 0]
t_linear, y_linear = ct.input_output_response(ss, time_vector, U, X0, method='LSODA')
fig, (ax1) = plt.subplots(1, 1, figsize=(8, 6))
ax1.plot(t_linear, y_linear[0], label='θ_m')
ax1.plot(t_linear, y_linear[1], label='θ_1')
ax1.plot(t_linear, y_linear[2], label='θ_l')
ax1.set_title("Posições angulares linearizadas")
ax1.legend()
ax1.grid()
fig, (ax2) = plt.subplots(1, 1, figsize=(8, 6))
ax2.plot(t_linear, y_linear[3], label='i_a')
ax2.plot(t_linear, y_linear[4], label='ω_m')
ax2.plot(t_linear, y_linear[5], label='ω_1')
ax2.plot(t_linear, y_linear[6], label='ω_l')
ax2.set_title("Corrente e velocidades angulares linearizadas")
ax2.legend()
ax2.grid()
plt.tight_layout()
plt.show()
