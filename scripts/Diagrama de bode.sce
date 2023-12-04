s = poly(0,"s")

// DEFINIÇÃO DAS CONSTANTES 

I_m = 3*(10^-4) //#Inércia do motor em kg.m^2
I_1 = 3*(10^-4) //#Inércia do motor em kg.m^2
I_2 = 1.5 //#Inércia da engrenagem em kg.m^2
I_L = 50 //#Inércia da carga em kg.m^2
k_L = 5*(10^7) //#Rigidez torsional da carga em Nm/rad
k_m = 1.2*(10^3) //#Rigidez torsional do motor em Nm/rad
c_L = 10 //#coeficiente de anmortecimento do rolamento do munhão em Nms/rad
c_m = 0.3 //#coeficiente de amortecimento do motor
K = 60 //#constante de torque do motor em Nm/Amp  #TESTAR COM 60
R_a = 0.8 //#Resistência elétrica do motor em Ohms
L_a = 1*(10^-3) //#Indutância do motor em henry
Np = 14 //#número de dentes do pinhão
Ne = 408 //#número de dentes do munhão
n = Ne/Np //e/Np # Relação de transmissão

// DEFINIÇÃO DAS MATRIZES 

A = [0, 0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 0, 1, 0;
    -k_m/I_m, k_m/I_m, 0, -c_m/I_m, c_m/I_m, 0, K/I_m;
    (k_m*n^2)/(I_1*n^2 + I_2), -(k_m*n^2 + k_L)/(I_1*n^2 + I_2), (k_L*n)/(I_1*n^2 + I_2), (c_m*n^2)/(I_1*n^2 + I_2), -(c_m*n^2 + c_L)/(I_1*n^2 + I_2), (c_L*n)/(I_1*n^2 + I_2), 0;
    0, k_L/(n*I_L), -k_L/I_L, 0, c_L/(n*I_L), -c_L/I_L, 0;
    0, 0, 0, -K/L_a, 0, 0, -R_a/L_a]
    
B1 = [0;
     0;
     0;
     0;
     0;
     (1/I_L);
     0]

B2 = [0;
     0;
     0;
     0;
     0;
     0;
    (1/L_a)]

D = [0]

Css1 = [0, 0, 1, 0, 0, 0, 0]
Css2 = [0, 0, 0, 0, 0, 1, 0]

//CÁLCULO DOS POLOS 

Phi = (s*eye(7,7)-A)
dem = detr(Phi)

// FUNÇÕES DE TRANFERENCIA

mat11 = [Phi, -B1; Css1, D]
num11 = detr(mat11)
G11 = num11/dem 
disp("Função de Transferência G11:",G11)
G11 = syslin('c',G11)

mat12 = [Phi, -B2; Css1, D]
num12 = detr(mat12)
G12 = num12/dem 
disp("Função de Transferência G12:",G12)
G12 = syslin('c',G12)

mat21 = [Phi, -B1; Css2, D]
num21 = detr(mat21)
G21 = num21/dem 
disp("Função de Transferência G21:",G21)
G21 = syslin('c',G21)

mat22 = [Phi, -B2; Css2, D]
num22 = detr(mat22)
G22 = num22/dem 
disp("Função de Transferência G22:",G22)
G22 = syslin('c',G22)

//DIAGRAMAS DE BODE

scf(1)
bode(G21)
xtitle("Diagrama de Bode - Entrada Q e Saída dot(θl) ")

scf(2)
bode(G22)
xtitle("Diagrama de Bode - Entrada e_a e Saída dot(θl) ")

scf(3)
bode(G11)
xtitle("Diagrama de Bode - Entrada Q e Saída θl")


scf(4)
bode(G12)
xtitle("Diagrama de Bode - Entrada e_a e Saída θl")

// POLOS DO SISTEMA 
disp("Equação Característica:",dem)
disp("Polos do Sistema:",roots(dem))

scf(5)
plzr(G11);
xtitle("Posição dos Polos no Plano Complexo")
xlabel("Eixo Real")
ylabel("Eixo Imaginário")
a = gca()
a.box="on"
a.data_bounds=[-1200, -120000; 1200, 120000]


//ANÁLISE DO CRITÉRIO DE ROUTH


[r , num] = routh_t(dem)
disp("Coeficientes Tabela de Rough:",r)
disp("Número de Troca de Sinais:",num)
if num==0
   disp("System is stable")
else
   mprintf("There is %g sign changes in entries of first column.\nTherefore, system is unstable.", num)
end




