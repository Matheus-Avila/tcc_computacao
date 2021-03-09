import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
sns.set()

T_final = 7 #
h_t = 0.001

L = 100  # 
h_x = 1

# Parametros
chi = 15  # [4, 55]
tau = 1  # [0.001, 1]
epsilon = 0.5  # [0.5, 1.5]
beta = 1  # [0.2, 1]
delta = 1  # [0, 1]
r = 6 # [1, 6]

# IC
# Macrofagos
mac_atual = np.zeros((int(L/h_x), int(L/h_x)))


for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if (i-int(L/h_x)/2)**2 + (j-int(L/h_x)/2)**2 < 20:
            mac_atual[i][j] = 1

# mac_atual[int(0.45*L/h_x):int(0.55*L/h_x)+1, int(0.45*L/h_x):int(0.55*L/h_x)+1] = 1
mac_anterior = np.copy(mac_atual)
# Citocinas pro-inflamatorias
cit_atual = np.zeros((int(L/h_x), int(L/h_x)))
cit_anterior = np.zeros((int(L/h_x), int(L/h_x)))
# Ol destruidos
olide_atual = np.zeros((int(L/h_x), int(L/h_x)))
olide_anterior = np.zeros((int(L/h_x), int(L/h_x)))

#BC
u_a = 0#Neumann cima
u_b = 0#Neumann direita
u_c = 0#Neumann baixo
u_d = 0#Neumann esquera


t = np.linspace(0, T_final, int(T_final/h_t))
x = np.linspace(0, L, int(L/h_x))
tam = len(x)
steps = len(t)

p = int(steps/100)

sol_tempo_d = []
sol_tempo_d.append(olide_anterior)
sol_tempo_m = []
sol_tempo_m.append(mac_anterior)

#Funcoes para escolha de Chi(m)

for k in range(steps):
    for i in range(tam):
        for j in range(tam):
            d = olide_anterior[i][j]
            c = cit_anterior[i][j]
            m = mac_anterior[i][j]
            # condição de contorno de Neumman macrofagos
            mac_ipj = mac_anterior[i+1][j] if i != tam-1 else m - 2*h_x*u_c
            mac_imj = mac_anterior[i-1][j] if i != 0 else m - 2*h_x*u_b
            mac_ijp = mac_anterior[i][j+1] if j != tam-1 else m - 2*h_x*u_c
            mac_ijm = mac_anterior[i][j-1] if j != 0 else m - 2*h_x*u_d
            
            # condição de contorno de Neumman oligodendrocitos destruidos
            olide_ipj = olide_anterior[i+1][j] if i != tam-1 else d - 2*h_x*u_c
            olide_imj = olide_anterior[i-1][j] if i != 0 else d - 2*h_x*u_b
            olide_ijp = olide_anterior[i][j+1] if j != tam-1 else d - 2*h_x*u_c
            olide_ijm = olide_anterior[i][j-1] if j != 0 else d - 2*h_x*u_d

            # condição de contorno de Neumman citocinas
            cit_ipj = cit_anterior[i+1][j] if i != tam-1 else c - 2*h_x*u_c
            cit_imj = cit_anterior[i-1][j] if i != 0 else c - 2*h_x*u_b
            cit_ijp = cit_anterior[i][j+1] if j != tam-1 else c - 2*h_x*u_c
            cit_ijm = cit_anterior[i][j-1] if j != 0 else c - 2*h_x*u_d
            
            #Dados da equacao macrofagos

            #Dependendo do quadrante vou fazer upwind ou downwind no eixo i ou eixo j

            #Decidindo qual combinacao usar no gradiente dos macrofagos
            
            gradiente_c_i = (cit_ipj - cit_imj)/(2*h_x)
            gradiente_c_j = (cit_ijp - cit_ijm)/(2*h_x)

            if gradiente_c_i > 0:
                up_wind_i = (m/(1 + m) - \
                mac_imj/(1 + mac_imj))/h_x
                gradiente_m_i = up_wind_i
            else:
                down_wind_i = (mac_ipj/(1 + mac_ipj) - \
                m/(1 + m))/h_x
                gradiente_m_i = down_wind_i
            if gradiente_c_j > 0:
                up_wind_j = (m/(1 + m) - \
                mac_ijm/(1 + mac_ijm))/h_x
                gradiente_m_j = up_wind_j
            else:
                down_wind_j = (mac_ijp/(1 + mac_ijp) - \
                m/(1 + m))/h_x
                gradiente_m_j = down_wind_j
            

            quimiotaxia_mac = chi*(gradiente_c_i*gradiente_m_i + gradiente_c_j*gradiente_m_j)
            difusao_mac = (mac_ipj + mac_imj - 4*m + mac_ijp + mac_ijm )/h_x**2
            reacao_mac = m*(1 - m)
            
            mac_atual[i][j] = m + h_t*(difusao_mac + reacao_mac - quimiotaxia_mac)

            #Dados da equacao citocinas
            difusao_cit = epsilon*(cit_ipj + cit_imj - 4*c + cit_ijp + cit_ijm )/h_x**2
            reacao_cit = delta*d - c + beta*m

            cit_atual[i][j] = c + h_t*(difusao_cit + reacao_cit)/tau

            #Dados da equacao oligodendrocitos destruidos
            olide_atual[i][j] = d + h_t*r*(m/(1 + m))*m*(1 - d)
            

    mac_anterior = np.copy(mac_atual)
    cit_anterior = np.copy(cit_atual)
    olide_anterior = np.copy(olide_atual)
    if k%p ==0 or k == steps-1:
            x_pts, y_pts = np.meshgrid(x, x)
            cp = plt.contourf(x_pts, y_pts,olide_anterior, 100)
            plt.title("Tempo: "+"{:.4f}".format(k*h_t))
            plt.colorbar(cp, label='Ols-destruidos')
            plt.contourf(x_pts, y_pts, olide_anterior,100)
            plt.savefig('figs/fig'+"{:.4f}".format(k*h_t)+'.png', dpi = 300)
            plt.clf()
    print('tempo: '+str(k*h_t))