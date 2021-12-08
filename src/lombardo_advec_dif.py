import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import chaospy as cp
import math
sns.set()

T_final = 1 #
h_t = 0.01

L = 100  # 
h_x = 1

# Parametros
chi = 15  # [4, 55]
tau = 1  # [0.001, 1]
epsilon = 0.5  # [0.5, 1.5]
beta = 1  # [0.2, 1]
delta = 1  # [0, 1]
r = 6 # [1, 6]

d_dc = 1 # difusao DC convencional(procurar na literatura)
d_da = 1 # difusao DC ativada(procurar na literatura)
d_t_cit = 1 # difusao t citotóxica(procurar na literatura)
d_anti = 1 # difusao anticorpo(procurar na literatura)
lamb_f_dc = 0.3 # taxa de anticorpos consumidos durante o processo de opsonização por DC convencional
lamb_f_da = 0.3 # taxa de anticorpos consumidos durante o processo de opsonização por DC ativada
lamb_f_m = 0.3 # taxa de anticorpos consumidos durante o processo de opsonização pela micróglia
b_d = 1 # taxa de ativacao de dc por odc(procurar na literatura)
r_dc = 1 # agressividade de dc convencional(procurar na literatura)
r_da = 9 # agressividade de dc ativada(procurar na literatura)
r_t = 1  # agressividade de t citotoxica(procurar na literatura)


# IC
# Macrofagos
mac_atual = np.zeros((int(L/h_x), int(L/h_x)))
# s = np.random.normal(0,1,(int(L/h_x), int(L/h_x)))
# x, y = np.random.multivariate_normal(0, 1, 5000).T
# plt.plot(x, y, 'x')
# plt.axis('equal')
# plt.show()
# distribution = cp.Normal(0,1)
# uloc = np.linspace(0, 1, int(L/h_x))
# xloc = distribution.inv(uloc)
# s = distribution.pdf(xloc).round(1)
#print(s)
# print(xloc.round(int(L/(h_x*2))))
for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if (i-int(L/h_x)/2)**2 + (j-int(L/h_x)/2)**2 < 5:
            mac_atual[i][j] = 0.2 # 0.1*math.exp(-1*( (i-int(L/h_x)/2)**2 + (j-int(L/h_x)/2)**2)/s[i,j])

# Dendríticas convencionais
dendritica_conv_atual = np.ones((int(L/h_x), int(L/h_x)))
dendritica_conv_anterior = np.ones((int(L/h_x), int(L/h_x)))

# Dendríticas ativadas
dendritica_ativ_atual = np.zeros((int(L/h_x), int(L/h_x)))
dendritica_ativ_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# Macrofagos
mac_anterior = np.copy(mac_atual)

# Citocinas pro-inflamatorias
cit_atual = np.zeros((int(L/h_x), int(L/h_x)))
cit_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# Ol destruidos
olide_atual = np.zeros((int(L/h_x), int(L/h_x)))
olide_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# T citotóxica
t_cito_atual = 0.1*np.ones((int(L/h_x), int(L/h_x)))
t_cito_anterior = 0.1*np.ones((int(L/h_x), int(L/h_x)))

# anticorpo
anticorpo_atual = 0.1*np.zeros((int(L/h_x), int(L/h_x)))
anticorpo_anterior = 0.1*np.zeros((int(L/h_x), int(L/h_x)))

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

for k in range(steps):
    for i in range(tam):
        for j in range(tam):
            d = olide_anterior[i][j]
            c = cit_anterior[i][j]
            m = mac_anterior[i][j]
            dc = dendritica_conv_anterior[i][j]
            da = dendritica_ativ_anterior[i][j]
            f = anticorpo_anterior[i][j]
            t_cito = t_cito_anterior[i][j]

            # condição de contorno de Neumman macrofagos
            mac_ipj = mac_anterior[i+1][j] if i != tam-1 else m - 2*h_x*u_c
            mac_imj = mac_anterior[i-1][j] if i != 0 else m - 2*h_x*u_b
            mac_ijp = mac_anterior[i][j+1] if j != tam-1 else m - 2*h_x*u_c
            mac_ijm = mac_anterior[i][j-1] if j != 0 else m - 2*h_x*u_d

            # condição de contorno de Neumman dc convencional
            dc_ipj = dendritica_conv_anterior[i+1][j] if i != tam-1 else dc - 2*h_x*u_c
            dc_imj = dendritica_conv_anterior[i-1][j] if i != 0 else dc - 2*h_x*u_b
            dc_ijp = dendritica_conv_anterior[i][j+1] if j != tam-1 else dc - 2*h_x*u_c
            dc_ijm = dendritica_conv_anterior[i][j-1] if j != 0 else dc - 2*h_x*u_d

            # condição de contorno de Neumman de ativadas
            da_ipj = dendritica_ativ_anterior[i+1][j] if i != tam-1 else da - 2*h_x*u_c
            da_imj = dendritica_ativ_anterior[i-1][j] if i != 0 else da - 2*h_x*u_b
            da_ijp = dendritica_ativ_anterior[i][j+1] if j != tam-1 else da - 2*h_x*u_c
            da_ijm = dendritica_ativ_anterior[i][j-1] if j != 0 else da - 2*h_x*u_d
            
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
            
            # condição de contorno de Neumman t citotóxicas
            t_cito_ipj = t_cito_anterior[i+1][j] if i != tam-1 else t_cito - 2*h_x*u_c
            t_cito_imj = t_cito_anterior[i-1][j] if i != 0 else t_cito - 2*h_x*u_b
            t_cito_ijp = t_cito_anterior[i][j+1] if j != tam-1 else t_cito - 2*h_x*u_c
            t_cito_ijm = t_cito_anterior[i][j-1] if j != 0 else t_cito - 2*h_x*u_d

            # condição de contorno de Neumman anticorpos
            f_ipj = anticorpo_anterior[i+1][j] if i != tam-1 else f - 2*h_x*u_c
            f_imj = anticorpo_anterior[i-1][j] if i != 0 else f - 2*h_x*u_b
            f_ijp = anticorpo_anterior[i][j+1] if j != tam-1 else f - 2*h_x*u_c
            f_ijm = anticorpo_anterior[i][j-1] if j != 0 else f - 2*h_x*u_d
            
            #Dados da equacao macrofagos

            #Dependendo do gradiente das citocinas vou fazer upwind ou downwind no eixo i ou eixo j

            #Decidindo qual combinacao usar no gradiente dos macrofagos
            
            gradiente_c_i = (cit_ipj - cit_imj)/(2*h_x)
            gradiente_c_j = (cit_ijp - cit_ijm)/(2*h_x)

            if gradiente_c_i > 0:
                gradiente_m_i = (m/(1 + m) - mac_imj/(1 + mac_imj))/h_x
                gradiente_dc_i = (dc/(1+dc) -dc_imj/(1+dc_imj))/h_x
                gradiente_da_i = (da/(1+da) -da_imj/(1+da_imj))/h_x
                gradiente_t_i = (t_cito/(1+t_cito) -t_cito_imj/(1+t_cito_imj))/h_x
            else:
                gradiente_m_i = (mac_ipj/(1 + mac_ipj) - m/(1 + m))/h_x
                gradiente_dc_i = (dc_ipj/(1 + dc_ipj) - dc/(1 + dc))/h_x
                gradiente_da_i = (da_ipj/(1 + da_ipj) - da/(1 + da))/h_x
                gradiente_t_i = (t_cito_ipj/(1 + t_cito_ipj) - t_cito/(1 + t_cito))/h_x
            if gradiente_c_j > 0:
                gradiente_m_j = (m/(1 + m) - mac_ijm/(1 + mac_ijm))/h_x
                gradiente_dc_j = (dc/(1 + dc) - dc_ijm/(1 + dc_ijm))/h_x
                gradiente_da_j = (da/(1 + da) - da_ijm/(1 + da_ijm))/h_x
                gradiente_t_j = (t_cito/(1 + t_cito) - t_cito_ijm/(1 + t_cito_ijm))/h_x
            else:
                gradiente_m_j = (mac_ijp/(1 + mac_ijp) - m/(1 + m))/h_x
                gradiente_dc_j = (dc_ijp/(1 + dc_ijp) - dc/(1 + dc))/h_x
                gradiente_da_j = (da_ijp/(1 + da_ijp) - da/(1 + da))/h_x
                gradiente_t_j = (t_cito_ijp/(1 + t_cito_ijp) - t_cito/(1 + t_cito))/h_x

            quimiotaxia_mac = chi*(gradiente_c_i*gradiente_m_i + gradiente_c_j*gradiente_m_j)
            difusao_mac = (mac_ipj + mac_imj - 4*m + mac_ijp + mac_ijm )/h_x**2
            reacao_mac = m*(1 - m)
            
            mac_atual[i][j] = m + h_t*(difusao_mac + reacao_mac - quimiotaxia_mac)
            
            #DC convencional
            quimiotaxia_dc = chi*(gradiente_c_i*gradiente_dc_i + gradiente_c_j*gradiente_dc_j)
            difusao_dc = d_dc*(dc_ipj + dc_imj - 4*dc + dc_ijp + dc_ijm )/h_x**2

            dendritica_conv_atual[i][j] = dc + h_t*(difusao_dc - quimiotaxia_dc - b_d*(1-d)*dc)

            #DC ativada
            quimiotaxia_da = chi*(gradiente_c_i*gradiente_da_i + gradiente_c_j*gradiente_da_j)
            difusao_da = d_da*(da_ipj + da_imj - 4*da + da_ijp + da_ijm )/h_x**2

            dendritica_ativ_atual[i][j] = da + h_t*(difusao_da - quimiotaxia_da + b_d*(1-d)*dc)

            #T citotóxica
            quimiotaxia_t_cito = chi*(gradiente_c_i*gradiente_t_i + gradiente_c_j*gradiente_t_j)
            difusao_t_cito = d_t_cit*(t_cito_ijm + t_cito_ijp - 4*t_cito + t_cito_imj + t_cito_ipj)/h_x**2
            t_cito_atual[i][j] = t_cito + h_t*(difusao_t_cito - quimiotaxia_t_cito) #termo de migração faltando!!

            # anticorpo
            difusao_anticorpo = d_anti*(f_ipj + f_imj - 4*f + f_ijp + f_ijm)
            anticorpo_atual[i][j] = f + h_t*(difusao_anticorpo - lamb_f_da*f*dc - lamb_f_da*f*da - lamb_f_m*f*m)#termo de migração faltando!!

            #Dados da equacao citocinas
            difusao_cit = epsilon*(cit_ipj + cit_imj - 4*c + cit_ijp + cit_ijm )/h_x**2
            reacao_cit = delta*d - c + beta*m

            cit_atual[i][j] = c + h_t*(difusao_cit + reacao_cit)/tau

            #Dados da equacao oligodendrocitos destruidos
            olide_atual[i][j] = d + h_t*(r*(m/(1 + m))*m*(1 - d) + r_dc*(dc/(1 + dc))*dc*(1 - d) + r_da*(da/(1 + da))*da*(1 - d) + r_t*(t_cito/(1 + t_cito))*t_cito*(1 - d))
            

    mac_anterior = np.copy(mac_atual)
    cit_anterior = np.copy(cit_atual)
    olide_anterior = np.copy(olide_atual)
    dendritica_conv_anterior = np.copy(dendritica_conv_atual)
    dendritica_ativ_anterior = np.copy(dendritica_ativ_atual)
    t_cito_anterior = np.copy(t_cito_atual)
    anticorpo_anterior = np.copy(anticorpo_atual)

    if k%p ==0 or k == steps-1:
            x_pts, y_pts = np.meshgrid(x, x)
            
            #results odc
            cp = plt.contourf(x_pts, y_pts,olide_anterior, 100)
            plt.title("Tempo: "+"{:.4f}".format(k*h_t))
            plt.colorbar(cp, label='ODC-destruidos')
            plt.contourf(x_pts, y_pts, olide_anterior,100)
            plt.savefig('../results/oligodendrocitos/fig'+"{:.4f}".format(k*h_t)+'.png', dpi = 300)
            plt.clf()
            
            #results macrofagos
            cp = plt.contourf(x_pts, y_pts,mac_anterior, 100)
            plt.title("Tempo: "+"{:.4f}".format(k*h_t))
            plt.colorbar(cp, label='Macrófagos')
            plt.contourf(x_pts, y_pts, mac_anterior,100)
            plt.savefig('../results/macrofagos/fig'+"{:.4f}".format(k*h_t)+'.png', dpi = 300)
            plt.clf()
            
            #results dc convencional
            cp = plt.contourf(x_pts, y_pts,dendritica_conv_anterior, 100)
            plt.title("Tempo: "+"{:.4f}".format(k*h_t))
            plt.colorbar(cp, label='DC- convencional')
            plt.contourf(x_pts, y_pts, dendritica_conv_anterior,100)
            plt.savefig('../results/dendriticas_convencionais/fig'+"{:.4f}".format(k*h_t)+'.png', dpi = 300)
            plt.clf()

            #results dc ativada
            cp = plt.contourf(x_pts, y_pts,dendritica_ativ_anterior, 100)
            plt.title("Tempo: "+"{:.4f}".format(k*h_t))
            plt.colorbar(cp, label='DC- ativada')
            plt.contourf(x_pts, y_pts, dendritica_ativ_anterior,100)
            plt.savefig('../results/dendriticas_ativadas/fig'+"{:.4f}".format(k*h_t)+'.png', dpi = 300)
            plt.clf()

            #results citocina
            cp = plt.contourf(x_pts, y_pts,cit_anterior, 100)
            plt.title("Tempo: "+"{:.4f}".format(k*h_t))
            plt.colorbar(cp, label='Citocina')
            plt.contourf(x_pts, y_pts, cit_anterior,100)
            plt.savefig('../results/citocinas/fig'+"{:.4f}".format(k*h_t)+'.png', dpi = 300)
            plt.clf()

            #results T citotóxica
            cp = plt.contourf(x_pts, y_pts,t_cito_anterior, 100)
            plt.title("Tempo: "+"{:.4f}".format(k*h_t))
            plt.colorbar(cp, label='T citotóxica')
            plt.contourf(x_pts, y_pts, t_cito_anterior,100)
            plt.savefig('../results/t_citotoxica/fig'+"{:.4f}".format(k*h_t)+'.png', dpi = 300)
            plt.clf()

            #results Anticorpo
            cp = plt.contourf(x_pts, y_pts,anticorpo_anterior, 100)
            plt.title("Tempo: "+"{:.4f}".format(k*h_t))
            plt.colorbar(cp, label='Anticorpo')
            plt.contourf(x_pts, y_pts, anticorpo_anterior,100)
            plt.savefig('../results/anticorpos/fig'+"{:.4f}".format(k*h_t)+'.png', dpi = 300)
            plt.clf()
    print('tempo: '+str(k*h_t))