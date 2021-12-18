import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import chaospy as cp
sns.set()

gradiente = lambda ponto_anterior, ponto_posterior: quimiotaxia(ponto_posterior) - quimiotaxia(ponto_anterior)
quimiotaxia = lambda ponto_atual: ponto_atual/(1+ponto_atual)
f_func = lambda populacao: populacao*populacao/(1+populacao)

T_final = 7 #
h_t = 0.01

L = 100  # 
h_x = 1

# Parametros
chi = 15  # [4, 55]
tau = 1  # [0.001, 1]
epsilon = 0.5  # [0.5, 1.5]
r_m = 6 # [1, 6]

d_dc = 1 # difusao DC convencional(procurar na literatura)
d_da = 1 # difusao DC ativada(procurar na literatura)
d_t_cit = 1 # difusao t citotóxica(procurar na literatura)
d_anti = 1 # difusao anticorpo(procurar na literatura)
lamb_f_dc = 0.3 # taxa de anticorpos consumidos durante o processo de opsonização por DC convencional
lamb_f_da = 0.3 # taxa de anticorpos consumidos durante o processo de opsonização por DC ativada
lamb_f_m = 0.3 # taxa de anticorpos consumidos durante o processo de opsonização pela micróglia
b_d = 1 # taxa de ativacao de dc por odc(procurar na literatura)
r_dc = 1 # agressividade de dc convencional(procurar na literatura)
r_t = 1  # agressividade de t citotoxica(procurar na literatura)

d_c_homeostase = 0.5 # valor de homeostase de células dendríticas
mu = 0.3 #Valor de homeostase de células dendríticas
alpha_d = 0.4 #Taxa de migração de DC ativadas para o linfonodo
gamma_at = 0.2 #Taxa de migração de anticorpos para o tecido
gamma_tcito = 0.2 #Taxa de migração de T citotoxica para o tecido
da_linf = 0.4 #DC ativadas no linfonodo
at_linf = 0.1 #Anticorpos no linfonodo 
t_cito_linf = 0.1 #T citotoxica no linfonodo 

# IC
# Macrofagos
# Macrofagos
mac_anterior = np.zeros((int(L/h_x), int(L/h_x)))
for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if (i-int(L/h_x)/2)**2 + (j-int(L/h_x)/2)**2 < 5:
            mac_anterior[i][j] = 0.2

# T citotóxica
t_cito_anterior = 0.3*t_cito_linf*np.ones((int(L/h_x), int(L/h_x)))

# Ol destruidos
olide_anterior = 0*np.ones((int(L/h_x), int(L/h_x)))

# anticorpo
anticorpo_anterior = 0.3*at_linf*np.ones((int(L/h_x), int(L/h_x)))

# Dendríticas convencionais
dendritica_conv_anterior = 0.8*d_c_homeostase*np.ones((int(L/h_x), int(L/h_x)))

# Dendríticas ativadas
dendritica_ativ_anterior = 0.1*da_linf*np.ones((int(L/h_x), int(L/h_x)))

mac_atual = np.zeros((int(L/h_x), int(L/h_x)))
t_cito_atual = np.zeros((int(L/h_x), int(L/h_x)))
olide_atual = np.zeros((int(L/h_x), int(L/h_x)))
anticorpo_atual = np.zeros((int(L/h_x), int(L/h_x)))
dendritica_conv_atual = np.zeros((int(L/h_x), int(L/h_x)))
dendritica_ativ_atual = np.zeros((int(L/h_x), int(L/h_x)))

#BC
bc_neumann_cima = 0#Neumann cima
bc_neumann_direita = 0#Neumann direita
bc_neumann_baixo = 0#Neumann baixo
bc_neumann_esquerda = 0#Neumann esquerda


t = np.linspace(0, T_final, int(T_final/h_t))
x = np.linspace(0, L, int(L/h_x))
tam = len(x)
steps = len(t)

intervalo_figs = int(steps/100)

for k in range(steps):
    for i in range(tam):
        for j in range(tam):
            oligo_destr = olide_anterior[i][j]
            m = mac_anterior[i][j]
            dc = dendritica_conv_anterior[i][j]
            da = dendritica_ativ_anterior[i][j]
            f = anticorpo_anterior[i][j]
            t_cito = t_cito_anterior[i][j]

            # condição de contorno de Neumman macrofagos
            mac_ipj = mac_anterior[i+1][j] if i != tam-1 else m - 2*h_x*bc_neumann_baixo
            mac_imj = mac_anterior[i-1][j] if i != 0 else m - 2*h_x*bc_neumann_cima
            mac_ijp = mac_anterior[i][j+1] if j != tam-1 else m - 2*h_x*bc_neumann_direita
            mac_ijm = mac_anterior[i][j-1] if j != 0 else m - 2*h_x*bc_neumann_esquerda

            # condição de contorno de Neumman dc convencional
            dc_ipj = dendritica_conv_anterior[i+1][j] if i != tam-1 else dc - 2*h_x*bc_neumann_baixo
            dc_imj = dendritica_conv_anterior[i-1][j] if i != 0 else dc - 2*h_x*bc_neumann_cima
            dc_ijp = dendritica_conv_anterior[i][j+1] if j != tam-1 else dc - 2*h_x*bc_neumann_direita
            dc_ijm = dendritica_conv_anterior[i][j-1] if j != 0 else dc - 2*h_x*bc_neumann_esquerda

            # condição de contorno de Neumman de ativadas
            da_ipj = dendritica_ativ_anterior[i+1][j] if i != tam-1 else da - 2*h_x*bc_neumann_baixo
            da_imj = dendritica_ativ_anterior[i-1][j] if i != 0 else da - 2*h_x*bc_neumann_cima
            da_ijp = dendritica_ativ_anterior[i][j+1] if j != tam-1 else da - 2*h_x*bc_neumann_direita
            da_ijm = dendritica_ativ_anterior[i][j-1] if j != 0 else da - 2*h_x*bc_neumann_esquerda
            
            # condição de contorno de Neumman oligodendrocitos destruidos ODC é EDO!!! Nao tem BC!!!! Desse jeito ta certo?
            olide_ipj = olide_anterior[i+1][j] if i != tam-1 else oligo_destr
            olide_imj = olide_anterior[i-1][j] if i != 0 else oligo_destr
            olide_ijp = olide_anterior[i][j+1] if j != tam-1 else oligo_destr
            olide_ijm = olide_anterior[i][j-1] if j != 0 else oligo_destr

            # condição de contorno de Neumman t citotóxicas
            t_cito_ipj = t_cito_anterior[i+1][j] if i != tam-1 else t_cito - 2*h_x*bc_neumann_baixo
            t_cito_imj = t_cito_anterior[i-1][j] if i != 0 else t_cito - 2*h_x*bc_neumann_cima
            t_cito_ijp = t_cito_anterior[i][j+1] if j != tam-1 else t_cito - 2*h_x*bc_neumann_direita
            t_cito_ijm = t_cito_anterior[i][j-1] if j != 0 else t_cito - 2*h_x*bc_neumann_esquerda

            # condição de contorno de Neumman anticorpos
            f_ipj = anticorpo_anterior[i+1][j] if i != tam-1 else f - 2*h_x*bc_neumann_baixo
            f_imj = anticorpo_anterior[i-1][j] if i != 0 else f - 2*h_x*bc_neumann_cima
            f_ijp = anticorpo_anterior[i][j+1] if j != tam-1 else f - 2*h_x*bc_neumann_direita
            f_ijm = anticorpo_anterior[i][j-1] if j != 0 else f - 2*h_x*bc_neumann_esquerda            

            #Dependendo do gradiente dos ODCs vou fazer upwind ou downwind no eixo i ou eixo j

            #Decidindo qual combinacao usar no gradiente das células com quimiotaxia
            
            gradiente_odc_i = (olide_ipj - olide_imj)/(2*h_x)
            gradiente_odc_j = (olide_ijp - olide_ijm)/(2*h_x)

            if gradiente_odc_i > 0:
                gradiente_m_i = gradiente(m, mac_imj)/h_x
                gradiente_dc_i = gradiente(dc, dc_imj)/h_x
                gradiente_da_i = gradiente(da, da_imj)/h_x
                gradiente_t_i = gradiente(t_cito, t_cito_imj)/h_x
            else:
                gradiente_m_i = gradiente(mac_ipj, m)/h_x
                gradiente_dc_i = gradiente(dc_ipj, dc)/h_x
                gradiente_da_i = gradiente(da_ipj, da)/h_x
                gradiente_t_i = gradiente(t_cito_ipj, t_cito)/h_x
            if gradiente_odc_j > 0:
                gradiente_m_j = gradiente(m, mac_ijm)/h_x
                gradiente_dc_j = gradiente(dc, dc_ijm)/h_x
                gradiente_da_j = gradiente(da, da_ijm)/h_x
                gradiente_t_j = gradiente(t_cito, t_cito_ijm)/h_x
            else:
                gradiente_m_j = gradiente(mac_ijp, m)/h_x
                gradiente_dc_j = gradiente(dc_ijp, dc)/h_x
                gradiente_da_j = gradiente(da_ijp, da)/h_x
                gradiente_t_j = gradiente(t_cito_ijp, t_cito)/h_x

            #Dados da equacao macrofagos
            quimiotaxia_mac = chi*(gradiente_odc_i*gradiente_m_i + gradiente_odc_j*gradiente_m_j)
            difusao_mac = (mac_ipj + mac_imj - 4*m + mac_ijp + mac_ijm )/h_x**2
            reacao_mac = m*(1 - m)
            
            mac_atual[i][j] = m + h_t*(difusao_mac + reacao_mac - quimiotaxia_mac)
            
            #T citotóxica
            quimiotaxia_t_cito = chi*(gradiente_odc_i*gradiente_t_i + gradiente_odc_j*gradiente_t_j)
            difusao_t_cito = d_t_cit*(t_cito_ijm + t_cito_ijp - 4*t_cito + t_cito_imj + t_cito_ipj)/h_x**2

            t_cito_atual[i][j] = t_cito + h_t*(difusao_t_cito - quimiotaxia_t_cito + gamma_tcito*(t_cito_linf - t_cito)) 

            #Dados da equacao oligodendrocitos destruidos 
            olide_atual[i][j] = oligo_destr + h_t*((r_m + lamb_f_m*f)*f_func(m)*(1 - oligo_destr) + r_dc*f_func(dc)*oligo_destr + r_t*f_func(t_cito)*(1 - oligo_destr))

            # anticorpo
            difusao_anticorpo = d_anti*(f_ipj + f_imj - 4*f + f_ijp + f_ijm)

            anticorpo_atual[i][j] = f + h_t*(difusao_anticorpo - lamb_f_dc*f*(1-oligo_destr)*f_func(dc) - lamb_f_m*f*(1-oligo_destr)*f_func(m) + gamma_at*(at_linf - f))

            #DC convencional
            quimiotaxia_dc = chi*(gradiente_odc_i*gradiente_dc_i + gradiente_odc_j*gradiente_dc_j)
            difusao_dc = d_dc*(dc_ipj + dc_imj - 4*dc + dc_ijp + dc_ijm )/h_x**2
            reacao_dc = mu*oligo_destr*(d_c_homeostase- dc)

            dendritica_conv_atual[i][j] = dc + h_t*(reacao_dc + difusao_dc - quimiotaxia_dc - b_d*oligo_destr*dc)

            #DC ativada
            difusao_da = d_da*(da_ipj + da_imj - 4*da + da_ijp + da_ijm )/h_x**2

            dendritica_ativ_atual[i][j] = da + h_t*(difusao_da + b_d*oligo_destr*dc - alpha_d*(da_linf-dc)) 

            if(mac_atual[i][j]>1 or mac_atual[i][j]<0):
                print('tempo: '+str(k*h_t))
                print("Mac:"+str(mac_atual[i][j]))

            if(t_cito_atual[i][j]>1 or t_cito_atual[i][j]<0):
                print('tempo: '+str(k*h_t))
                print("t_cito:"+str(t_cito_atual[i][j]))

            if(olide_atual[i][j]>1 or olide_atual[i][j]<0):
                print('tempo: '+str(k*h_t))
                print("olide:"+str(olide_atual[i][j]))

            if(anticorpo_atual[i][j]>1 or anticorpo_atual[i][j]<0):
                print('tempo: '+str(k*h_t))
                print("anticorpo:"+str(anticorpo_atual[i][j]))

            if(dendritica_conv_atual[i][j]>1 or dendritica_conv_atual[i][j]<0):
                print('tempo: '+str(k*h_t))
                print("dendritica_conv_atual:"+str(dendritica_conv_atual[i][j]))

            if(dendritica_ativ_atual[i][j]>1 or dendritica_ativ_atual[i][j]<0):
                print('tempo: '+str(k*h_t))
                print("dendritica_ativ_atual:"+str(dendritica_ativ_atual[i][j]))
            
    olide_anterior = np.copy(olide_atual)
    dendritica_conv_anterior = np.copy(dendritica_conv_atual)
    dendritica_ativ_anterior = np.copy(dendritica_ativ_atual)
    t_cito_anterior = np.copy(t_cito_atual)
    anticorpo_anterior = np.copy(anticorpo_atual)
    mac_anterior = np.copy(mac_atual)

    if k%intervalo_figs ==0 or k == steps-1:
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