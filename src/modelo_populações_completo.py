import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import time
import os
from linfonodo import diferential

sns.set()
os.system("python3 ../criaDiretorios.py")

gradiente = lambda ponto_anterior, ponto_posterior, valor_maximo: quimiotaxia(ponto_posterior, valor_maximo) - quimiotaxia(ponto_anterior, valor_maximo)
quimiotaxia = lambda ponto_atual, valor_maximo: ponto_atual/(valor_maximo + ponto_atual)
f_func = lambda populacao, valor_maximo: populacao*populacao/(valor_maximo + populacao)


T_final = 1# Dia
h_t = 0.001

L = 10  # Comprimento da malha
h_x = 0.1

t = np.linspace(0, T_final, int(T_final/h_t))
x = np.linspace(0, L, int(L/h_x))
tam = len(x)
steps = len(t)

num_figuras = 5
intervalo_figs = int(steps/num_figuras)

def verifica_cfl(difusao_mic, difusao_dc, difusao_da, quimiotaxia_dc, quimiotaxia_mic):
    if(difusao_mic*h_t/h_x**2 < 1/4 and difusao_dc*h_t/h_x**2 < 1/4 and difusao_da*h_t/h_x**2 < 1/4 and quimiotaxia_dc*h_t/h_x < 1 and quimiotaxia_mic*h_t/h_x < 1):
        return True
    else:
        return False

V_BV = 0
V_LV = 0

theta_BV = np.zeros((int(L/h_x), int(L/h_x)))
for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if i > L/h_x*.9 and j > L/h_x*.9:
            theta_BV[i][j] = 1
            V_LV += 1

theta_LV = np.zeros((int(L/h_x), int(L/h_x)))
for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if i < L/h_x*.5 and j < L/h_x*.5:
            theta_LV[i][j] = 1
            V_BV += 1

V_LN = 160

# IC
# microglia
mic_media = 350
mic_anterior = np.zeros((int(L/h_x), int(L/h_x)))

for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if (i-int(L/h_x)/2)**2 + (j-int(L/h_x)/2)**2 < 20:
            mic_anterior[i][j] = mic_media/3.0

# T citotóxica
t_cito_anterior = np.zeros((int(L/h_x), int(L/h_x)))
# Ol destruidos
olide_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# anticorpo
anticorpo_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# Dendríticas convencionais
dc_media = 5
dendritica_conv_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# Dendríticas ativadas
dendritica_ativ_anterior = np.zeros((int(L/h_x), int(L/h_x)))

#***********************Declaracao das matrizes que vao guardar os valores do passo de tempo atual*********************

mic_atual = np.zeros((int(L/h_x), int(L/h_x)))
t_cito_atual = np.zeros((int(L/h_x), int(L/h_x)))
olide_atual = np.zeros((int(L/h_x), int(L/h_x)))
anticorpo_atual = np.zeros((int(L/h_x), int(L/h_x)))
dendritica_conv_atual = np.zeros((int(L/h_x), int(L/h_x)))
dendritica_ativ_atual = np.zeros((int(L/h_x), int(L/h_x)))

# Modelo linfonodo
estable_B = 8.4*10**-4
linfonodo_eqs = np.zeros(5)
linfonodo_eqs[0]= 0    # Dendritic cells
linfonodo_eqs[1]= 0.2  # Cytotoxic T cells
linfonodo_eqs[2]= 0.4  # Helper T cells
linfonodo_eqs[3]= estable_B    # B cells
linfonodo_eqs[4]= 0    # Antibodies

#Valores das populaçoes que migram que estão em contato com os vasos sanguineos ou linfaticos
DendriticasTecido = 0
AnticorposTecido = 0
TcitotoxicaTecido = 0

for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if theta_LV[i][j] == 1:
            DendriticasTecido += dendritica_ativ_anterior[i][j]
        if theta_BV[i][j] == 1:
            AnticorposTecido += anticorpo_anterior[i][j]
            TcitotoxicaTecido += t_cito_anterior[i][j]

#**********************Funcao print dos resultados*************************

populationTitle = {
    "odc": "Oligodendrócitos destruídos",
    "microglia": "Microglia",
    "dc": "Dendrítica convencional",
    "da": "Dendrítica ativada",
    "tke": "T citotóxica",
    "anticorpo": "Anticorpos"
}

def printMesh(time, population, type):

    x_pts, y_pts = np.meshgrid(x, x)
    max_population = np.max(population)
    if max_population == 0:
        max_population += 1
    levels = np.linspace(0, max_population, 10)

    cp = plt.contourf(x_pts, y_pts,population, levels=levels)
    plt.title(populationTitle[type])
    plt.xlabel("Milímetros")
    plt.ylabel("Milímetros")
    plt.colorbar(cp, label="Concentração (células/$mm^2$)")
    plt.savefig('../results/'+type+'/fig'+'{:.4f}'.format(time*h_t)+'.png', dpi = 300)
    plt.clf()


parameters = {
    "chi": 0.298*60*2 , # Quimioatracao(a mesma para todas as celulas por enquanto). valor por Dia
    "D_mic": 60*24*6.6*10**-5, # Difusao da microglia. valor por Dia
    "mu_m": 60*24*3*10**-6, # Taxa de ativação da microglia. valor por Dia
    "r_m": 60*24*3.96*10**-6, # intensidade dos danos causados pela microglia valor por Dia

    "d_dc": 60*24*6.6*10**-5, # difusao DC convencional(procurar na literatura)
    "d_da": 60*24*6.6*10**-5, # difusao DC ativada(procurar na literatura)
    "d_t_cit": 60*24*6.6*10**-5, # difusao t citotóxica(procurar na literatura)
    "d_anti": 60*24*6.6*10**-5, # difusao anticorpo(procurar na literatura)
    "lamb_f_m": 60*24*3.96*10**-6, # taxa de anticorpos consumidos durante o processo de opsonização pela micróglia
    "b_d": 0.001, # taxa de ativacao de dc por odc destruidos(procurar na literatura)
    "r_dc": 0.001, # taxa de coleta de odc destruidos pelas DCs (procurar na literatura)
    "r_t": 0.1 , # agressividade de t citotoxica(procurar na literatura)

    "mu_dc": 60*24*3*10**-6, #Taxa de producao de células dendríticas (procurar na literatura)
    "gamma_D": 0.001, #Taxa de migração de DC ativadas para o linfonodo (procurar na literatura)
    "gamma_F": 0.0003, #Taxa de migração de anticorpos para o tecido (procurar na literatura)
    "gamma_T": 0.2, #Taxa de migração de T citotoxica para o tecido (procurar na literatura)

    "t_cito_media": 37,
    "dc_media": dc_media,
    "mic_media": mic_media,
    "odc_media": 400,


    "alpha_T_h": 0.01 ,
    "alpha_T_c": 0.5,
    "alpha_B": 1,
    "b_T": 0.017,
    "b_rho": 10**5,
    "b_rho_b": 6.02*10**3,
    "rho_T": 2,
    "rho_B": 16,
    "rho_F": 5.1*10**4,
    "estable_T_h": 8.4*10**-3,
    "estable_B": estable_B,
    "estable_T_c": 8.4*10**-3,
    "DendriticasTecido": DendriticasTecido,
    "AnticorposTecido": AnticorposTecido,
    "TcitotoxicaTecido": TcitotoxicaTecido,
    "V_LV": V_LV,
    "V_BV": V_BV,
    "V_LN": V_LN
}
if not verifica_cfl(parameters["D_mic"], parameters["d_dc"], parameters["d_da"], parameters["chi"], parameters["chi"]):
    print("Falhou cfl!!!")
    exit(1)

#BC
bc_neumann_cima = 0
bc_neumann_direita = 0
bc_neumann_baixo = 0
bc_neumann_esquerda = 0

DL_vetor = np.zeros(steps)
TL_c_vetor = np.zeros(steps)
TL_h_vetor = np.zeros(steps)
B_vetor = np.zeros(steps)
FL_vetor = np.zeros(steps)

DL_vetor[0] = linfonodo_eqs[0]
TL_c_vetor[0] = linfonodo_eqs[1]
TL_h_vetor[0] = linfonodo_eqs[2]
B_vetor[0] = linfonodo_eqs[3]
FL_vetor[0] = linfonodo_eqs[4]

printMesh(0,olide_anterior, "odc")
printMesh(0,mic_anterior, "microglia")
printMesh(0,dendritica_conv_anterior, "dc")
printMesh(0,dendritica_ativ_anterior, "da")
printMesh(0,t_cito_anterior, "tke")
# printMesh(0,anticorpo_anterior, "anticorpo")

#Inicio da contagem do tempo
tic = time.perf_counter()

for k in range(1,steps):
    dy = diferential(linfonodo_eqs, parameters)
    DL_atual = linfonodo_eqs[0] + h_t*dy[0]
    TL_c_atual = linfonodo_eqs[1] + h_t*dy[1]
    TL_h_atual = linfonodo_eqs[2] + h_t*dy[2]
    B_atual = linfonodo_eqs[3] + h_t*dy[3]
    FL_atual = linfonodo_eqs[4]# + h_t*dy[4]
    
    for i in range(tam):
        for j in range(tam):
            oligo_destr = olide_anterior[i][j]
            microglia = mic_anterior[i][j]
            dc = dendritica_conv_anterior[i][j]
            da = dendritica_ativ_anterior[i][j]
            # anticorpo = anticorpo_anterior[i][j]
            t_cito = t_cito_anterior[i][j]
            
            # condição de contorno de Neumman microglia
            mic_ipj = mic_anterior[i+1][j] if i != tam-1 else microglia - 2*h_x*bc_neumann_baixo
            mic_imj = mic_anterior[i-1][j] if i != 0 else microglia - 2*h_x*bc_neumann_cima
            mic_ijp = mic_anterior[i][j+1] if j != tam-1 else microglia - 2*h_x*bc_neumann_direita
            mic_ijm = mic_anterior[i][j-1] if j != 0 else microglia - 2*h_x*bc_neumann_esquerda

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
            
            # ponto fantasma oligodendrocitos destruidos ODC Nao é contorno
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
            # f_ipj = anticorpo_anterior[i+1][j] if i != tam-1 else anticorpo - 2*h_x*bc_neumann_baixo
            # f_imj = anticorpo_anterior[i-1][j] if i != 0 else anticorpo - 2*h_x*bc_neumann_cima
            # f_ijp = anticorpo_anterior[i][j+1] if j != tam-1 else anticorpo - 2*h_x*bc_neumann_direita
            # f_ijm = anticorpo_anterior[i][j-1] if j != 0 else anticorpo - 2*h_x*bc_neumann_esquerda            

            #Dependendo do gradiente dos ODCs vou fazer upwind ou downwind no eixo i ou eixo j

            #Decidindo qual combinacao usar no gradiente das células com quimiotaxia
            
            gradiente_odc_i = (olide_ipj - olide_imj)/(2*h_x)
            gradiente_odc_j = (olide_ijp - olide_ijm)/(2*h_x)

            if gradiente_odc_i > 0:#@@@@ TODO Fazer teste sem anticorpo. Depois voltar para quando so tinha DC-mic-ods-DA
                gradiente_m_i = gradiente(microglia, mic_imj, parameters["mic_media"])/h_x
                gradiente_dc_i = gradiente(dc, dc_imj,parameters["dc_media"])/h_x
                gradiente_t_i = gradiente(t_cito, t_cito_imj, parameters["t_cito_media"])/h_x
            else:
                gradiente_m_i = gradiente(mic_ipj, microglia, parameters["mic_media"])/h_x
                gradiente_dc_i = gradiente(dc_ipj, dc, parameters["dc_media"])/h_x
                gradiente_t_i = gradiente(t_cito_ipj, t_cito, parameters["t_cito_media"])/h_x
            if gradiente_odc_j > 0:
                gradiente_m_j = gradiente(microglia, mic_ijm, parameters["mic_media"])/h_x
                gradiente_dc_j = gradiente(dc, dc_ijm, parameters["dc_media"])/h_x
                gradiente_t_j = gradiente(t_cito, t_cito_ijm, parameters["t_cito_media"])/h_x
            else:
                gradiente_m_j = gradiente(mic_ijp, microglia, parameters["mic_media"])/h_x
                gradiente_dc_j = gradiente(dc_ijp, dc, parameters["dc_media"])/h_x
                gradiente_t_j = gradiente(t_cito_ijp, t_cito, parameters["t_cito_media"])/h_x

            #Dados da equacao microglia
            quimiotaxia_mic = parameters["chi"]*(gradiente_odc_i*gradiente_m_i + gradiente_odc_j*gradiente_m_j)
            difusao_mic = parameters["D_mic"]*(mic_ipj + mic_imj - 4*microglia + mic_ijp + mic_ijm )/h_x**2
            reacao_mic = parameters["mu_m"]*microglia*(parameters["mic_media"] - microglia)
            
            mic_atual[i][j] = microglia + h_t*(difusao_mic + reacao_mic - quimiotaxia_mic)

            #T citotóxica
            quimiotaxia_t_cito = 0#parameters["chi"]*(gradiente_odc_i*gradiente_t_i + gradiente_odc_j*gradiente_t_j)
            difusao_t_cito = parameters["d_t_cit"]*(t_cito_ijm + t_cito_ijp - 4*t_cito + t_cito_imj + t_cito_ipj)/h_x**2
            migracao_t_cito = theta_BV[i][j]*parameters["gamma_T"]*(TL_c_atual - t_cito)
            
            t_cito_atual[i][j] = t_cito + h_t*(difusao_t_cito - quimiotaxia_t_cito + migracao_t_cito)

            #Oligodendrocitos destruidos 
            fag_mic_ant = 0#parameters["lamb_f_m"]*anticorpo
            apoptose_tke = parameters["r_t"]*f_func(t_cito, parameters["t_cito_media"])*(parameters["odc_media"] - oligo_destr)
            olide_atual[i][j] = oligo_destr + h_t*((parameters["r_m"] + fag_mic_ant)*f_func(microglia, mic_media)*(parameters["odc_media"] - oligo_destr) + apoptose_tke)

            #Anticorpo
            #difusao_anticorpo = parameters["d_anti"]*(f_ipj + f_imj - 4*anticorpo + f_ijp + f_ijm)
            #reacao_anticorpo = parameters["lamb_f_m"]*anticorpo*(parameters["odc_media"] - oligo_destr)*f_func(microglia, mic_media)
            #migracao_anticorpo = theta_BV[i][j]*parameters["gamma_F"]*(FL_atual - anticorpo)

            #anticorpo_atual[i][j] = anticorpo + h_t*(difusao_anticorpo - reacao_anticorpo + migracao_anticorpo)

            #DC convencional
            quimiotaxia_dc = 0# parameters["chi"]*(gradiente_odc_i*gradiente_dc_i + gradiente_odc_j*gradiente_dc_j)
            difusao_dc = parameters["d_dc"]*(dc_ipj + dc_imj - 4*dc + dc_ijp + dc_ijm )/h_x**2
            reacao_dc = parameters["mu_dc"]*oligo_destr*(parameters["dc_media"] - dc)
            ativacao_dc_da = parameters["b_d"]*oligo_destr*dc

            dendritica_conv_atual[i][j] = dc + h_t*(reacao_dc + difusao_dc - quimiotaxia_dc - ativacao_dc_da)
            
            #DA ativada
            difusao_da = parameters["d_da"]*(da_ipj + da_imj - 4*da + da_ijp + da_ijm)/h_x**2
            migracao_da = theta_LV[i][j]*parameters["gamma_D"]*(DL_atual - da)

            dendritica_ativ_atual[i][j] = da + h_t*(difusao_da + ativacao_dc_da + migracao_da)
            if microglia < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel microglia: " + str(microglia))
            if da < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel DA: " + str(da))
            if dc < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel dc: " + str(dc))
            if t_cito < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel t_cito: " + str(t_cito))
            # if anticorpo < 0:
            #     print("Tempo do Erro: " + str(k*h_t) + " - Variavel anticorpo: " + str(anticorpo))
            if oligo_destr < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel oligo_destr: " + str(oligo_destr))
            
    olide_anterior = np.copy(olide_atual)
    dendritica_conv_anterior = np.copy(dendritica_conv_atual)
    dendritica_ativ_anterior = np.copy(dendritica_ativ_atual)
    t_cito_anterior = np.copy(t_cito_atual)
    anticorpo_anterior = np.copy(anticorpo_atual)
    mic_anterior = np.copy(mic_atual)
    
    #Atualização da concentração das populações que migram.
    #Valores das populaçoes que migram que estão em contato com os vasos sanguineos ou linfaticos
    DendriticasTecido = 0
    AnticorposTecido = 0
    TcitotoxicaTecido = 0

    for i in range(int(L/h_x)):
        for j in range(int(L/h_x)):
            if theta_LV[i][j] == 1:
                DendriticasTecido += dendritica_ativ_anterior[i][j]
            if theta_BV[i][j] == 1:
                AnticorposTecido += anticorpo_anterior[i][j]
                TcitotoxicaTecido += t_cito_anterior[i][j]

    parameters["TcitotoxicaTecido"] = TcitotoxicaTecido
    parameters["DendriticasTecido"] = DendriticasTecido
    parameters["AnticorposTecido"] = AnticorposTecido

    linfonodo_eqs = [DL_atual, TL_c_atual, TL_h_atual, B_atual, FL_atual]
    DL_vetor[k] = DL_atual
    TL_c_vetor[k] = TL_c_atual
    TL_h_vetor[k] = TL_h_atual
    B_vetor[k] = B_atual
    FL_vetor[k] = FL_atual

    if k%intervalo_figs ==0 or k == steps-1:
        printMesh(k,olide_anterior, "odc")
        printMesh(k,mic_anterior, "microglia")
        printMesh(k,dendritica_conv_anterior, "dc")
        printMesh(k,dendritica_ativ_anterior, "da")
        printMesh(k,t_cito_anterior, "tke")
        # printMesh(k,anticorpo_anterior, "anticorpo")
        print("Tempo: "+ str(k*h_t))

#Fim da contagem do tempo
toc = time.perf_counter()

final_time = (toc - tic)/60

print("Tempo de execução: " + str(final_time) + " min")

#Transforma de dias para horas no plot
t = np.multiply(t,24)

plt.plot(t,DL_vetor)
plt.title("Dendríticas ativadas no linfonodo")
plt.xlabel("Tempo (horas)")
plt.ylabel("Concentração (células/$mm^2$)")
plt.savefig('../results/dc_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,TL_c_vetor)
plt.title("T citotóxicas no linfonodo")
plt.xlabel("Tempo (horas)")
plt.ylabel("Concentração (células/$mm^2$)")
plt.savefig('../results/t_cito_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,TL_h_vetor)
plt.title("T helper no linfonodo")
plt.xlabel("Tempo (horas)")
plt.ylabel("Concentração (células/$mm^2$)")
plt.savefig('../results/t_helper_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,B_vetor)
plt.title("Células B no linfonodo")
plt.xlabel("Tempo (horas)")
plt.ylabel("Concentração (células/$mm^2$)")
plt.savefig('../results/b_cell_linfonodo.png', dpi = 300)
plt.clf()

# plt.plot(t,FL_vetor)
# plt.title("Anticorpos no linfonodo")
# plt.xlabel("Tempo (horas)")
# plt.ylabel("Concentração (células/$mm^2$)")
# plt.savefig('../results/anticorpo_linfonodo.png', dpi = 300)
# plt.clf()

