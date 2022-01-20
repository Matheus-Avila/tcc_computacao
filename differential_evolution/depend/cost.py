import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()

def diferential(y, parameters):
    # DT Dendritic cells in the tissue
    # FL Average number of antibodies on tissue 
    V_LV = parameters[0]     # Volume of the domain areas in contact with lymph vessels
    V_LN = parameters[1]     # Volume of the lymph node
    V_BV = parameters[2]
    
    gamma_D = parameters[3] # Dendritic cells migration rate (Tissue -> Lymph node)
    alpha_T_c = parameters[4]  
    estable_T_c = parameters[5]
    gamma_T = parameters[6]    
    theta_BV = parameters[7]
    Tt_c = parameters[8] # T citotoxic from tissue

    b_T = parameters[9]    
    rho_T = parameters[10]    
    b_rho = parameters[11]    
    alpha_T_h = parameters[12]    
    estable_T_h = parameters[13]
    rho_B = parameters[14]    #* Replication rate of B cells
    alpha_B = parameters[15]  #* Replication rate of B cells 
    estable_B = parameters[16]# Estable B value
    
    rho_F = parameters[17]  # Antibodies production rate (by B cells) 
    gamma_F = parameters[18]# Antibodies migration rate (to tissue)
    DT = parameters[19]
    FT = parameters[20]
    b_rho_b = parameters[21]

    # When working with only one diferential function there must be created a null extra position
    dy = np.zeros(5)

    # Dendritic cells
    dy[0] = gamma_D * (DT - y[0]) * (V_LV / V_LN)
    # Cytotoxic T cells
    dy[1] = alpha_T_c * (estable_T_c - y[1]) - (gamma_T * theta_BV * (y[1] - Tt_c)) * (V_BV / V_LN)
    # Helper T cells
    dy[2] = b_T*(rho_T * y[2] * y[0] - y[2]*y[0]) - (b_rho * y[2] * y[0] * y[3]) + alpha_T_h * (estable_T_h - y[2])
    # B cells
    dy[3] = (b_rho_b * ((rho_B * y[2] * y[0]) - (y[2] * y[0] * y[3]))) + alpha_B * (estable_B - y[3])
    # Antibodies
    dy[4] = rho_F * y[3] - ((gamma_F * theta_BV * (y[4] - FT)) * (V_BV / V_LN))

    return dy

def cost(poi):
    print("Nova tentativa!!!")
    gradiente = lambda ponto_anterior, ponto_posterior, valor_maximo: quimiotaxia(ponto_posterior, valor_maximo) - quimiotaxia(ponto_anterior, valor_maximo)
    quimiotaxia = lambda ponto_atual, valor_maximo: ponto_atual/(valor_maximo + ponto_atual)
    f_func = lambda populacao, valor_maximo: populacao*populacao/(valor_maximo + populacao)

    T_final = 1 # Dia
    h_t = 0.05

    L = 10#25.8  # Comprimento da malha
    # L = 100
    h_x = 0.1# 0.05

    chi = poi[5]# 0.298*60*24  # Quimioatracao(a mesma para todas as celulas por enquanto). valor por Dia
    D_mac = poi[6] # 60*24*6.6*10**-5 # Difusao da microglia. valor por Dia
    mu_m = poi[7] # 60*24*3*10**-6 # Taxa de ativação da microglia. valor por Dia
    r_m = poi[8] # 60*24*3.96*10**-6 # intensidade dos danos causados pela microglia valor por Dia

    d_dc = poi[0] # difusao DC convencional(procurar na literatura)
    d_da = poi[1] # difusao DC ativada(procurar na literatura)
    d_t_cit = D_mac # difusao t citotóxica(procurar na literatura)
    d_anti = D_mac # difusao anticorpo(procurar na literatura)
    lamb_f_m = 7.14*10**-2 # taxa de anticorpos consumidos durante o processo de opsonização pela micróglia ( precisa converter!!!) 
    b_d = poi[2] # taxa de ativacao de dc por odc destruidos(procurar na literatura)
    r_dc = 0.1 # taxa de coleta de odc destruidos pelas DCs (procurar na literatura)
    r_t = 0.1  # agressividade de t citotoxica(procurar na literatura)

    mu_dc = poi[3] #Taxa de producao de células dendríticas (procurar na literatura)
    gamma_d = poi[4] #Taxa de migração de DC ativadas para o linfonodo (procurar na literatura)
    gamma_anticorpo = 0.43 #Taxa de migração de anticorpos para o tecido 
    gamma_tcito = 0.3 #Taxa de migração de T citotoxica para o tecido 

    t_cito_media = 37
    dc_media = 5 
    mac_media = 350
    odc_media = 400

    V_LV = 1
    V_LN = 1
    V_BV = 1
    alpha_T_h = 0.01 # tese barbara
    alpha_T_c = 0.5
    alpha_B = 1
    b_T = 0.017
    b_rho = 10**5
    b_rho_b = 6.02*10**3
    rho_T = 2
    rho_B = 16
    rho_F = 5.1*10**4
    estable_T_h = 8.4*10**-3
    estable_B = 8.4*10**-4
    estable_T_c = 8.4*10**-3
    theta_BV = 1


    # IC
    # Macrofagos
    mac_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    for i in range(int(L/h_x)):
        for j in range(int(L/h_x)):
            if (i-int(L/h_x)/2)**2 + (j-int(L/h_x)/2)**2 < 20:
                mac_anterior[i][j] = mac_media/3.0

    # T citotóxica
    t_cito_anterior = np.zeros((int(L/h_x), int(L/h_x)))
    # Ol destruidos
    olide_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    # anticorpo
    anticorpo_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    # Dendríticas convencionais
    dendritica_conv_anterior = dc_media*np.ones((int(L/h_x), int(L/h_x)))

    # Dendríticas ativadas
    dendritica_ativ_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    mac_atual = np.zeros((int(L/h_x), int(L/h_x)))
    t_cito_atual = np.zeros((int(L/h_x), int(L/h_x)))
    olide_atual = np.zeros((int(L/h_x), int(L/h_x)))
    anticorpo_atual = np.zeros((int(L/h_x), int(L/h_x)))
    dendritica_conv_atual = np.zeros((int(L/h_x), int(L/h_x)))
    dendritica_ativ_atual = np.zeros((int(L/h_x), int(L/h_x)))

    # Modelo linfonodo
    linfonodo_eqs = np.zeros(5)
    linfonodo_eqs[0]= 0    # Dendritic cells
    linfonodo_eqs[1]= 0.2  # Cytotoxic T cells
    linfonodo_eqs[2]= 0.4  # Helper T cells
    linfonodo_eqs[3]= 0    # B cells
    linfonodo_eqs[4]= 0    # Antibodies

    DT = np.sum(dendritica_ativ_anterior)
    FT = np.sum(anticorpo_anterior)
    Tt_c = np.sum(t_cito_anterior)

    parameters = np.zeros(22)

    parameters[0] = V_LV     # Volume of the domain areas in contact with lymph vessels
    parameters[1] = V_LN     # Volume of the lymph node
    parameters[2] = V_BV

    parameters[3] = gamma_d # Dendritic cells migration rate (Tissue -> Lymph node)
    parameters[4] = alpha_T_c  
    parameters[5] = estable_T_c
    parameters[6] = gamma_tcito    
    parameters[7] = theta_BV
    parameters[8] = Tt_c

    parameters[9] = b_T    
    parameters[10] = rho_T     
    parameters[11] = b_rho     
    parameters[12] = alpha_T_h     
    parameters[13] = estable_T_h 
    parameters[14] = rho_B     #* Replication rate of B cells
    parameters[15] = alpha_B   #* Replication rate of B cells 
    parameters[16] = estable_B # Estable B value

    parameters[17] = rho_F   # Antibodies production rate (by B cells) 
    parameters[18] = gamma_anticorpo # Antibodies migration rate (to tissue)
    parameters[19] = DT 
    parameters[20] = FT 
    parameters[21] = b_rho_b


    DL_atual = 0     
    TL_c_atual = 0 
    TL_h_atual = 0 
    B_atual = 0      
    FL_atual = 0     

    #BC
    bc_neumann_cima = 0
    bc_neumann_direita = 0
    bc_neumann_baixo = 0
    bc_neumann_esquerda = 0

    t = np.linspace(0, T_final, int(T_final/h_t))
    x = np.linspace(0, L, int(L/h_x))
    tam = len(x)
    steps = len(t)

    DL_vetor = np.zeros(steps)
    TL_c_vetor = np.zeros(steps)
    TL_h_vetor = np.zeros(steps)
    B_vetor = np.zeros(steps)
    FL_vetor = np.zeros(steps)
    custo = 0
    for k in range(1,steps):
        dy = diferential(linfonodo_eqs, parameters)
        DL_atual = linfonodo_eqs[0] + h_t*dy[0]
        TL_c_atual = linfonodo_eqs[1] + h_t*dy[1]
        TL_h_atual = linfonodo_eqs[2] + h_t*dy[2]
        B_atual = linfonodo_eqs[3] + h_t*dy[3]
        FL_atual = linfonodo_eqs[4] + h_t*dy[4]

        for i in range(tam):
            for j in range(tam):
                oligo_destr = olide_anterior[i][j]
                microglia = mac_anterior[i][j]
                dc = dendritica_conv_anterior[i][j]
                da = dendritica_ativ_anterior[i][j]
                anticorpo = anticorpo_anterior[i][j]
                t_cito = t_cito_anterior[i][j]
                
                # condição de contorno de Neumman macrofagos
                mac_ipj = mac_anterior[i+1][j] if i != tam-1 else microglia - 2*h_x*bc_neumann_baixo
                mac_imj = mac_anterior[i-1][j] if i != 0 else microglia - 2*h_x*bc_neumann_cima
                mac_ijp = mac_anterior[i][j+1] if j != tam-1 else microglia - 2*h_x*bc_neumann_direita
                mac_ijm = mac_anterior[i][j-1] if j != 0 else microglia - 2*h_x*bc_neumann_esquerda

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
                f_ipj = anticorpo_anterior[i+1][j] if i != tam-1 else anticorpo - 2*h_x*bc_neumann_baixo
                f_imj = anticorpo_anterior[i-1][j] if i != 0 else anticorpo - 2*h_x*bc_neumann_cima
                f_ijp = anticorpo_anterior[i][j+1] if j != tam-1 else anticorpo - 2*h_x*bc_neumann_direita
                f_ijm = anticorpo_anterior[i][j-1] if j != 0 else anticorpo - 2*h_x*bc_neumann_esquerda            

                #Dependendo do gradiente dos ODCs vou fazer upwind ou downwind no eixo i ou eixo j

                #Decidindo qual combinacao usar no gradiente das células com quimiotaxia
                
                gradiente_odc_i = (olide_ipj - olide_imj)/(2*h_x)
                gradiente_odc_j = (olide_ijp - olide_ijm)/(2*h_x)

                if gradiente_odc_i > 0:
                    gradiente_m_i = gradiente(microglia, mac_imj, mac_media)/h_x
                    gradiente_dc_i = gradiente(dc, dc_imj,dc_media)/h_x
                    gradiente_t_i = gradiente(t_cito, t_cito_imj, t_cito_media)/h_x
                else:
                    gradiente_m_i = gradiente(mac_ipj, microglia, mac_media)/h_x
                    gradiente_dc_i = gradiente(dc_ipj, dc, dc_media)/h_x
                    gradiente_t_i = gradiente(t_cito_ipj, t_cito, t_cito_media)/h_x
                if gradiente_odc_j > 0:
                    gradiente_m_j = gradiente(microglia, mac_ijm, mac_media)/h_x
                    gradiente_dc_j = gradiente(dc, dc_ijm, dc_media)/h_x
                    gradiente_t_j = gradiente(t_cito, t_cito_ijm, t_cito_media)/h_x
                else:
                    gradiente_m_j = gradiente(mac_ijp, microglia, mac_media)/h_x
                    gradiente_dc_j = gradiente(dc_ijp, dc, dc_media)/h_x
                    gradiente_t_j = gradiente(t_cito_ijp, t_cito, t_cito_media)/h_x


                #Dados da equacao macrofagos
                quimiotaxia_mac = chi*(gradiente_odc_i*gradiente_m_i + gradiente_odc_j*gradiente_m_j)
                difusao_mac = D_mac*(mac_ipj + mac_imj - 4*microglia + mac_ijp + mac_ijm )/h_x**2
                reacao_mac = mu_m*microglia*(mac_media - microglia)
                
                mac_atual[i][j] = microglia + h_t*(difusao_mac + reacao_mac - quimiotaxia_mac)

                #T citotóxica
                quimiotaxia_t_cito = chi*(gradiente_odc_i*gradiente_t_i + gradiente_odc_j*gradiente_t_j)
                difusao_t_cito = d_t_cit*(t_cito_ijm + t_cito_ijp - 4*t_cito + t_cito_imj + t_cito_ipj)/h_x**2
                migracao_t_cito = gamma_tcito*(TL_c_atual - t_cito)
                
                t_cito_atual[i][j] = t_cito + h_t*(difusao_t_cito - quimiotaxia_t_cito + migracao_t_cito)

                #Oligodendrocitos destruidos 
                olide_atual[i][j] = oligo_destr + h_t*((r_m + lamb_f_m*anticorpo)*f_func(microglia, mac_media)*(odc_media - oligo_destr) + r_dc*f_func(dc, dc_media)*oligo_destr + r_t*f_func(t_cito, t_cito_media)*(odc_media - oligo_destr))

                #Anticorpo
                difusao_anticorpo = d_anti*(f_ipj + f_imj - 4*anticorpo + f_ijp + f_ijm)
                reacao_anticorpo = -lamb_f_m*anticorpo*(odc_media - oligo_destr)*f_func(microglia, mac_media)
                migracao_anticorpo = gamma_anticorpo*(FL_atual - anticorpo)

                anticorpo_atual[i][j] = anticorpo + h_t*(difusao_anticorpo + reacao_anticorpo + migracao_anticorpo)

                #DC convencional
                quimiotaxia_dc = chi*(gradiente_odc_i*gradiente_dc_i + gradiente_odc_j*gradiente_dc_j)
                difusao_dc = d_dc*(dc_ipj + dc_imj - 4*dc + dc_ijp + dc_ijm )/h_x**2
                reacao_dc = mu_dc*oligo_destr*(dc_media- dc)

                dendritica_conv_atual[i][j] = dc + h_t*(reacao_dc + difusao_dc - quimiotaxia_dc - b_d*oligo_destr*dc)
                
                #DC ativada
                difusao_da = d_da*(da_ipj + da_imj - 4*da + da_ijp + da_ijm )/h_x**2

                dendritica_ativ_atual[i][j] = da + h_t*(difusao_da + b_d*oligo_destr*dc + gamma_d*(DL_atual - da))

                if( mac_atual[i][j]<0):
                    print('tempo: '+str(k*h_t))
                    print("Mac:"+str(mac_atual[i][j]))
                    custo = 10
                    return custo

                if(t_cito_atual[i][j]<0):
                    print('tempo: '+str(k*h_t))
                    print("t_cito:"+str(t_cito_atual[i][j]))
                    custo = 10
                    return custo

                if(olide_atual[i][j]<0):
                    print('tempo: '+str(k*h_t))
                    print("olide:"+str(olide_atual[i][j]))
                    custo = 10
                    return custo

                if(anticorpo_atual[i][j]<0):
                    print('tempo: '+str(k*h_t))
                    print("anticorpo:"+str(anticorpo_atual[i][j]))
                    custo = 10
                    return custo

                if(dendritica_conv_atual[i][j]<0):
                    print('tempo: '+str(k*h_t))
                    print("dendritica_conv_atual:"+str(dendritica_conv_atual[i][j]))
                    custo = 10
                    return custo

                if(dendritica_ativ_atual[i][j]<0):
                    print('tempo: '+str(k*h_t))
                    print("dendritica_ativa_atual:"+str(dendritica_ativ_atual[i][j]))
                    custo = 10
                    return custo
                

        olide_anterior = np.copy(olide_atual)
        dendritica_conv_anterior = np.copy(dendritica_conv_atual)
        dendritica_ativ_anterior = np.copy(dendritica_ativ_atual)
        t_cito_anterior = np.copy(t_cito_atual)
        anticorpo_anterior = np.copy(anticorpo_atual)
        mac_anterior = np.copy(mac_atual)
        DT = np.sum(dendritica_ativ_anterior)
        FT = np.sum(anticorpo_anterior)
        Tt_c = np.sum(t_cito_anterior)
        parameters[8] = Tt_c
        parameters[19] = DT
        parameters[20] = FT

        linfonodo_eqs = [DL_atual, TL_c_atual, TL_h_atual, B_atual, FL_atual]
        DL_vetor[k] = DL_atual
        TL_c_vetor[k] = TL_c_atual
        TL_h_vetor[k] = TL_h_atual
        B_vetor[k] = B_atual
        FL_vetor[k] = FL_atual
    
    return custo