# --- IMPORT DEPENDENCIES ------------------------------------------------------+

import numpy as np
from depend.cost import cost
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, NonlinearConstraint
import os
import time
# --- MAIN ---------------------------------------------------------------------+


# --- CONSTANTS ----------------------------------------------------------------+

cost_func = cost                                  #Cost function
bounds = [(0.000001,0.5),(0.000001,0.5),(0.000001,0.5),(0.000001,0.5),(0.000001,0.5)] 
popsize = 1000                                               #Population size
# Mutation factor [0,2]
mutate = 0.7
# Recombination rate [0,1]
recombination = 0.5
# Max number of generations (maxiter)
maxiter = 8

# --- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    
    cwd = os.getcwd()
    saida = open(cwd+"/differential_evolution/relatorio_DE.txt", "a")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: ' + str(popsize) +
                     '\nNumber of generations: ' + str(maxiter) + '\n')
    tic = time.perf_counter()
    sol_pat = differential_evolution(cost_func, bounds,maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination)
    tac = time.perf_counter()
    # sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
    saida.writelines(f"Tempo de execução: {tac - tic:0.4f} segundos")
    saida.writelines(
        '\nCusto do melhor conjunto de parametros: ' + str(sol_pat.fun) + '\n\n')
    saida.writelines(
        '\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
    print('\nCusto do melhor conjunto de parametros: ' +
            str(sol_pat.fun) + '\n\n')
    print('\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
    saida.close()

    
# --- END ----------------------------------------------------------------------+
