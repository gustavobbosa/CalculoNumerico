from ..Sistemas import metodossl as sl
import numpy as np

def interpolação_polinomial_sistema(dados, metodo=sl.eliminacao_gauss,**kwargs):
    grau = len(dados) + 1

    A = sl.ut.array([[dados[i][0]**j for j in range(len(dados))] for i in range(len(dados))])
    
    b = sl.ut.vetor_para_matriz([dados[i][1] for i in range(len(dados))])
    
    alfas = list(reversed(metodo(A,b,**kwargs)))
    
    g = np.poly1d(alfas)
    
    return g