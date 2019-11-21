'''
Métodos para efetuar interpolação polinomial de um conjunto de dados


UTILIZANDO O PACOTE Interpolacao.metodos:

As funções fazem interpolação polinomial dos pontos contidos na lista de entrada dados. Os pontos estão na forma de tuple.

Então a entrada deve ser uma lista de tuples como o exemplo:

EXEMPLO:
    dados = [ (0,1), (1,2.5), (2.2,3.7)]

    # pontos x=0, y=1
    #             x=1, y=2.5
    #             x=2.2, y=3.7
    #
    
Crie o tuple e execute a função.

'''

from ..Sistemas import metodos as sl
from . import auxiliar as aux
import numpy as np


def sistema(dados, metodo=sl.eliminacao_gauss,**kwargs):
    
    grau = len(dados) - 1

    A = sl.aux.matriz([[dados[i][0]**j for j in range(grau+1)] for i in range(grau+1)])
    
    b = sl.aux.matriz([dados[i][1] for i in range(len(dados))])
    
    alfas = list(reversed(metodo(A,b,**kwargs)))
    
    return np.poly1d(alfas)
    
def lagrange( dados ):
    
    pn = np.poly1d([0])
    n = len(dados) - 1
    
    for k in range(n + 1):
        
        Lk = np.poly1d([1])
        
        for j in range(n + 1):
            
            if j != k: Lk = Lk * np.poly1d([1, -dados[j][0]]) /  (dados[k][0] - dados[j][0] )
            
        pn += dados[k][1] * Lk
        
    return pn
    
def newton ( dados, imprimir=False ):
    
    global tabela_dd
    
    '''
Entrada tal qual a doctring do módulo CalculoNumerico.Interpolacao.metodos (este arquivo)
    
Saída: polinômio de uma dimensão do numpy.
    '''
    
    n = len(dados) - 1
    pn = dados[0][1] # primeiro termo do polinômio é f(x₀)
    
    for k in range (n): # para cada termo restante do polinômio
        
        # Calcula a diferença dividida dos pontos entre 0 e k+1. Executando aux.diferencas_divididas no modo 1
        dif_div = aux.diferencas_divididas(dados[0:k+2],imprimir=imprimir)
        
        # Executa o produtório (x-x₁)·(x-x₂)· ... ·(x-xₖ)
        prod = np.poly1d([1])
        for j in range(k+1):
            prod = prod * np.poly1d([1,-dados[j][0]])
        
        pn += dif_div * prod # Soma o elemento no polinômio final
            
    return pn
   
def quadrados_discreto (dados,grau,metodo=sl.eliminacao_gauss,**kwargs):
    
    m = len(dados)
    if m <= grau: raise ValueError("m deve ser maior que n")
    
    A,b = [],[]
    func = aux.dados_para_funcao(dados)
    g = aux.base_polinomial(grau + 1)
    
    for i in range(grau + 1):
        
        bi = sum( [ y * g[i](x) for x,y in dados ] )
        b.append(bi)
        
        Ai = []
        for j in range(grau + 1):
            Aij = sum( [ g[i](x) * g[j](x) for x,y in dados ] )
            Ai.append(Aij)
        A.append(Ai)
    
    alfas = metodo(A,b,**kwargs)
    
    return sum( [ float(alfas[i]) * g[i] for i in range(grau + 1) ] ) 