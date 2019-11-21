# coding: utf-8

import numpy as np

def gera_dados(X,Y):
    '''Entrada: duas listas, uma com os X e uma com os Y. Essa função pareia cada X com cada Y e retorna uma lista com os tuples de pontos em ordem crescente de X.'''
    return sorted( list( zip(X,Y) ) )   
    
def verifica_interpolacao (dados,polinomio,tol=0.01,imprimir=False):
    '''Verifica se o polinômio (tem que ser executável) interpola os pontos em dados, com tolerância tol (padrão 1%)'''
    
    for x,y in dados:

        if imprimir: print(f"{polinomio(x)} ?= {y}")
        if y ==0:
            if polinomio(x) > tol: return False
        else:
            if abs (1 - polinomio(x)/y) > tol: return False            
    return True

def diferencas_divididas(dados,imprimir=False):
    '''
Recebe como entrada uma lista de dados como ensinado em CalculoNumerico.Interpolacao.metodos

Retorna uma matriz do numpy que representa a tabela de diferenças divididas.
O elemento [-1,0] é a diferença dividida da lista toda.
    '''
    
    ord_max = len(dados)
    
    # Criando uma tabela vazia
    tabela_dd = np.zeros((ord_max,ord_max))
    
    # Encontrando as dd de ordem zero, ou seja, só o valor da função
    for i in range(ord_max):
        tabela_dd[0,i] = dados[i][1]
        if imprimir: print(tabela_dd,'\n')
        
    # Determinando para todas outras ordens a partir das anteriores
    for ordem in range(1,ord_max):
        # não é pra preencher toda a linha, mas um a menos que a linha anterior.
        for i in range(ord_max - ordem): 
            numerador = tabela_dd[ordem - 1, i] - tabela_dd[ordem - 1, i + 1] 
            denominador = dados[i][0] - dados[i+ordem][0]
            tabela_dd[ordem,i] = numerador/denominador
            if imprimir: print(tabela_dd,'\n')
            
    return tabela_dd
    
def dados_para_funcao (dados):
    '''
Converte uma lista de tuples em uma função f(i) = j para cada tuple (i,j)
    '''

    dados_dict = {}
    
    for x,y in dados: dados_dict.update({f"{np.float64(x)}":y})
    
    def func(x):
        return dados_dict[f"{np.float64(x)}"]
    
    return func

def base_polinomial(n):
    '''
Retorna uma lista com polinômios puros de grau 0 a n
    '''
    base = []
    for i in range(n+1):
        base.append(np.poly1d([1] + i * [0]))
    
    return base



def diferencas_divididas_recursivo(dados, func=None, imprimir=False):
    '''
FUNÇÃO DEPRECADA, não consegui fazer funcionar pra dar a tabela inteira

Calcula diferença dividida da lista de dados. Dois modos de uso:
    
    Modo 1: dados é uma lista de tuples e não há argumento para func.
    Modo 2: dados é uma lista de números e func deve ser um objeto executável (função, polinômio numpy, etc.)
    '''
    
    if imprimir: impressao = f'f{dados} = ' 
    # Verifica se há apenas um elemento em dados, senão a função se autorreferencia até o infinito e dá RecursionError.
    if len(dados) == 1:
        if func == None: return dados[0][1]
        else: return func(dados[0])
    
    # Senão, função se executa com listas menores
    else:
        numerador = diferencas_divididas(dados[0:-1],func=func,imprimir=imprimir) - diferencas_divididas(dados[1:],func=func,imprimir=imprimir)
        
        if func == None: denominador = dados[0][0] - dados[-1][0]
        else: denominador = dados[0] - dados[-1]
        
        resultado = numerador/denominador
        
        if imprimir: 
            impressao += f'{resultado}'
            print(impressao)
        
        return resultado