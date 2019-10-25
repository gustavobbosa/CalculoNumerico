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

def diferencas_divididas(dados, func=None):
    '''Calcula diferença dividida da lista de dados. Dois modos de uso:
    
    Modo 1: dados é uma lista de tuples e não há argumento para func.
    Modo 2: dados é uma lista de números e func deve ser um objeto executável (função, polinômio numpy, etc.)
    '''
    
    if len(dados) == 1:
        if func == None: return dados[0][1]
        else: return func(dados[0])
    
    else:
        numerador = diferencas_divididas(dados[0:-1],func=func) - diferencas_divididas(dados[1:],func=func)
        
        if func == None: denominador = dados[0][0] - dados[-1][0]
        else: denominador = dados[0] - dados[-1]
        
        return numerador / denominador
    
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