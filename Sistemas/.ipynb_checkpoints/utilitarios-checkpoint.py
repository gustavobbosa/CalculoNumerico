'''Utilitários para matrizes'''
import numpy as np

def array(*args, **kwargs):
    
    '''Faz matrizes com float64 no lugar de int32. 
Para mais detalhes, ver numpy.array'''
    
    kwargs.setdefault("dtype", np.float64)
    return np.array(*args, **kwargs)

def verifica_quadrada(lista):
    
    '''Retorna a raiz do tamannho da lista. Levanta um erro caso o tamanho da lista não seja um número quadrado.'''
    
    n = len(lista)**(1/2)
    
    if not n.is_integer():
        raise ValueError("Tamanho da lista não é um quadrado")
        
    return int(n)

def vetor_para_matriz (vetor,coluna=True):
    
    '''Verifica se o vetor é uma matriz. Se não for, retorna a matriz correspondente.
Por padrão será uma matriz coluna, mas isso pode ser sobrescrito'''
    
    try:
        # Se ele for uma lista, esse try vai dar erro.
        vetor[0][0] = vetor[0][0]
        return vetor
    
    except:
        if coluna: return array([vetor]).transpose()
        else: return array([vetor])

def faz_matriz_quadrada(lista):
    
    '''Tranforma uma lista com um tamanho quadrado em uma matriz quadrada.'''
    
    n = verifica_quadrada(lista)
    
    quadrada = array(lista).reshape(n,n)    
    return quadrada

def imprime_matriz(matriz):
    
    '''Imprime uma matriz linha por linha'''
    
    matriz = vetor_para_matriz(matriz,coluna=False)
    saida = '\n'
    
    for linha in matriz:
        
        string_linha = '|'
        
        for elemento in linha:
            string_linha += f'\t{elemento:.3e}'
        
        string_linha += '\t|'
        
        saida += string_linha + '\n'
        
    print (saida)

def verifica_solução_sistema (A,x,b,tol=0.01,imprimir=False):
    
    '''

Verifica se x é solução de Ax = b, levando em conta um erro de ε. 
Para cada equação do sistema, avalia o membro esquerdo e compara com o direito.
Se, para todas equações, o erro relativo for menor que a tolerância (tol), retorna verdadeiro.
Opção imprimir mostra cada avaliação enquanto nenhuma for falsa.
    
    '''
    
    n = len(A)
    
    for i in range(n):
        lado_esquerdo = 0
        
        for j in range(n):
            lado_esquerdo += A[i][j] * x[j]
            
        if imprimir: print(lado_esquerdo,'?=',b[i])
        
        if abs( 1 - lado_esquerdo/b[i]) > tol:
            return False
        
    return True
