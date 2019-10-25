'''Funções auxiliares para os métodos'''
import numpy as np

# ###### Utilitários para matrizes #######

def matriz(*args, **kwargs):
    
    '''Faz matrizes com float64 no lugar de int32. 
Para mais detalhes, ver numpy.array'''
    
    kwargs.setdefault("dtype", np.float64)
    return np.array(*args, **kwargs)

def faz_matriz_quadrada(lista):
    
    '''Tranforma uma lista com um tamanho quadrado em uma matriz quadrada.'''
    
    n = len(lista)**(1/2)
    
    if not n.is_integer():
        raise ValueError("Tamanho da lista não é um quadrado")
    
    return matriz(lista).reshape(n,n)    
    
def verifica_solução_sistema (A,b,x,tol=0.01,imprimir=False):
    
    '''

Verifica se x é solução de Ax = b, levando em conta certa toleerância (tol). 
Para cada equação do sistema, avalia o membro esquerdo e compara com o direito.
Se, para todas equações, o erro relativo for menor que a tol, retorna verdadeiro.
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

# ###### Funções auxiliares dos métodos #######

def pivoteamento_parcial(A,b,k,debug = False):
    
    '''Escolhe o melhor pivô dentre os candidatos da coluna k e faz a troca de linhas na matriz e no vetor in place.'''
    
    if debug: print(f"Pivoteamento etapa {k}:\n",A,b)
    n = len(A)
    
    # Selecionando os candidatos na forma de tuples (índice, elemento)
    candidatos = [(i,A[i][k]) for i in range(k,n)]
    
    # Determina o vencedor ao organizá-los em ordem crescente do módulo do elemnto
    candidatos.sort( key = lambda x: abs(x[1]) )
    
    # Vencedor é o índice da linha do maior elemento
    vencedor = candidatos[-1][0]
    
    # Se o vencedor já está na posição k, não precisa trocar nada
    if vencedor == k: return
    
    # Trocando as linhas da matriz. Tá escrito de um jeito zoado pra impedir que
    # a variável linha_original mude com A[k].
    linha_original = [el for el in A[k]]
    A[k] = A[vencedor]
    A[vencedor] = linha_original
    
    # Trocando os elementos no vetor
    el_original = b[k]
    b[k] = b[vencedor]
    b[vencedor] = el_original
    
    if debug: print(f"Saída:\n",A,b)

def pivoteamento_total(A,b,k,x0):
    pass

def escalonador_de_sistemas(A,b,metodo=pivoteamento_parcial,imprimir=False,**kwargs):
        
    '''Transforma a matriz A em uma matriz triangular superior. Possibilidade de escolher o método de pivoteamento'''
    
    if imprimir: print("\n## Escalonando sistema... ##\n")
     
    n = len(A)
    
    # Limpando cada coluna k em cada etapa k
    for k in range(n):
        
        # Método de pivoteamento. Tem que ser "in place" pra mexer nas linhas e colunas da matriz.
        metodo(A,b,k,**kwargs)
        pivo = A[k][k]
        
        # Processando cada linha i
        for i in range(k+1,n):
            # determinando o multiplicador para a linha i
            multiplicador = A[i][k] / pivo
            # processando cada elemento aᵢⱼ e bᵢ
            b[i] -= multiplicador * b[k]
            for j in range(n): 
                A[i][j] -= multiplicador * A[k][j]
    
    if imprimir: print(f"A = {A}\nb={b}\n## Fim do escalonamento ##")

def criterio_das_linhas(A,imprimir=False):
    
    '''Lê os alfas da matriz A e vê se eles são menores que 1. Pode imprimir cada alfa'''
    
    n = len(A)
    
    saida = True
    
    for i in range(n):
        
        soma = 0
        
        # Soma só os elementos que não estejam na diagonal
        for j in range(n): 
            if not j==i: soma += abs(A[i][j])
        
        alfa = soma / abs(A[i][i])
        
        if imprimir: print(f"α{i+1} = {alfa}")
        
        # Se algum alfa for maior que 1, retorna falso. Caso contrário, termina o for e retorna True
        if alfa >= 1: saida = False
        
    return saida

def criterio_sassenfeld(A,imprimir=False):
    
    '''Lê os betas da matriz A e vê se eles são menores que 1. Pode imprimir os betas'''
    
    n = len(A)
    
    betas = []
    
    saida = True
    
    for i in range(n):
        
        soma = 0
        
        for j in range(n):
            if j<i: soma += betas[j]*abs(A[i][j])
            if j>i: soma += abs(A[i][j])
        
        beta = soma / abs(A[i][i])
        
        if imprimir: print(f"β{i+1} = {beta}")
        
        betas.append(beta)
        
        # Se algum beta for maior que 1, retorna falso. Caso contrário, termina o for e retorna True
        if beta >= 1: saida = False
        
    return saida

def produz_c_e_g (A,b,imprimir=False,debug=False):
    
    '''Retorna um tuple com a matriz C e o vetor g de um sistema linear.'''
    
    if imprimir: print("\n## Fazendo C e G ##\n")
    n, C, g = len(A), [], []
    
    # ### Verifica as merdas ####
    zero_na_diagonal = True
    passos = 0
    
    # verificando zeros na diagonal
    while zero_na_diagonal:    
        for i in range(n):
            if A[i][i] == 0:
                
                if imprimir: print(f"Elemento {i+1} da diagonal é zero", A)
                if debug: print(f"Sistema antes de mexer no passo {passos}\nA={A}\nb={b}")
                
                temp = [ el for el in A[i]]
                A[i] = A[i-1]
                A[i-1] = temp
                temp = b[i]
                b[i] = b[i-1]
                b[i-1] = temp
                
                if debug: print(f"Sistema depois de mexer no passo {passos}\nA={A}\nb={b}")
                
                passos += 1
                break
            zero_na_diagonal = False
            
        if passos > 10*n: raise Exception("Não é possível evitar zeros na diagonal")
            
    if imprimir:
        print(f"Matrizes prontas para fazer C e g:\nC={C}\ng={g}")
    # ####### Fim da verificação de merdas ###########
    
    # para cada linha da matriz, ou seja, para cada equação
    for k in range(n):
        
        linha = []
        
        # para cada elemento da linha 
        for j in range(n):
            
            # Colocar os zeros ou as frações
            if j == k: linha.append(0)
            else: linha.append(- A[k][j] / A[k][k])
        
        # Colocar a linha pronta na matriz
        C.append(linha)
        
        # Mexendo no vetor
        g.append(b[k] / A[k][k])
    
    if imprimir: print("C e g produzidos:\nC={C}\ng={g}\n## Função produz_c_e_g terminada ##")
        
    return matriz(C), matriz(g)
