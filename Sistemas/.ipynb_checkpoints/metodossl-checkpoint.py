from . import utilitarios as ut
from . import auxiliar as aux

###############################################################################

def eliminacao_gauss(A,b,metodo=aux.pivoteamento_parcial,imprimir = False):
    
    '''Retorna uma lista com os xn.'''
    
    if imprimir: print("\n## Começando eliminação de Gauss ##\n")
        
    if not type(b[0]) == list:
        b = ut.vetor_para_matriz(b)
    A = ut.array(A)
    aux.escalonador_de_sistemas(A,b,metodo,imprimir)
    
    if imprimir:
        print("Matrizes escolonadas:")
        ut.imprime_matriz(A)
        ut.imprime_matriz(b)
    
    n = len(A)
    respostas = [0 for _ in range(n)]
    i,j = n,n
    
    while i > 0:
        
        soma = 0
        while j > i:
            soma += A[i-1][j-1] * respostas[j-1]
            j -= 1
            
        respostas[i-1] = (b[i-1][0] - soma) / A[i-1][i-1]
        
        j = n
        i -= 1
    
    if imprimir: print("\n## Eliminação finalizada ##\n")
    return respostas

##############################################################

def gauss_jacobi (A,b,x0,epsilon=0.05,imprimir=False,debug=False,**kwargs):
    
    '''Executa o método de Gauss-Jacobi até a precisão epsilon e retorna um dicionário com os xn.'''
    
    if imprimir: 
        print("\n## Executando Gauss-Jacobi ##\n")
        ut.imprime_matriz(A)
        ut.imprime_matriz(b)
        
    # Fazendo C e g
    C,g = aux.produz_c_e_g(A,b,imprimir=imprimir,**kwargs)
    
    n = len(A)
    
    k = 0
    # processo iterativo
    while True:
        
        if imprimir:
            print(k)
            ut.imprime_matriz(x0)
        
        if debug: print(f"DEBUG\nC = {C}\nx0 = {x0}\ng = {g}\n")
        x1 = C.dot(x0) + g
        if debug: print("x1 = ",x1)
            
        d = []
        
        for i in range(n): d.append(abs(x1[i]-x0[i]))
        if debug: print(f"vetor d é {d}")
        
        if sorted(d)[-1] < epsilon:
            if imprimir:
                print(k+1)
                ut.imprime_matriz(x1)
            return ut.vetor_para_matriz(x1)
        
        if k > 1000: raise ValueError("Não pára!")
        
        k += 1
        x0 = x1
        
############################################################################################3

def gauss_seidel (A,b,x0,epsilon=0.05,imprimir=False):
    
    '''Executa o método de Gauss-Seidel até a precisão epsilon e retorna um dicionário com os xn.'''
    
    # Fazendo C e g
    C,g = aux.produz_c_e_g(A,b)
    
    n = len(A)
    
    k = 0
    # processo iterativo
    while True:
        
        if imprimir:
            print(k)
            ut.imprime_matriz(x0)
        
        ####### Isso muda tudo!!! ###########
        # x1 = C.dot(x0) + g
        x_temp = ut.array([el for el in x0])
        
        for i in range(n):
            x_i = C[i].dot(x_temp) + g[i]
            x_temp[i] = x_i
            
        x1 = x_temp
        #####################################
        
        d = []
        
        for i in range(n): d.append(abs(x1[i]-x0[i]))
            
        if sorted(d)[-1] < epsilon:
            if imprimir:
                print(k+1)
                ut.imprime_matriz(x1)
            return ut.ver_é_matriz(x1)
        
        # Caso dê merda e a parada vá ao infinito
        if k > 10000: raise ValueError("Não pára!")
        
        k += 1
        
        x0 = x1
