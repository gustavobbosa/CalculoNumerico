from . import auxiliar as aux

# ################################################################

def eliminacao_gauss(A,b,metodo=aux.pivoteamento_parcial,imprimir = False):
    
    '''Retorna uma lista com os xn.'''
    
    if imprimir: print("\n## Começando eliminação de Gauss ##\n")
    
    A, b = aux.matriz(A), aux.matriz(b)
    aux.escalonador_de_sistemas(A,b,metodo,imprimir)
    n = len(A)
    respostas = [0 for _ in range(n)]
    i,j = n,n
    
    while i > 0:
        
        soma = 0
        while j > i:
            soma += A[i-1][j-1] * respostas[j-1]
            j -= 1
            
        respostas[i-1] = (b[i-1] - soma) / A[i-1][i-1]
        
        j = n
        i -= 1
    
    if imprimir: print("x=",respostas,"\n## Eliminação finalizada ##\n")
    return respostas

# #############################################################

def gauss_jacobi (
    A,b,x0,
    limite=1000,
    epsilon=0.05,
    imprimir=False,
    debug=False
    ):
    
    '''Executa o método de Gauss-Jacobi até a precisão epsilon e retorna um dicionário com os xn.'''
    
    if imprimir: print(f"\n## Executando Gauss-Jacobi ##\nA={A}\nb={b}")
    
    # Fazendo C e g
    C,g = aux.produz_c_e_g(A,b,imprimir,debug)
    
        # Verifica critério das linhas
    if aux.criterio_das_linhas(A,imprimir=imprimir) == False:
        if input("Critério das linhas não cumprido. Continuar? (s/n) ") == "n": 
            raise Exception("Operação cancelada pelo usuário")
    
    n = len(A)
    
    k = 0
    # processo iterativo
    while True:
        
        if imprimir: print(f"Etapa {k}:\nx = {x0}")
        
        if debug: print(f"DEBUG\nC = {C}\nx0 = {x0}\ng = {g}\n")
        x1 = C.dot(x0) + g
        if debug: print("x1 = ",x1)
        
        d = []
        
        for i in range(n): d.append(abs(x1[i]-x0[i]))
        if debug: print(f"Vetor d é {d}")
        
        if sorted(d)[-1] < epsilon:
            if imprimir: print(f"Etapa {k} (final):\nx = {x1}")
            return aux.matriz(x1)
        
        if k > 1000: raise ValueError("Não pára!")
        
        k += 1
        x0 = x1
        
# ##########################################################

def gauss_seidel (
    A,b,x0,
    limite=1000,
    epsilon=0.05,
    imprimir=False,
    debug=False
    ):
    
    '''Executa o método de Gauss-Seidel até a precisão epsilon e retorna um dicionário com os xn.'''
    
    if imprimir: print(f"\n## Executando Gauss-Seidel ##\nA={A}\nb={b}")
    
    # Fazendo C e g
    C,g = aux.produz_c_e_g(A,b,imprimir,debug)
    
    # Verifica critério das linhas
    if aux.criterio_sassenfeld(A,imprimir=imprimir) == False:
        if input("Critério de Sassenfeld não cumprido. Continuar? (s/n) ") == "n": 
            raise Exception("Operação cancelada pelo usuário")
    
    n = len(A)
    
    k = 0
    # processo iterativo
    while True:
        
        if imprimir: print(f"Etapa {k}:\nx = {x0}")
        
        if debug: print(f"DEBUG\nC = {C}\nx0 = {x0}\ng = {g}\n")
        x_temp = aux.matriz([el for el in x0])
        
        for i in range(n):
            x_i = C[i].dot(x_temp) + g[i]
            x_temp[i] = x_i
            
        x1 = x_temp
        if debug: print("x1 = ",x1)
        
        d = []
        
        for i in range(n): d.append(abs(x1[i]-x0[i]))
        if debug: print(f"Vetor d é {d}")
        
        if sorted(d)[-1] < epsilon:
            if imprimir: print(f"Etapa {k} (final):\nx = {x1}")
            return aux.matriz(x1)
        
        if k > limite: raise ValueError("Não pára!")
        
        k += 1
        x0 = x1
