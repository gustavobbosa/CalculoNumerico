# coding: utf-8

import sympy as sym

def bisseccao(func,a,b,epsilon=0.05,imprimir=True):
    '''
Computa a raiz no intervalo [a,b] de uma função de um argumento e saída numérica. Retorna a raiz.

epsilon é o limite do tamanho do intervalo [a,b].
Pode imprimir as etapas na tela.
    '''
    if b <= a: raise ValueError("B deve ser maior que A!")
    if epsilon <= 0: raise ValueError("Precisão deve ser positiva!")
    k = 1
    
    if b - a < epsilon: # Verifica condição de parada
            return a
    
    if imprimir: print("Iteração    x\t\tf(x)\t\ta\t\tb\t\tb-a\n") # Cabeçalho da impressão
    
    while True:

        m = func(a)

        x = (a+b)/2 # Computa o ponto médio
        
        if imprimir: print(f"{k:2}\t{x:.8f}\t{func(x):.4e}\t{a:.8f}\t{b:.8f}\t{b-a:.4e}") # imprime iteração
        
        if b - a < epsilon: # Verifica condição de parada
            return x
        
        if m*func(x)>0: # Determina qual ponto excluir
            a = x
        else:
            b = x
        
        k += 1
		
###################################################################

# Tem duas precisões, uma pro intervalo no domínio e um pra imagem
def posicao_falsa(func,a,b,epsilon1=0.05,epsilon2=0.05,imprimir=True):
    
    '''
Coomputa a raiz no intervalo [a,b] de uma função de um argumento e saída numérica. Retorna a raiz.

epsilon1 é o limite do tamanho do intervalo [a,b].
epsilon2 é a precisão do valor absoluto da função.
    '''
    if b <= a: raise ValueError("B deve ser maior que A!")
    if (epsilon1 <= 0) or (epsilon2 <= 0) : raise ValueError("Precisão deve ser positiva!")
    k = 1
    
    # Verifica precisões antes de começar
    if b - a < epsilon1:
        return a
    if abs(func(a)) < epsilon2:
        return a
    if abs(func(b)) < epsilon2:
        return b
    
    if imprimir: print("Iteração    x\t\tf(x)\t\ta\t\tb\t\tb-a\n") # Imprime cabeçalho
    
    while True:
        
        fa = func(a)
        fb = func(b)
        
        # A diferença do peso é aqui
        x = ( a * abs(fb) + b * abs(fa) ) / ( abs(fb) + abs(fa) )
        
        if imprimir: print(f"{k:2}\t{x:.8f}\t{func(x):.4e}\t{a:.8f}\t{b:.8f}\t{b-a:.4e}")

        # Verifica parada
        if b - a < epsilon1:
            return x
        
        # Também se verifica precisão na imagem 
        if abs(fa) < epsilon2:
            return a
        if abs(fb) < epsilon2:
            return b
        
        # Determina ponto a ser excluído
        if fa*func(x)>0:
            a = x
        else:
            b = x
        
        k += 1
		
#######################################################################

def ponto_fixo(func,phi,x0,epsilon1=1e-5,epsilon2=1e-5,imprimir=True,k_max=200,f_max=1000000):
    
    '''
Calcula raiz de func usando o método do ponto fixo e a função phi. Retorna a raiz.

requer ponto inicial x0
epsilon1: precisão da imagem
epsilon2: precisão do domínio
k_max e f_max: Limites para evitar que o código trave. k_max é numero máximo de passos e f_max é valor máximo da imagem.
    '''
    
    x = x0
    
    ponto_fixo.k = 0 # k é atributo para poder ser lido após execução da função (Questão 3)
    
    if imprimir: print("Iteração\t\t\t\tx\t\t\t\tf(x)")
    
    while True:
        
        f = func(x)
        p = phi(x)
        
        # Verifica parada
        if abs(f) < epsilon1:
            print(f"\nPronto! f(x) = {f} < {epsilon1} = epsilon1")
            return x
        if abs(x - p) < epsilon2:
            print(f"\nPronto! |x-φ(x)| = {x-p} < {epsilon2} = epsilon2")
            return p
        
        # Evitar que o código trave
        if abs(f) > f_max or ponto_fixo.k > k_max: 
            print("Saiu do controle!")
            return
        
        ponto_fixo.k += 1
        
        if imprimir: print (f"{ponto_fixo.k}\t\t\t\t{x:.8f}\t\t\t\t{f:.8f}")

        x = p # phi(x) se torna o novo x
        
#########################################################################

def newton_raphson(func,x0,deriv=None,**kwargs):
    
    '''
Calcula a raiz de func pelo método de newton. 
Pode-se usar como entrada a função e sua derivada como funções do python, ou só a função como expressão sympy.

Usa os mesmos parâmetros que a função ponto_fixo
    '''
    
    # Se não houver derivada na entrada, calculá-la a partir de func. Só possível se func for expressão sympy.
    if deriv == None: 
        if 'sympy' in str(type(func)):
            deriv = func.diff()
        else:
            raise TypeError("Falta a derivada e não consigo calculá-la!")
    
    # Tranformar as expressões sympy em funções, caso seja necessário
    if 'sympy' in str(type(func)):
        def func_exec(t):
            return float(func.subs(func.free_symbols.pop(),t))
    else:
        func_exec = func
    if 'sympy' in str(type(deriv)):
        def deriv_exec(t):
            return float(deriv.subs(deriv.free_symbols.pop(),t))
    else:
        deriv_exec = deriv
    
    # Definir phi(x) de acordo com o método de Newton Raphson
    def phi(x):
        return float(x-func_exec(x)/deriv_exec(x))
    resposta = ponto_fixo(func_exec,phi,x0,**kwargs)
    
    newton_raphson.k = ponto_fixo.k # k como atributo pra poder pegar o número final de passos caso necessário (Questão 3)
    return resposta
        
#################################################################################

def secante(func,x0,x1,epsilon1=1e-5,epsilon2=1e-5,imprimir=True):
    
    if abs(func(x0)) < epsilon1:
        return x0
    
    if (abs(func(x1)) < epsilon1) or (abs( x1 - x0 ) < epsilon2):
        return x1

    k = 1
    
    if imprimir: print(f"Iteração\t\t\t\tx\t\t\t\tf(x)")
    
    while True:
        
        f1 = func(x1)
        f0 = func(x0)
        
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        f2 = func(x2)
        
        if imprimir: print(f"{k}\t\t\t\t{x2:.10f}\t\t\t{f2:.5e}")
        
        if (abs(f2) < epsilon1) or (abs ( x2 - x1 ) < epsilon2):
            return x2
        
        x0 = x1
        x1 = x2
        
        k += 1