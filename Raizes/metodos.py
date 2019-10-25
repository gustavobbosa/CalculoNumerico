def bisseccao(func,a,b,epsilon=0.05,imprimir=False):
    
    k = 1
    
    if b - a < epsilon:
            return a
    
    if imprimir: print("Iteração    x\t\tf(x)\t\ta\t\tb\t\tb-a\n")
    
    while True:

        m = func(a)

        x = (a+b)/2
        
        if imprimir: print(f"{k:2}\t{x:.8f}\t{func(x):.4e}\t{a:.8f}\t{b:.8f}\t{b-a:.4e}")
        
        if b - a < epsilon:
            return x
        
        if m*func(x)>0:
            a = x

        else:
            b = x
        
        k += 1
		
###################################################################

# Tem duas precisões, uma pro intervalo no domínio e um pra imagem
def posicao_falsa(func,a,b,epsilon1=0.05,epsilon2=0.05,imprimir=False):
    
    k = 1
    
    if b - a < epsilon1:
        return a
    
    if abs(func(a)) < epsilon2:
        return a
        
    if abs(func(b)) < epsilon2:
        return b
    
    if imprimir: print("Iteração    x\t\tf(x)\t\ta\t\tb\t\tb-a\n")
    
    while True:
        
        fa = func(a)
        fb = func(b)
        
        # A diferença do peso é aqui
        x = ( a * abs(fb) + b * abs(fa) ) / ( abs(fb) + abs(fa) )
        
        if imprimir: print(f"{k:2}\t{x:.8f}\t{func(x):.4e}\t{a:.8f}\t{b:.8f}\t{b-a:.4e}")

        if b - a < epsilon1:
            return x
        
        # Também se verifica precisão na imagem 
        if abs(fa) < epsilon2:
            return a
        
        if abs(fb) < epsilon2:
            return b
        
        if fa*func(x)>0:
            a = x
        
        else:
            b = x
        
        k += 1
		
#######################################################################

def ponto_fixo(func,phi,x0,epsilon1=1e-5,epsilon2=1e-5,imprimir=True,k_max=200,f_max=1000000):
    
    x = x0
    
    k = 0
    
    if imprimir: print("Iteração\t\t\t\tx\t\t\t\tf(x)")
    
    while True:
        
        f = func(x)
        p = phi(x)
        
        if abs(f) < epsilon1:
            
            print(f"\nPronto! f(x) = {f} < {epsilon1} = epsilon1")
            return x
        
        if abs(x - p) < epsilon2:
            
            print(f"\nPronto! |x-φ(x)| = {x-p} < {epsilon2} = epsilon2")
            return p
        
        if abs(f) > f_max or k > k_max: 
            
            print("Saiu do controle!")
            return
        
        k += 1
        
        if imprimir: print (f"{k}\t\t\t\t{x:.8f}\t\t\t\t{f:.8f}")

        x = p
        
#########################################################################

def newton_raphson(func,deriv,x0,**kwargs):
    
    '''Função e derivada como funções normais do python, não precisa de bibliotecas externas.
Opção imprimir lista os valores de x e f(x) em cada etapa'''
    
    def phi(x):
        return x-func(x)/deriv(x)
    
    return metodo_ponto_fixo(func,phi,x0,**kwargs)
        
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