import numpy as np
import matplotlib.pyplot as plt

'''
criar função que sera usada
    ela é continua?
        qual o dominio
Plotar funçõa
Tabela de Valores
Raiz Numerica
    qual a derivada
'''


class RaizNumerica:
    '''
    Variaveis: dominio:                 dominio, 
               função:                  funcao,
               derivada da função:      derivFuncao,
               função phi:              phi, 
               derivada de phi:         derivPhi,
               valor inicial:           v,
               tolerancia:              tol,
               epsolon:                 eps,
               limite inferoir:         a,
               limite superior:         b,
               quantidade de valores:   n,
               limite de iterações:     max_iter

    Metódos: 
        Dados da função:
            Mostrar Função: exibe o grafica da função a partir do dominio definido
            Tabelar valores: são tabelados valoers que estão definidos 
                             entre um limite inferior e um limite superior e 
                             pode ser definido uma quantidade de valores entre eles.
        Métodos numéricos: 
            são os metodos utilizados para obtenção de raizes numericas dentro de um intervalo dado, 
            cada um com suas particularidades e convergencias.
            
            bisseccao,
            ponto_fixo,
            falsa_posicao,
            newton-raphson,
            secante   
    '''
    
    
    def __init__(self,dominio = None, funcao = None, derivFuncao = None, phi = None, derivPhi = None):
        
        self.dominio = dominio
        self.funcao  = funcao
        self.derivFuncao = derivFuncao
        self.phi = phi
        self.derivPhi = derivPhi
        
    
    
    def mostrarFuncao(self):
        x = self.dominio
        y = np.array(list(map(self.funcao, self.dominio)))
        plt.plot(x,y)
        plt.axhline(y=0, c = 'k',linewidth = 0.5)
        plt.axvline(x=0,c ='k',linewidth = 0.5)
        plt.show()
    
    
    def tabelarValores(self,a,b,n):
        x = np.linspace(a,b,n)
        y = list(map(self.funcao,x))
        for i in range(len(y)):
            print(f'{i+1}\t x: {x[i]} \t y: {y[i]}')
    
    
    def bisseccao(self,a,b,tol):
        
        if (self.funcao(a))*(self.funcao(b))<0:
            c = (a + b)/2
            iteracoes = 1
            while abs(self.funcao(c)) > tol:
                if (self.funcao(c))*(self.funcao(b))<0:
                    a = c
                else: 
                    b = c
                c = (a + b)/2
                iteracoes += 1
            print()
            print("Método: Bissecção.")
            print("x: ",c,"f(x)",(self.funcao(c)))
            print("nº de iterações: ",iteracoes)
            print()           
        else:
            print("nao ha raizes no intervalo dado")
    
    
    
    def ponto_fixo(self,a,b, eps):
        v = (a+b)/2
        if (np.abs(self.derivPhi(a)) >= np.abs(self.derivPhi(b))) & (self.derivPhi(a) < 1):
            if np.abs(self.phi(a)) < 1:
                print("Método: Ponto Fixo")
                print('phi converge')
        elif np.abs(self.derivPhi(b)) < 1:
            if np.abs(self.phi(a)) < 1:
                print("Método: Ponto Fixo")
                print('phi converge')
            x = self.phi(v)
            iteracoes = 1
            while abs(self.funcao(x)) > eps:
                v = self.phi(x)
                x = self.phi(v)
                
                #print('x:',x,'f(x): ',f(x))
                iteracoes += 1
                if self.funcao(x) < eps:
                    break
            print("x: ",x,"f(x)",self.funcao(x))
            print("nº de iterações: ",iteracoes)
            print()
        else:
            print()
            print('phi não converge')
            print()

    def falsa_posicao(self,a,b,tol):
        if (self.funcao(a))*(self.funcao(b))<0:
            x = (a*self.funcao(b)-b*self.funcao(a))/(self.funcao(b) - self.funcao(a))
            iteracoes = 1
            while abs(self.funcao(x)) > tol:
                if (self.funcao(x))*(self.funcao(b))<0:
                    a = x
                else: 
                    b = x
                x = (a*self.funcao(b)-b*self.funcao(a))/(self.funcao(b) - self.funcao(a))
                iteracoes += 1
            print()
            print("Método:  Falsa Posição.")
            print("x: ",x,"f(x)",(self.funcao(x)))
            print("nº de iterações: ",iteracoes)
            print()
        else:
            print()
            print("nao ha raizes no intervalo dado")
            print()

    def newton_raphson(self,v, eps):
        x = v - (self.funcao(v)/self.derivFuncao(v))
        iteracoes = 1
        while self.funcao(x) > eps:
            x = v - (self.funcao(v)/self.derivFuncao(v))
            v = x
            #print('x:',x,'f(x): ',f(x))
            iteracoes += 1
            if self.funcao(x) < eps:
                break
        print()    
        print("Método: Newton-Raphson")
        print("x: ",x,"f(x)",self.funcao(x))
        print("nº de iterações: ",iteracoes)
        print()

    def secante(self,x_0,x_1, eps):
        iteracoes = 1
        while np.abs(self.funcao(x_1)) > eps:
            x = (x_0*self.funcao(x_1) - x_1*self.funcao(x_0))/(self.funcao(x_1) - self.funcao(x_0))
            x_0 = x_1
            x_1 = x
            iteracoes +=1
        print()    
        print("Método:  Secante.")
        print("x: ",x,"f(x)",(self.funcao(x)))
        print("nº de iterações: ",iteracoes)
        print()


def polyG3(x = None):
    return x**3 - 9*x  +3

def derivPolyG3(x = None):
    return 3*x**2 -9

def phiPolyG3(x = None):
    return (x**3+3)/9
def derivPhiPolyG3(x = None):
    return x**2/3

r = RaizNumerica(dominio = np.linspace(-2,2,1000), funcao=polyG3,derivFuncao=derivPolyG3,phi = phiPolyG3, derivPhi= derivPhiPolyG3)
#r.mostrarFuncao()
#r.bisseccao(a = 0.32, b = 0.34, tol = 10/10**6)
r.ponto_fixo(a = 0.3, b = 0.4, eps = 5/10**4)
r.falsa_posicao(a = 0.32,b = 0.34, tol = 1/10**6)
r.newton_raphson(v = 1, eps = 5/10**4)
r.secante(x_0 = 0.32, x_1 = 0.34,eps = 5/(10**4))

