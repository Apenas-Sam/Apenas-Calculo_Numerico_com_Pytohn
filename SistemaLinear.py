import numpy as np

class SistemaLinear:
    def __init__(self, matriz = None , vetor_constante = None, dimensao = None, tipo = 'C'):
        self.matriz = matriz #Recebe a matriz dos coeficientes 
        self.vetor_constante = vetor_constante # Recebe o vetor de constantates, tambem conhecido como vetor b
        self.dimensao = matriz.shape #Recebe uma tupla-ordenada com as dimensões do sistema
        self.tipo = tipo #Rece o tipo de sistema, se a instancia possui a matriz de coeficientes e o vetor constante (C), ou se estao numa mesma matriz ,que recebe o nome de matriz extendida (E)
        if np.any(vetor_constante) == None and tipo != 'E': # Caso nao seja informado o vetor constante é criado um vetor com zeros no lugar
            self.vetor_constante = np.zeros(shape=(self.matriz.shape[0],1))
        if tipo == 'E':
            self.matriz = matriz[:,:-1] #Recebe a matriz dos coeficientes 
            self.vetor_constante = matriz[:,-1].reshape((dimensao[0],1)) # Recebe o vetor de constantates, tambem conhecido como vetor b    
        
        
            

    def cramer(self):    
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Método: Cramer <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
        print()
        print(f'Sistema:')
        print(f'{self.matriz}') #Exibe a matriz de coeficientes do sistema
        print()
        if (np.linalg.det(self.matriz) != 0):
            matriz_transicao = np.copy(self.matriz) #Esta matriz_transição recebe uma copia da matriz original para nao alterar seus dados
            dets = [np.linalg.det(matriz_transicao)] #Armazena o valor das determinantes que são necessarias para a reolução do sistema. A determinante da matriz dos coeficientes já está adicionada
            vetor_incognitas = np.zeros_like(self.vetor_constante) #Criamos um vetor com as dimensões do vetor constante
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Determinantes <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
            for i in range(0,self.dimensao[0],1): #Loop responsável pelo calculo e armazenamento das determinantes
                matriz_transicao[:,i] = self.vetor_constante.reshape((self.dimensao[0],)) #Substitui a coluna da matriz pelo vetor constante na coluna de x_i 
                dets.append(np.linalg.det(matriz_transicao)) #Armazena o calculo da determinante
                print()
                print(f'Dx[{i+1}]: {dets[i+1]}')#Exibe o resultado do determinante calculado
                print()
                matriz_transicao[:,i] = np.copy(self.matriz[:,i]) # Retorna o valor original da culuna subistituida pelo vetor constante
            print()
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Soluções <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
            print()  
            for s in range(0,self.dimensao[0],1): #Loop responsavel por armazenar e calcular os valores de x
                vetor_incognitas[s] = dets[s+1]/dets[0] #Armazena os resultados de x
                print()
                print(f'x[{s+1}]: {vetor_incognitas[s]}')#Exibe os valroes de x
            print()
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim do Método: Cramer <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
            return vetor_incognitas #retorna os valores de x, podendo ser utilizados ou armazenados em uma variavel
        else:
            return print("determinante é igual a 0. Não há solução para o sistema, utilizando este método.")

    def triangular_matriz(self):
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Triangulando Sistema <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
        print()
        print(f'Sistema:')
        print(f'{self.matriz}') #Exibe a matriz dos coeficientes 
        print()
        print(f'Vetor Solução:')
        print(f'{self.vetor_constante}') #Exibe o vetor constante
        print()
        for k in range(0,self.matriz.shape[0]-1,1): #Loop responsável pelas iterações necessarias para trinagular a matriz
            #print(f'k: {k + 1}')
            for i in range(k + 1,self.matriz.shape[0],1): #Loop responsável pelos calculos das tranformações da matriz
                m = self.matriz[i,k]/self.matriz[k,k] #constante utilizada na multiplicação da linha anterior a que será subitraida. Exemplo L_2 = L_2 - m*L_1
                self.matriz[i,k] = 0 #Zeramos o elemento abaixo da diagonal principal
                #print(f'i: {i+1}')
                for j in range(k + 1,self.matriz.shape[0],1):
                    self.matriz[i,j] = self.matriz[i,j] - (m*(self.matriz[k,j])) #Realiza um proceso de operação de um elemento da matriz com outro que se encontra uma linha abaixo. EX: a_12 = a_12 - m*a_22
                    #print(f'j: {j+1}')
                print(f'Triangulando.....')
                print(f'{self.matriz}') #Exibe a matriz triangulada, dependendo do passo, não estará completamente triangulada.
                print()
                self.vetor_constante[i,0] = self.vetor_constante[i,0] - (m*(self.vetor_constante[k,0]))  #Realiza um proceso de operação de um elemento do vetor constante com outro que se encontra uma linha abaixo. EX: b_2 = b_2 - m*b_1          
                print(f'Vetor Solução:')
                print(f'{self.vetor_constante}')#Exibe o vetor constante
                print()
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim da Triangulação <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
    
    def eliminacao_gaussiana(self):
        self.triangular_matriz() #É necessário triangular o sistema para utilizar esse método
        print()
        print()
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Método: Eliminação Gaussiana <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
        print()
        print(f'Sistema:')
        print(f'{self.matriz}')#Exibe a matriz
        print()
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Soluções <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
        print()
        x = np.zeros_like(self.vetor_constante)#Vetor incógnita, ou vetor x
        x[-1 ,0] = self.vetor_constante[-1,0]/self.matriz[-1,-1] #O ultimao valor é xn = b_n/ann, a partir disso é possivel substituir o valor de xn nas outras linhas para obter x_(n-1) até x_1
        for i in range(self.dimensao[0]-1,-1,-1): #Loop responsavel por armazenar x_(n-1) até x_1. Como x_1 é na verdade x_0, range possui intervalo aberto [a,b[ então é preciso ir até -1 para i ser igual a 0 e preencher x  
            soma = 0
            for j in range(i+1,self.dimensao[0],1): #Loop responsavel por realizar os calculos de x
                soma += self.matriz[i,j]*x[j,0] #Somando a_ij*x_j de i entre [n-1,-1[ com j = i +1 entre [i+1,n]
            x[i] = (self.vetor_constante[i] - soma)/matriz[i,i] #Calulo de x, x_(n-1) = (b_(n-1) - soma)/a_(n-1,n-1)
        print('vetor x:')
        print(f'{x}') #Exibindo x
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim do Método <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<') 
        return x   

    def convergencia_linhas(self):
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Teste de convergencia: Critéro de linha <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
        print()
        alfas = np.zeros_like(self.vetor_constante)#Vetor que armazena os valores de alfa     
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Alfas <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
        print()
        for i in range(self.dimensao[0]):#Loop para armazenar os valores de alfa
            soma = 0
            for j in range(self.dimensao[0]-1):#Loop responsavel pela soma dos elementos da linha i                               
                if  (i==0):
                    soma += np.abs(self.matriz[i,j+1])#Caso a estejamos na peimeira linha não poderemos utilizar seu primeira elemento, pois ele faz parte da diagonal principal. Entao somamos a partir do proximo elemento    
                elif (i == j):
                    soma += np.abs(self.matriz[i,j+1])#Caso o nº da linha coincidir com o nº da coluna entao, ou seja, caso seja elemento da diagonal principal, somamos o proximo elemento da linha
                else:
                    soma += np.abs(self.matriz[i,j])#Para qualquer outro caso realize a soma do elementos da linha
            alfas[i,0] = soma/np.abs(self.matriz[i,i])#Armazena o valor de alfa na linha i até n
            print(f'alfa[{i+1}]: {alfas[i,0]}')#Exibe os valores do vetor alfa        
        print()
        max_alpha = np.max(alfas) #Armazeno o maior valor de alfa em módulo
        print(f'Alfa máximo: {max_alpha}') #Exibe o maior valor de alfa em módulo
        if max_alpha < 1: #Se o maior valor de alfa em módulo for menor do que 1 ele converge, e então, o sistema pode ser resolvido com o método de Gauss-Jacobi
            print("O método de Gauss-Jacobi pode ser usado neste Sistema Linear")
            return True
        else: #Caso dele nao convergir o método de Gauss-Jacobi nao convirgirá para as soluções desse sistema
            print()
            print("O método não converge para este sistema")
            print("O método de Gauss-Jacobi pode ser usado neste Sistema Linear")
            print()
    
    def gauss_jacobi(self,tol = 10**(-3),iteracoes=30):
        if self.convergencia_linhas() == True: #Caso o sistema convirja suas soluções podem ser calculadas a partir desse método iterativo
            iteracao = 1
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Método: Gauss-Jacobi <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
            print()
            print(f'Sistema:')
            print(f'{self.matriz}')# Exibe o sistema
            print()
            x = np.zeros_like(self.vetor_constante)#Irá receber as aproxiamações do método
            for i in range(0,self.dimensao[0],1):#Loop responsável por nos dar o valor de um 'chute'inicial para x que segue a regra: x_i = b_i/a_ii
                x[i,0] = self.vetor_constante[i,0]/self.matriz[i,i]         
            x_next = np.zeros(shape = x.shape)#Irá receber os valores das novas aproximações de x
            while iteracao <= iteracoes: #Este método será iterado até um valor limite de iterações que pode ser alterado ná hora de utilizar o método(da classe Sistema Linear), ou utilizar o limite pradrao do método
                print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> {iteracao}ª iteração.  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                print()
                print(f'Nº máximo de iterações {iteracoes}') #Exibe o valor máximo de iterações
                print(f'{iteracao}ª iteração de x') #Exibe em qual iteração de x estamos
                print(f'{x}') #Exibe o valor de x      
                print()
                print()
                print(f'Nova aproximação:')
                print()
                for i in range(self.matriz.shape[0]): #Loop responsável por atribuir os valores para nossa proxima aproxima para x, que se chama: x_next
                    soma = self.vetor_constante[i,0] #Primeiro é adicionado o valor de b_i que será utilizado para a proxima apraximação de x^(k)_i, ou seja x^(k+1)
                    for j in range(self.matriz.shape[1]):#Loop responsável por realizar o somatório
                        if (i != j): #Caso o elemento da matriz dos coeficientes não pertença a diagonal principal realize a soma
                            soma -= self.matriz[i,j]*x[j,0]#Ele decrementa o valor da soma a_ij*x^(k)_(j), pois esta variavél soma é na verdade a seguinte espressão: b_i - Somatorio de a_ij*x^(k)_j quando i é diferente de j
                    x_next[i,0] = soma/self.matriz[i,i]#É adicionado a próxima valor de x^(k+1)_j que é (b_i - Somatorio de a_ij*x^(k)_j quando i é diferente de j) dividido por a_ii
                    
                    print(f'x_next[{i+1}] = {x_next[i,0]}')#exibe o vetor de x^(k+1)_i
                erro = np.max(np.abs(x_next-x)) #Calcula o erro em módolo entre x^(k+1) e x^(k) e seleciona o maior valor
                print()
                print(f'Dados sobre a iteração ')
                print(f'erro: {erro}') #Exibe o erro
                print()
                print(f'x_next:')
                print(x_next)##exibe o vetor de x^(k+1)
                print()
                if erro < tol: #Caso o erro seja menor do que a tolerância, temos a solução do sistema sendo x^(k + 1)
                    print()
                    print('Solunção final')
                    print(x_next)#Exibe x^(k+1)
                    print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim desta iteração <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                    print()
                    print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim do Método <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                    break
                else:#Caso o erro seja maior do que a tolerância, temos que x assumirá os valores de  x^(k+1), e assim até ele satisfazer a condição de parada ou atingir o limite de iterações
                    x = np.copy(x_next) #x recebe os valores de x^(k+1) desta iteração
                    print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim desta iteração <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                    print()
                iteracao += 1 #Acrescentar mais um ao nosso iterador
            else: #Caso o sistema não atenda a convegencia de linhas
                print()
                print("O método não converge para este sistema")
                print("O método de Gauss-Jacobi pode ser usado neste Sistema Linear")
                print()
                


    def convergencia_sassenfeld(self):
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Teste de convergencia: Critério de Sassenfeld <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
        print()
        betas = np.zeros_like(self.vetor_constante)#Armazena os valores de betas
        print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Betas <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')    
        print()
        for i in range(self.dimensao[0]):#Loop responsavel por adicionar os valores de beta_i no vetor betas
            soma = 0
            for j in range(self.dimensao[1]):#Loop responsavel pelo calculo dos betas 
                if (i != j ) and (i == 0) or (i<j):#Caso i seja diferente de j e a_ij com i e j são iguais ao primeiro elemento da diagonal principal ou i menor que j faça:  
                    soma += (np.abs(self.matriz[i,j])) #Somando |a_ij|
                elif (i!=j) and (i != 0): 
                    soma += (np.abs(self.matriz[i,j]))*betas[j]#Somando |a_ij|*beta_j
            betas[i,0] = soma/(np.abs(self.matriz[i,i]))#Adicionando b_i no vetor betas na posição i
            print(f'beta[{i+1}]: {betas[i,0]}')#Exibindo beta_i        
        print()
        max_beta = np.max(betas)
        print(f'Beta máximo: {max_beta}')#Exibe o beta máximo
        print()
        if max_beta < 1: #Caso o beta máximo seja menor que 1 então o sistema converge 
            print("O método de Gauss-Seidel pode ser usado neste Sistema Linear")
            print()
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim do Método <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
            return True 
        else:
            print()
            print("O método não converge para este sistema")
            print("O método de Gauss-Seidel pode ser usado neste Sistema Linear")
            print()
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim do Método <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
            print()
    
    def gauss_seidel(self,iteracoes=30,tol = 10**(-3)):
        if self.convergencia_sassenfeld() == True: #Caso o sistema convirja então a solução do sistema pode ser calculada com este método 
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Método: Gauss-Seidel <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
            print()
            print(f'Sistema:')
            print(f'{self.matriz}')#Exibe a Matriz
            print()
            iteracao = 1 #Setando nosso iterador para 1, pois este será usado no laço while  
            x = np.zeros_like(self.vetor_constante)#Irá receber as aproxiamações do método
            for i in range(0,self.matriz.shape[0],1):#Loop responsável por nos dar o valor de um 'chute'inicial para x que segue a regra: x_i = b_i/a_ii
                x[i,0] = self.vetor_constante[i,0]/self.matriz[i,i] #x_i = b_i/a_ii
            x_next = np.zeros_like(x)#Irá receber os valores das novas aproximações de x
            
            while iteracao <= iteracoes: #Este método será iterado até um valor limite de iterações que pode ser alterado ná hora de utilizar o método(da classe Sistema Linear), ou utilizar o limite pradrao do método      
                print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> {iteracao}ª iteração.  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                print()
                print(f'Nº máximo de iterações {iteracoes}') #Exibe o valor máximo de iterações
                print(f'{iteracao}ª iteração de x')#Exibe em qual iteração de x estamos
                print(f'{x}') #Exibe o valor de x       
                print()
                print()
                print(f'Nova aproximação:')
                print()            
                for i in range(self.matriz.shape[0]):#Loop responsavel por adicionar x^(k+1)_i
                    soma = self.vetor_constante[i,0] #Primeiro é adicionado o valor de b_i que será utilizado para a proxima apraximação de x^(k)_i, ou seja x^(k+1)
                    for j in range(self.matriz.shape[1]): #Loop responsável por realizar os somatórios
                        if (i > j): #Caso i seja maior que j será decrementado a_ij*x^(k+1) 
                            soma -= self.matriz[i,j]*x_next[j,0] #Entende-se por - somatoriio de i > j de a_ij com a aproximação de x^(k+1)_j, caso já haja um x^(k+1)_j nesta iteração
                        elif (i < j): #Caso i seja menor que j será decrementado a_ij*x^(k) 
                            soma -= self.matriz[i,j]*x[j,0]#Entende-se por - somatoriio de i > j de a_ij com a aproximação de x^(k)_j, para ques seja utilizado x^(k+1)_j na próxima iteração
                    x_next[i,0] = soma/self.matriz[i,i] # entende-se por (b_i  - somatoriio de i > j de a_ij com a aproximação de x^(k+1)_j - somatoriio de i > j de a_ij com a aproximação de x^(k)_j)/a_ii
                    print(f'x_next[{i+1}] = {x_next[i,0]}')#exibe o vetor de x^(k+1)_i            
                erro = np.max(np.abs(x_next - x))#Calcula o erro em módolo entre x^(k+1) e x^(k) e seleciona o maior valor
                print()
                print(f'Dados sobre a iteração ')
                print(f'erro: {erro}') 
                print()
                print(f'x_next:')
                print(x_next)
                print()
                if erro > tol:#Caso o erro seja maior do que a tolerância, temos que x assumirá os valores de  x^(k+1), e assim até ele satisfazer a condição de parada ou atingir o limite de iterações
                    x = np.copy(x_next)#x recebe os valores de x^(k+1) desta iteração
                    print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim desta iteração <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                    print()
                    iteracao += 1
                else:
                    print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Solução <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                    print(x_next)#Exibe x^(k+1)
                    print()
                    print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Fim do Método <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                    break
           
