from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Definindo a Classe Foguete


class Foguete:
 # Construtor da Classe
    def __init__(self, mi, cd, nome):            # Parâmetros de Entrada do Construtor
        self.mi = mi                             # Massa inicial [Kg]
        self.cd = cd                             # Coeficiente de Arrasto
        self.nome = nome                         # Nome do arquivo de empuxo
        self.importarArquivo()                   # Assim que a classe for criada, o arquivo com a curva já vai ser importado

    def importarArquivo(self):                   # Método
        tabelaEmpuxoArquivo = open(self.nome)    # Importando o arquivo da curva

        self.tabelaTempo = []                    # Criando a tabela que vai receber os valores da curva
        self.tabelaEmpuxo = []                   # Criando a tabela que vai receber os valores da curva

        for linha in tabelaEmpuxoArquivo:        # Percorrendo o arquivo e adicionando a um dicionário
            tempo, empuxo = linha.split()
            self.tabelaEmpuxo[float(tempo)] = float(empuxo)  # Para o tempo dado, estamos pegando o empuxo
            self.tabelaTempo.append(float(tempo))            # Para cada linha do arquivo, a primeira informação é armazenada na tabela do de tempo
            self.tabelaEmpuxo.append(float(empuxo))          # Para cada linha do arquivo, a segunda informação é armazenada na tabela do de empuxo
        
    def thrust(self, tempo):
        print(self.tabelaEmpuxo[tempo])


    def plotEmpuxo(self):
        plt.plot(self.tabelaTempo,self.tabelaEmpuxo)
        plt.axis([0, 2, 0, 600])
        plt.show()

    def interpolador(self):
        listaTempo = np.arange(0,1.8923,step=0.1)

        for tempo in listaTempo: 
            cubic_interp = interp1d(listaTempo,self.tabelaEmpuxo, kind="cubic")
            self.cubic_results = cubic_interp(tempo)
            
        print(cubic_interp)


