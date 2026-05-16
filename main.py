"Objetivo: O objetivo do trabalho é fazer um filtro que filtre um sinal de encoder, ou seja, o livre de ruido \
    "
    
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

"Variaveis globais"




class analogFilt:
    """ Implementação do filtro analógico. Topologia Sallen Key, usando sistema de segunda ordem
    A ideia é utilizar o filtro passa altas de butterworth para filtrar um sinal de  (aceleração), antes que seja efetuada a integração numérica. A entrada do filtro será um arquivo binário de 6k bytes, 2k bytes para cada eixo"""
    def __init__(self,fc, fs, n, signal):
        self.fc = fc
        self.fs = fs
        self.t = np.linspace(0, n, self.fs, endpoint=False)
        self.signal = signal
        self.W0 = 2*np.pi * self.fc
        self.R = 0
        self.C = 10**(-6)
        self.RC = 0
        "PARA BUTTERWORTH"
        self.Q = np.sqrt(2)/2
        self.order = 4
        
        
        "Executa as funções"
        self.fn_calculate()
        self.fn_printAll()
        self.fn_SOS()
        
    def fn_calculate(self):
        self.RC = 1/(self.W0)
        self.R = int(self.RC / self.C)
        return True
    
    def fn_SOS(self):
        sos = signal.butter(self.order, self.fc, btype='high', fs=self.fs, output='sos')
        return sos

    def fn_SOS_TF(self):
        lista = self.fn_SOS()
        for i, order in enumerate(lista):
            a_z = order[3:6]
            b_z = order[0:3] 
            b_s, a_s = signal.bilinear(b_z, a_z, fs=2*self.fs)
            "Normalizo"
            b_s = b_s / a_s[0]
            a_s = a_s / a_s[0]
            self.fn_plot_TF(b_s, a_s, i)
            
    def fn_plot_TF(self,b,a, n_sessao):
        print(f"\n[Seção {n_sessao}] H_{n_sessao}(s) =")
        print(f"     {b[0]:.4e}*s^2 + {b[1]:.4e}*s + {b[2]:.4e}")
        print(f"    --------------------------------------------------")
        print(f"     {a[0]:.4e}*s^2 + {a[1]:.4e}*s + {a[2]:.4e}")
    def fn_printAll(self):
        print("=======================================================")
        print("========================================================")
        print("====RELATÓRIO DE COMPONENTES E LISTA DE MATERIAIS=======")
        print("========================================================")
        print("========================================================")
        print(f"R={self.R} ohm")
        print("========================================================")
        print(f"C={self.C} farad")
        print("========================================================")
        self.fn_SOS_TF()
    
def fn_signal():
    signal = []
    return signal
    
if __name__ == "__main__":
    "Frequencia de corte utilizada é igual a 10hz"
    fc = 10
    "O numero de amostras é igual a 1024"
    n = 1024
    "A frequência de amostragem é igual a 3600Hz, frequencia do adxl345"
    fs = 3600
    "O sinal é adquirido através de uma função"
    filtro = analogFilt(fc=fc,fs=fs,n=n,signal=fn_signal)
   