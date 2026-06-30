# Objetivo: O objetivo do trabalho é fazer um filtro que filtre um sinal de encoder, ou seja, o livre de ruido

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as spy_signal  # Renomeado para evitar conflito com a variável de sinal
from scipy.fft import fft, fftfreq
"Variaveis globais"

class analogFilt:
    """ Implementação do filtro analógico. Topologia Sallen Key, usando sistema de segunda ordem
    A ideia é utilizar o filtro passa altas de butterworth para filtrar um sinal de  (aceleração), antes que seja efetuada a integração numérica. A entrada do filtro será um arquivo binário de 6k bytes, 2k bytes para cada eixo"""
    def __init__(self, fc, fs, n, signal_data):
        self.fc = fc
        self.fs = fs
        
        # CORREÇÃO 1: O vetor de tempo deve ir de 0 até o tempo total (n / fs) com exatamente 'n' pontos
        self.t = np.arange(n) / self.fs
        
        self.signal = signal_data
        self.order = 4
        self.W0 = 2*np.pi * self.fc
        self.R = {"primeiro_estagio":[0,0],
                  "segundo_estagio":[0,0]}
        self.C = 0.00001
        self.n = 1024
        self.RC = np.zeros(int(self.order/2))

        self.filtered_signal = None
        "PARA BUTTERWORTH, 4 ordem"
        self.Q = np.zeros(int(self.order/2))
        
        "Ganhos"
        self.K = np.zeros(int(self.order/2))
        self.K[0] = 1
        self.Rf = np.zeros(int(self.order/2))
        self.Rg = np.zeros(int(self.order/2))
        self.Rf[0] = 0
        self.Rg[0] = 0
        
        
        
        "Executa as funções"
        self.fn_calculate()
        self.fn_SOS()
        self.fn_printAll()
        
        
        self.fn_filter_signal()
        
        self.fn_plot_bode()
        self.fn_plot_signals()
        self.fn_plot_fft()
        
    def fn_calculate(self):
        """ Calcula os componentes do filtro Sallen-Key para Butterworth de 4ª ordem """
        
        # Para Butterworth de 4ª ordem, os fatores Q são:
        # Q1 = 0.5412 (para o primeiro estágio)
        # Q2 = 1.3066 (para o segundo estágio)
        # Estes valores vêm das raízes do polinômio de Butterworth
        
        if self.order == 4:
            # Fatores Q para Butterworth de 4ª ordem
            self.Q[0] = 0.5412  
            self.Q[1] = 1.3066  
            
            "Calculo dos resistores"
            
            self.R["primeiro_estagio"][0] = 1/(2*self.Q[0]*self.W0*self.C)
            self.R["primeiro_estagio"][1] = 2*self.Q[0]/(self.W0*self.C)
            
            
            self.R["segundo_estagio"][0] = 1/(2*self.Q[1]*self.W0*self.C)
            self.R["segundo_estagio"][1] = 2*self.Q[1]/(self.W0*self.C)
    
    def fn_SOS(self):
        sos = spy_signal.butter(self.order, self.fc, btype='high', fs=self.fs, output='sos')
        return sos

    def fn_SOS_TF(self):
        lista = self.fn_SOS()
        for i, order in enumerate(lista):
            a_z = order[3:6]
            b_z = order[0:3] 
            b_s, a_s = spy_signal.bilinear(b_z, a_z, fs=2*self.fs)
            "Normalizo"
            b_s = b_s / a_s[0]
            a_s = a_s / a_s[0]
            "Posso pegar os valores de resitores a partir daqui"
            self.fn_plot_TF(b_s, a_s, i)

    
        
    def fn_filter_signal(self):
        """ Passa o sinal de entrada pela matriz de seções de segunda ordem (SOS) """
        sos = self.fn_SOS()
        # Usando sosfiltfilt para evitar qualquer atraso de fase na integração posterior
        self.filtered_signal = spy_signal.sosfiltfilt(sos, self.signal)
        
    
        
    def fn_plot_TF(self, b, a, n_sessao):
        print(f"\n[Seção {n_sessao+1}] H_{n_sessao+1}(s) =")
        print(f"     {b[0]:.4e}*s^2 + {b[1]:.4e}*s + {b[2]:.4e}")
        print(f"-------------------------------------------------------")
        print(f"     {a[0]:.4e}*s^2 + {a[1]:.4e}*s + {a[2]:.4e}")
        print(f"-------------------------------------------------------")
        
    def fn_plot_bode(self):
        """ Calcula e plota o Diagrama de Bode em escala SEMILOG do filtro completo """
        sos = self.fn_SOS()
        
        # Calcula a resposta em frequência do sistema SOS completo
        w, h = spy_signal.sosfreqz(sos, worN=5000, fs=self.fs)
        
        # Converte a magnitude para dB e a fase para graus
        magnitude_db = 20 * np.log10(np.maximum(abs(h), 1e-5))
        fase_graus = np.degrees(np.unwrap(np.angle(h)))
        
        # Criando a figura do Diagrama de Bode
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
        
        # --- Sub-bloco 1: Magnitude (Semilogx) ---
        ax1.semilogx(w, magnitude_db, color='blue', linewidth=2, label='Filtro Passa-Altas')
        ax1.axvline(self.fc, color='red', linestyle='--', label=f'Corte: {self.fc} Hz')
        ax1.axhline(-3, color='black', linestyle=':', alpha=0.7, label='-3 dB')
        ax1.set_title(f'Diagrama de Bode (Semilog) - Butterworth Ordem {self.order}')
        ax1.set_ylabel('Magnitude (dB)')
        ax1.grid(True, which='both', linestyle='--', alpha=0.5)
        ax1.legend(loc='lower right')
        ax1.set_ylim([-60, 5])
        
        # --- Sub-bloco 2: Fase (Semilogx) ---
        ax2.semilogx(w, fase_graus, color='purple', linewidth=2)
        ax2.axvline(self.fc, color='red', linestyle='--')
        ax2.set_xlabel('Frequência (Hz)')
        ax2.set_ylabel('Fase (Graus)')
        ax2.grid(True, which='both', linestyle='--', alpha=0.5)
        
        plt.xlim([1, self.fs / 2])
        plt.tight_layout()
        plt.show()
        
    def fn_plot_signals(self):
        """ Plota o gráfico comparativo do sinal no domínio do tempo """
        plt.figure(figsize=(10, 4))
        plt.plot(self.t, self.signal, label='Sinal Bruto (ADXL345)', color='gray', alpha=0.7)
        plt.plot(self.t, self.filtered_signal, label='Sinal Filtrado (Passa-Altas)', color='green', linewidth=2)
        
        plt.title('Comparativo do Sinal no Domínio do Tempo')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Aceleração / Amplitude')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plt.show()
        
    def fn_plot_fft(self):
        """ Calcula e plota a FFT do sinal bruto e do sinal filtrado """
        # 1. Calcular a transformada e as frequências correspondentes
        # Pegamos apenas a primeira metade dos pontos (frequências positivas)
        metade = self.n // 2
        freqs = fftfreq(self.n, 1/self.fs)[:metade]
        
        # FFT do sinal bruto (normalizada pelo número de pontos)
        fft_bruto = fft(self.signal)[:metade]
        magnitude_bruta = (2.0 / self.n) * np.abs(fft_bruto)
        
        # FFT do sinal filtrado (normalizada pelo número de pontos)
        fft_filtrado = fft(self.filtered_signal)[:metade]
        magnitude_filtrada = (2.0 / self.n) * np.abs(fft_filtrado)
        
        # 2. Plotar os espectros de frequência
        plt.figure(figsize=(10, 5))
        
        plt.plot(freqs, magnitude_bruta, label='Espectro Bruto', color='gray', alpha=0.7)
        plt.plot(freqs, magnitude_filtrada, label='Espectro Filtrado (Passa-Altas)', color='green', linewidth=1.8)
        
        # Linha vertical indicando onde o filtro de 10Hz atuou
        plt.axvline(self.fc, color='red', linestyle='--', label=f'Corte: {self.fc} Hz')
        
        plt.title('Transformada Rápida de Fourier (FFT)')
        plt.xlabel('Frequência (Hz)')
        plt.ylabel('Amplitude (Magnitude)')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()
        
        # Foco inicial na banda de interesse (de 0 a 300Hz, por exemplo) para ver melhor o corte de 10Hz
        # Se quiser ver até Nyquist (1800Hz), mude para [0, self.fs/2]
        plt.xlim([0, 300]) 
        
        plt.tight_layout()
        plt.show()    
        
    def fn_printAll(self):
        print("=======================================================")
        print("========================================================")
        print("====RELATÓRIO DE COMPONENTES E LISTA DE MATERIAIS=======")
        print("========================================================")
        k = 0
        for stage in ["primeiro_estagio", "segundo_estagio"]:
            print("===================================================")
            print(stage)
            print("===================================================")
            for j in [0,1]:
                print("===================================================")
                print(f"{j} resistor = {self.R[stage][j]} ohms")
                print("===================================================")
            print(f"FATOR DE QUALIDADE = {self.Q[k]}")
            print("===================================================")
            k += 1
        print(f"VALOR DO CAPACITOR = {self.C}")       
        print("===================================================")
        self.fn_SOS_TF()
        
    
    
def fn_signal():
    return []

def fn_signal_binario(caminho_arquivo="dados_sensor.bin"):
    """ Extração dos dados contidos no arquivo binário de 6KB (3 eixos x 1024 amostras int16) """
    try:
        dados_puros = np.fromfile(caminho_arquivo, dtype=np.int16)
        dados_eixos = dados_puros.reshape(3, 1024)
        return dados_eixos[0, :], dados_eixos[1, :], dados_eixos[2, :]
    except FileNotFoundError:
        print(f"Arquivo {caminho_arquivo} não encontrado. Gerando sinal de teste simulado.")
        t_sim = np.arange(1024) / 3200
        sinal_teste = 5.0 * np.sin(2 * np.pi * 2 * t_sim) + 1.5 * np.sin(2 * np.pi * 30 * t_sim)
        return sinal_teste, sinal_teste, sinal_teste

if __name__ == "__main__":
    fc = 10
    n = 1024
    fs = 3200
    
    # Busca os dados reais do arquivo gerado pelo ventilador
    sinal_x, sinal_y, sinal_z = fn_signal_binario("dados_Motor_grande_20260408_211131.bin")
    
    # CORREÇÃO 3: Passando 'sinal_x' (o array numérico de 1024 posições) e mudando o nome do parâmetro para signal_data
    filtro = analogFilt(fc=fc, fs=fs, n=n, signal_data=sinal_x)