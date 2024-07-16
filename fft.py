import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def read_csv_file(file_path):
    """CSV 파일에서 데이터를 읽어 NumPy 배열로 반환하는 함수"""
    try:
        data = pd.read_csv(file_path, header=None)
        return data[0].values
    except FileNotFoundError:
        print(f"파일을 찾을 수 없습니다: {file_path}")
        return None
    except pd.errors.EmptyDataError:
        print("빈 CSV 파일입니다.")
        return None
    except pd.errors.ParserError:
        print("CSV 파일을 파싱하는 동안 오류가 발생했습니다.")
        return None

def perform_fft(data):
    """데이터에 대해 FFT를 수행하고 진폭 스펙트럼을 반환하는 함수"""
    fft_result = np.fft.fft(data)
    fft_amplitude = np.abs(fft_result)
    return fft_result, fft_amplitude

def plot_fft(fft_amplitude):
    """FFT 진폭 스펙트럼을 시각화하는 함수"""
    plt.figure(figsize=(10, 6))
    plt.plot(fft_amplitude)
    plt.title('FFT Amplitude Spectrum')
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()

def main():
    # 사용자로부터 CSV 파일 경로 입력 받기
    file_path = input("CSV 파일 경로를 입력하세요: ").strip()
    
    # CSV 파일에서 데이터 읽기
    data = read_csv_file(file_path)
    
    if data is not None:
        # FFT 수행
        fft_result, fft_amplitude = perform_fft(data)
        
        # 결과 출력
        print("FFT 결과 (복소수):")
        print(fft_result)
        print("\nFFT Amplitude:")
        print(fft_amplitude)
        
        # 결과 시각화
        plot_fft(fft_amplitude)

if __name__ == "__main__":
    main()

