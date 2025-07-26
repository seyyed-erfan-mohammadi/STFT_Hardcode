# 🧠 EEG Time-Frequency Analysis & Topographic Mapping


This project presents a comprehensive analysis pipeline for EEG time-frequency decomposition and topographical visualization, implemented entirely in MATLAB. Using a sample EEG dataset, it explores and compares three widely used methods for extracting spectral power over time: **Short-Time Fourier Transform (STFT)**, **Complex Morlet Wavelet Convolution**, and **Filter-Hilbert Transform**. Each method is implemented from scratch, offering transparent and educational insight into how time-frequency representations are computed.

The STFT is used to calculate frequency-specific power across all scalp electrodes, enabling the creation of detailed topographical maps. These maps highlight power in the theta band (\~6 Hz) and alpha band (\~10 Hz) at two key time points (150 ms and 700 ms), revealing how brain rhythms evolve in both time and space. The topographical plots are generated using EEGLAB's `topoplot` function to offer intuitive visualizations of the underlying brain dynamics.

To further investigate method-specific differences, the power over time is extracted from the F1 electrode using all three techniques at 10 Hz. The resulting time-series plots demonstrate the effects of temporal smoothing and spectral leakage inherent to each approach. Special care is taken to tune parameters, such as wavelet cycles and STFT window size and overlap, to align the temporal resolution across methods as much as possible.

This code provides a hands-on framework for learning and comparing core concepts in EEG signal processing, and it serves as a visual tool for understanding the strengths and limitations of different time-frequency analysis techniques.

## 📊 Topographical EEG Mapping

We compute the STFT of each EEG electrode and visualize scalp power maps at:

| Frequency Band | Time Points    |
| -------------- | -------------- |
| Theta (∼6 Hz)  | 150 ms, 700 ms |
| Alpha (∼10 Hz) | 150 ms, 700 ms |

Plots are created using EEGLAB's `topoplot` function, revealing spatial distributions of brain oscillations across the scalp.

---

## 🔍 Time-Frequency Method Comparison

At electrode **F1**, we extract power over time using:

* ✅ **STFT** – Short-time windowed FFT with high overlap
* ✅ **Wavelet** – Morlet wavelet convolution (7 cycles)
* ✅ **Hilbert** – FIR filter followed by Hilbert transform

All methods target **10 Hz** activity. We adjust STFT parameters (300 ms window, 98% overlap) to approximate wavelet-like smoothing.
