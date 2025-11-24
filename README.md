# Primer Design App: Stem-Loop RT-qPCR for miRNA

A specialized web application for designing **stem-loop pulsed reverse transcription (RT) primers** and specific **PCR primer pairs** for microRNA (miRNA) quantification. This tool implements thermodynamic modeling to ensure high specificity and efficiency in miRNA detection.

## 🔬 Scientific Background

MicroRNAs are short (~22 nt) non-coding RNAs that play critical roles in gene regulation. Due to their short length, standard PCR methods are insufficient. This application designs primers for the **Stem-Loop RT-qPCR** method, widely regarded as the gold standard for miRNA quantification due to its high sensitivity and specificity.

### The Method
1.  **Stem-Loop RT Primer**: A specifically designed hairpin oligonucleotide binds to the 3' end of the mature miRNA. The stem-loop structure stabilizes the interaction and prevents binding to precursor molecules (pre-miRNA/pri-miRNA), enhancing specificity.
2.  **Reverse Transcription**: The miRNA is reverse transcribed into cDNA.
3.  **qPCR**: The cDNA is amplified using a specific **Forward Primer** (matching the miRNA sequence) and a universal or specific **Reverse Primer** (binding to the stem-loop structure).

## 🧮 Algorithmic Implementation

The application uses rigorous thermodynamic calculations to predict primer-template stability and secondary structures.

### 1. Thermodynamic Modeling
Melting temperatures ($T_m$) and Gibbs free energy ($\Delta G$) changes are calculated using the **Nearest-Neighbor (NN) Model** with parameters from **SantaLucia et al. (1998)**.
*   **$\Delta G$ (Gibbs Free Energy)**: Determines the stability of the primer-template duplex. More negative values indicate stronger binding.
*   **$T_m$ (Melting Temperature)**: Calculated correcting for salt ($Na^+$, $Mg^{2+}$) and oligonucleotide concentrations using the entropy ($\Delta S$) and enthalpy ($\Delta H$) of adjacent base pairs.

### 2. Design Constraints & Heuristics
*   **3' End Stability**: The algorithm specifically optimizes the $\Delta G$ of the last 5 nucleotides at the 3' end to prevent non-specific priming (the "GC clamp" effect).
*   **Hairpin Stability**: The stem-loop structure is verified to ensure it forms a stable hairpin at the RT temperature (typically 42°C) but unfolds during the PCR denaturation step (95°C).
*   **Secondary Structure Check**: (Optional) Integration with **ViennaRNA** allows for folding predictions to avoid internal hairpins or self-dimers in the primers themselves.

## ✨ Features

*   **Automated Stem-Loop Design**: Generates RT primers with customizable stem and loop sequences.
*   **Thermodynamic Filtering**: Filters candidates based on $\Delta G$, $T_m$, and GC content windows.
*   **Batch Processing**: Design primers for multiple miRNA sequences simultaneously.
*   **miRNA Database**: Built-in lookup for standard miRBase entries.
*   **Visual Feedback**: Color-coded visualization of the primer-miRNA hybridization complex.

## 🚀 Quickstart

### Local Installation
```bash
git clone https://github.com/Rusklass/primer_design_app.git
cd primer_design_app

# Create virtual environment
python -m venv .venv
# Windows
.venv\Scripts\activate
# Linux/Mac
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run the application
flask run
```
Open [http://127.0.0.1:5000](http://127.0.0.1:5000) in your browser.

### Deployment
This app is configured for easy deployment on **Render**.
1.  Fork this repo.
2.  Create a new **Web Service** on Render.
3.  Connect your repo.
4.  Render will automatically use the `render.yaml` blueprint.

## 📚 References

1.  **SantaLucia, J. Jr. (1998).** "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." *Proceedings of the National Academy of Sciences*, 95(4), 1460-1465.
2.  **Chen, C., et al. (2005).** "Real-time quantification of microRNAs by stem-loop RT-PCR." *Nucleic Acids Research*, 33(20), e179.

## ✅ License
Copyright © 2025 Ruslan Klassen. All rights reserved.
