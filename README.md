# Primer Design App: Two-tailed RT-qPCR for miRNA

A specialized web application for designing **Two-tailed RT primers** for highly specific microRNA (miRNA) quantification. This tool implements the thermodynamic modeling described by **Androvic et al. (2017)** to ensure high sensitivity and specificity.

## 🔬 Scientific Background

MicroRNAs are short (~22 nt) non-coding RNAs. Standard RT-qPCR methods often struggle with specificity due to the short length of the target. This application implements the **Two-tailed RT-qPCR** approach, which offers superior specificity by using a novel RT primer design.

### The Method: Two-tailed RT-qPCR
Unlike standard stem-loop primers, the **Two-tailed RT primer** contains two hemiprobes separated by a hairpin structure:
1.  **3' Hemiprobe**: Binds to the 3' end of the miRNA.
2.  **5' Hemiprobe**: Binds to the 5' end (or a central region) of the miRNA.
3.  **Hairpin Tether**: A stable stem-loop structure connects the two hemiprobes, preventing them from binding to each other or forming non-specific products.

This "clamp" mechanism ensures that the primer only extends when *both* hemiprobes hybridize to the target miRNA, significantly increasing specificity and allowing discrimination of homologous miRNAs with single-nucleotide differences.

## 🧮 Algorithmic Implementation

The application uses rigorous thermodynamic calculations to design the two hemiprobes and the connecting hairpin.

### 1. Thermodynamic Modeling
Melting temperatures ($T_m$) and Gibbs free energy ($\Delta G$) changes are calculated using the **Nearest-Neighbor (NN) Model** with parameters from **SantaLucia et al. (1998)**.
*   **$\Delta G$ (Gibbs Free Energy)**: Used to optimize the binding stability of both the 3' and 5' hemiprobes.
*   **$T_m$ (Melting Temperature)**: Calculated correcting for salt ($Na^+$, $Mg^{2+}$) and oligonucleotide concentrations.

### 2. Design Constraints
*   **Hemiprobe Balance**: The algorithm searches for optimal lengths for both the 3' and 5' hemiprobes to ensure balanced binding stability within a user-defined $\Delta G$ window.
*   **Structure Check**: The hairpin linker is designed to be stable at the RT temperature (42°C) but accessible for the polymerase.
*   **Specificity**: The split-probe design minimizes off-target binding to precursor miRNAs or other small RNAs.

## ✨ Features

*   **Two-tailed Primer Design**: Automated generation of RT primers with dual binding sites.
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

1.  **Androvic, P., Valihrach, L., Elling, J., Sjoback, R., & Kubista, M. (2017).** "Two-tailed RT-qPCR: a novel method for highly accurate miRNA quantification." *Nucleic Acids Research*, 45(15), e144. [https://doi.org/10.1093/nar/gkx588](https://doi.org/10.1093/nar/gkx588)
2.  **SantaLucia, J. Jr. (1998).** "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." *Proceedings of the National Academy of Sciences*, 95(4), 1460-1465.

## ✅ License
Copyright © 2025 Ruslan Klassen. All rights reserved.
