# Nucleon Decay SMEFT (ΔB=1)

This repository contains the code for the RG evolution of Baryon number violating (ΔB=1) Nucleon Decay in the Standard Model Effective Field Theory (SMEFT) framework. It uses Python and `scipy`'s ODE solver to evolve parameters across energy scales. If you use this code, please cite the paper: [arXiv:2511.06106](https://arxiv.org/abs/2511.06106).

## 📁 Project Structure
The main code is in the folder `d6py`. Various example are provided in the root directory as a notebook file to demonstrate how to use the package for different scenarios. The code for nucleon decay is provided in the three notebooks (.ipynb) for different cases as was given in the paper: [arXiv:2511.06106](https://arxiv.org/abs/2511.06106).
```
NucleonDecaySMEFT/
├── beta_functions.py      # Contains the different beta functions
├── d6py                 # Main package directory
│   ├── __init__.py        # Initializes the package
│   ├── BetaSMEFT.txt    # The default beta functions
│   ├── d6py.py         # Main module for RG evolution
│   ├── default_initial_values.py # Default initial values for Wilson coefficients
│   ├── Redundant_Elements.txt # List of Redundant operators
├── results/                # Directory to store results
├── README.md               # Project documentation
├── requirements.txt        # Required Python packages
├── SM_gauge_coupling.ipynb # Example notebook for evolving SM gauge couplings
├── Untitled-1.ipynb     # Example notebook for Nucleon Decay (Case 1)
├── Untitled-2.ipynb     # Example notebook for Nucleon Decay (Case 2)
└── Untitled-3.ipynb     # Example notebook for Nucleon Decay (Case 3)
```

## ⚙️ Features

- Solves renormalization group equations (RGEs) for:  
  - Yukawa couplings: `yu`, `yd`  
  - Gauge couplings: `g1`, `g2`, `g3`  
  - Wilson coefficients  
- Customizable Beta functions via text files

---

## 🚀 Getting Started
**Clone the repository:**
```bash
git clone https://github.com/yourusername/NucleonDecaySMEFT.git
cd NucleonDecaySMEFT
```

### 🛠️ Requirements

Install all dependencies with:

```bash
pip install -r requirements.txt
```

## 🚀 How to Use

A full explanation of all modules and functions is available in the [Wiki](https://github.com/rp-winter/Nucleon-Decay-SMEFT/wiki).

For starters, see the example notebook `SM_gauge_coupling.ipynb` to understand how to evolve the SM gauge couplings. The other notebooks `Untitled-1.ipynb`, `Untitled-2.ipynb`, and `Untitled-3.ipynb` demonstrate the nucleon decay scenarios for the different cases discussed in the paper [arXiv:2511.06106](https://arxiv.org/abs/2511.06106).