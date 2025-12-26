# AGL-Score (Accelerated)

**Fast, C++ accelerated calculation of Multiscale Weighted Colored Subgraphs (MWCGS) for protein-ligand complexes.**

This package calculates algebraic graph features (AGL-Score) for molecular binding interfaces. It uses a hybrid C++/Python architecture (via Pybind11 and Eigen) to achieve orders of magnitude speedup over pure Python implementations.

---

## 🛠 Prerequisites

Before installing the Python package, you **must** have the following C++ libraries installed on your system.

### 1. OpenMP (Required for Parallelism)
* **macOS:**
    ```bash
    brew install libomp
    ```
* **Ubuntu/Debian:**
    ```bash
    sudo apt-get install libomp-dev
    ```

### 2. Chemfiles (Required for PDB/Mol2 Parsing)
The easiest way to install Chemfiles is via Conda/Mamba, which places the headers and libraries where the build system can find them.
```bash
conda install -c conda-forge chemfiles
conda install xgboost scikit-learn numpy pandas