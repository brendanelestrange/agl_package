# AGL-Score (Accelerated)

**Fast, C++ accelerated calculation of Multiscale Weighted Colored Subgraphs (MWCGS) for protein-ligand complexes.**

This package calculates algebraic graph features (AGL-Score) for molecular binding interfaces. It utilizes a hybrid C++/Python architecture (via Pybind11 and Eigen) to achieve orders of magnitude speedup compared to pure Python implementations.

---

## 🛠 Prerequisites

Before building the package, you **must** ensure the following C++ libraries and environments are set up.

### 1. System Dependencies (OpenMP)
OpenMP is required for parallelization support in the C++ layer.

* **macOS:**
    ```bash
    brew install libomp
    ```
* **Ubuntu/Debian:**
    ```bash
    sudo apt-get install libomp-dev
    ```

### 2. Conda Environment (Chemfiles & Python Libs)
The easiest way to handle the dependency on `Chemfiles` (required for PDB/Mol2 parsing) is via Conda/Mamba.

```bash
# Create and activate the environment
conda create -n agl_fast
conda activate agl_fast

# Install C++ dependencies
conda install -c conda-forge chemfiles

# Install Python dependencies
conda install xgboost scikit-learn numpy pandas
