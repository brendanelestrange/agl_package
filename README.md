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
    conda create -n agl_fast
    conda activate agl_fast
    conda install -c conda-forge chemfiles
    conda install xgboost scikit-learn numpy pandas
    ```

### 3. Use ./make.sh
Ensure you're in the base directory and run the file. It is likely you will have to fiddle with the CMakeLists.txt to make it work.
    ```bash
    chmod +x ./make.sh
    ./make.sh
    ```

### 4. Use it for Feature Generation
Now we can use this the same way we use AGL-Score. You can generate the features using get_agl_features.py in the src directory. Please get the data from PDBBind. Here's a sample usage: 
    ```bash
     python get_agl_features.py -k 112 -c 12 -m Adjacency -f '../csv_data_file/PDBbindv2016_GeneralSet.csv' -dd '../data/2016' -fd '../Features'
    ```
