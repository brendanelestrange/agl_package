#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <chemfiles.hpp>
#include <Eigen/Dense>
#include <omp.h>

namespace py = pybind11;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::string;

const std::vector<std::string> PROTEIN_TYPES = {
    "C", "N", "O", "S"
};

const std::vector<std::string> LIGAND_ELEMS = {
    "H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"
};

std::unordered_map<std::string, float> vwals_rad = {
    {"H",  1.20},
    {"C",  1.70},
    {"N",  1.55},
    {"O",  1.52},
    {"F",  1.47},
    {"P",  1.80},
    {"S",  1.80},
    {"Cl", 1.75},
    {"Br", 1.85},
    {"I",  1.98}
};

const float SIGMA = 0.2127; 

struct Point { float x, y, z; };

struct trajFile {
    chemfiles::Frame frame;
    chemfiles::Topology topology;
    auto get_positions() const { return frame.positions(); }
};

// --- HELPERS ---

float get_radius(const string& atom_type) {
    if (vwals_rad.count(atom_type)) return vwals_rad.at(atom_type);
    string first_char = atom_type.substr(0, 1);
    if (vwals_rad.count(first_char)) return vwals_rad.at(first_char);
    return 1.70f; 
}

VectorXd filter_pos(const VectorXd &values) {
    std::vector<double> temp_list;
    for (int i = 0; i < values.size(); i++) {
        if (values[i] > 1e-5) temp_list.push_back(values[i]);
    }
    return Eigen::Map<VectorXd>(temp_list.data(), temp_list.size());
}

float distance_calc(const Point& first, const Point& second) {
    double dx = first.x - second.x;
    double dy = first.y - second.y;
    double dz = first.z - second.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

float kernel(float dist, float vdw_rad, float kappa, float tau, const string& k_type) {
    float eta = tau * vdw_rad;
    char type_char = std::tolower(k_type[0]); // Check first char ('e', 'l', 'r', etc.)

    // 1. Exponential
    if (type_char == 'e') {
        return std::exp( -std::pow(dist / eta, kappa) );
    }
    // 2. Lorentz
    else if (type_char == 'l') {
        return 1.0f / (1.0f + std::pow(dist / eta, kappa));
    }
    // 3. Radial Basis (RBF)
    else if (type_char == 'r') {
        return std::exp( -std::pow(dist, 2) / (2.0f * std::pow(eta, 2)) );
    }
    // 4. Sigmoid
    else if (type_char == 's') {
        // 1 / tanh(kappa * d - eta)
        return 1.0f / std::tanh(kappa * dist - eta);
    }
    // 5. Polynomial
    else if (type_char == 'p') {
        float ratio = dist / eta;
        // (1 + 2x^3 - 3x^2)^kappa
        float poly = 1.0f + 2.0f * std::pow(ratio, 3) - 3.0f * std::pow(ratio, 2);
        return std::pow(poly, kappa);
    }

    // Default to Exponential if unknown
    return std::exp( -std::pow(dist / eta, kappa) );
}

string trim(const string& str) {
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first) return str;
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

std::map<string, vector<int>> get_atom_groups(const trajFile& file, bool is_ligand) {
    std::map<string, vector<int>> groups;
    for (size_t i = 0; i < file.frame.size(); i++) {
        string key;
        if (is_ligand) {
            string type = file.topology[i].type();
            size_t dot_pos = type.find('.');
            key = (dot_pos != string::npos) ? type.substr(0, dot_pos) : type;
        } else {
            key = trim(file.topology[i].name());
        }
        groups[key].push_back(i);
    }
    return groups;
}

vector<float> calculate_stats(const VectorXd &raw_eigens, float counts) {
    VectorXd eigens = filter_pos(raw_eigens);
    
    if (eigens.size() == 0) {
        vector<float> empty_vec(10, 0.0f);
        empty_vec[0] = counts; 
        return empty_vec;
    }

    float mean = eigens.mean();
    float sum = eigens.sum();
    float variance = (eigens.array() - mean).square().sum() / (eigens.size()); 
    float std_dev = std::sqrt(variance);

    float median_val;
    int n = eigens.size();
    VectorXd sorted = eigens;
    std::sort(sorted.data(), sorted.data() + sorted.size());

    if (n % 2 == 1) median_val = sorted(n / 2);
    else median_val = (sorted(n / 2 - 1) + sorted(n / 2)) / 2.0;

    return {
        counts,
        sum,
        (float)eigens.minCoeff(),
        (float)eigens.maxCoeff(),
        mean,
        median_val,
        std_dev,
        variance,
        (float)eigens.size(),
        (float)eigens.squaredNorm()
    };
}

// --- ANALYZER ---
vector<float> analyze_pair(const trajFile &pro_file, const vector<int>& pro_indices,
                           const trajFile &lig_file, const vector<int>& lig_indices,
                           string p_type, string l_type, float cutoff,
                           float kappa, float tau, string kernel_type,
                           const Point& global_min, const Point& global_max) {

    if (pro_indices.empty() || lig_indices.empty()) return vector<float>(10, 0.0f);

    auto pro_pos = pro_file.get_positions();
    auto lig_pos = lig_file.get_positions();

    vector<int> valid_pro_indices;
    for (int idx : pro_indices) {
        double x = pro_pos[idx][0];
        double y = pro_pos[idx][1];
        double z = pro_pos[idx][2];
        
        if (x > global_min.x - cutoff && x < global_max.x + cutoff &&
            y > global_min.y - cutoff && y < global_max.y + cutoff &&
            z > global_min.z - cutoff && z < global_max.z + cutoff) {
            valid_pro_indices.push_back(idx);
        }
    }

    float total_pair_count = (float)valid_pro_indices.size() * (float)lig_indices.size();

    if (valid_pro_indices.empty()) return vector<float>(10, 0.0f);

    int n_p = valid_pro_indices.size();
    int n_l = lig_indices.size();
    int N = n_p + n_l;

    float r_pro = get_radius(p_type);
    float r_lig = get_radius(l_type);
    float covalent_restriction = r_pro + r_lig + SIGMA;

    MatrixXd adj_mat(N, N);
    adj_mat.setZero();
    
    bool has_interaction = false;

    for (int i = 0; i < n_p; i++) {
        int p_real_idx = valid_pro_indices[i];
        Point p1 = { (float)pro_pos[p_real_idx][0], (float)pro_pos[p_real_idx][1], (float)pro_pos[p_real_idx][2] };

        for (int j = 0; j < n_l; j++) {
            int l_real_idx = lig_indices[j];
            Point p2 = { (float)lig_pos[l_real_idx][0], (float)lig_pos[l_real_idx][1], (float)lig_pos[l_real_idx][2] };
            
            float d_calc = distance_calc(p1, p2);

            if (d_calc <= cutoff && d_calc > covalent_restriction) {
                float eta = r_lig + r_pro;
                // PASS KERNEL TYPE HERE
                float weight = kernel(d_calc, eta, kappa, tau, kernel_type);
                
                int row = i;
                int col = n_p + j;
                
                adj_mat(row, col) = -weight; 
                adj_mat(col, row) = -weight; 
                has_interaction = true;
            }
        }
    }

    if (!has_interaction) {
        vector<float> zeros(10, 0.0f);
        zeros[0] = total_pair_count;
        return zeros;
    }

    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(adj_mat); 
    return calculate_stats(eigensolver.eigenvalues(), total_pair_count);
}

trajFile analyze_file(const std::string &file_name) {
    chemfiles::Trajectory trajectory(file_name);
    trajFile result;
    result.frame = trajectory.read(); 
    result.topology = result.frame.topology();
    return result;
}

// --- PYTHON API ---
vector<float> get_agl_scores(string pro_file, string lig_file, float cutoff, float kappa, float tau, string kernel_type) {
    chemfiles::set_warning_callback([](std::string){});

    trajFile pdb = analyze_file(pro_file);
    trajFile mol2 = analyze_file(lig_file);
    
    auto pro_groups = get_atom_groups(pdb, false);
    auto lig_groups = get_atom_groups(mol2, true);

    auto lig_pos = mol2.get_positions();
    Point global_min = {1e9, 1e9, 1e9};
    Point global_max = {-1e9, -1e9, -1e9};

    for (size_t i = 0; i < mol2.frame.size(); i++) {
        if (lig_pos[i][0] < global_min.x) global_min.x = lig_pos[i][0];
        if (lig_pos[i][1] < global_min.y) global_min.y = lig_pos[i][1];
        if (lig_pos[i][2] < global_min.z) global_min.z = lig_pos[i][2];
        if (lig_pos[i][0] > global_max.x) global_max.x = lig_pos[i][0];
        if (lig_pos[i][1] > global_max.y) global_max.y = lig_pos[i][1];
        if (lig_pos[i][2] > global_max.z) global_max.z = lig_pos[i][2];
    }

    vector<std::pair<string, string>> tasks;
    for (const string& p : PROTEIN_TYPES) {
        for (const string& l : LIGAND_ELEMS) {
            tasks.push_back({p, l});
        }
    }

    vector<float> global_features(tasks.size() * 10);

    py::gil_scoped_release release;

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < tasks.size(); i++) {
        string p_type = tasks[i].first;
        string l_elem = tasks[i].second;

        vector<int> p_idxs = pro_groups.count(p_type) ? pro_groups.at(p_type) : vector<int>{};
        vector<int> l_idxs = lig_groups.count(l_elem) ? lig_groups.at(l_elem) : vector<int>{};

        // PASS KERNEL TYPE
        vector<float> feats = analyze_pair(pdb, p_idxs, mol2, l_idxs, p_type, l_elem, cutoff, kappa, tau, kernel_type, global_min, global_max);

        size_t offset = i * 10;
        for (int k = 0; k < 10; k++) {
            global_features[offset + k] = feats[k];
        }
    }
    
    py::gil_scoped_acquire acquire;
    return global_features;
}

PYBIND11_MODULE(agl_cpp, m) {
    m.doc() = "AGL Score Calculator (C++ Accelerated)";

    m.def("get_agl_scores", &get_agl_scores, "Calculate Flat AGL Features",
          py::arg("protein_file"), 
          py::arg("ligand_file"),
          py::arg("cutoff") = 12.0,
          py::arg("kappa") = 6.0,
          py::arg("tau") = 4.0,
          py::arg("kernel_type") = "exponential"); // DEFAULT ARGUMENT
}