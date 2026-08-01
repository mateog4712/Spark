#ifndef PART_FUNC
#define PART_FUNC

#include "base_types.hh"
#include "matrix.hh"
#include "trace_arrow.hh"
#include "sparse_tree.hh"

#include "ViennaRNA/loops.hh"
#include "ViennaRNA/pair_mat.hh"
#include "ViennaRNA/params/io.hh"
#include <vector>

#define debug 0
#define INFover2 5000000 /* (INT_MAX/20) */

/*
 * If the global use_mfelike_energies flag is set, truncate doubles to int
 * values and cast back to double. This makes the energy parameters of the
 * partition (folding get_scaled_exp_params()) compatible with the mfe folding
 * parameters (get_scaled_exp_params()), e.g. for explicit partition function
 * computations.
 */
#define TRUNC_MAYBE(X) ((!pf_smooth) ? (double)((int)(X)) : (X))
/* Rescale Free energy contribution according to deviation of temperature from measurement conditions */
#define RESCALE_dG(dG, dH, dT) ((dH) - ((dH) - (dG)) * dT)

/*
 * Rescale Free energy contribution according to deviation of temperature from measurement conditions
 * and convert it to Boltzmann Factor for specific kT
 */
#define RESCALE_BF(dG, dH, dT, kT) (exp(-TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) * 10. / kT))

struct quatret_PF {
    cand_pos_t first;
    pf_t second;
    pf_t third;
    pf_t fourth;
    quatret_PF() {
        first = 1;
        second = 2;
        third = 3;
        fourth = 4;
    }
    quatret_PF(cand_pos_t x, pf_t y, pf_t z, pf_t w) {
        first = x;
        second = y;
        third = z;
        fourth = w;
    }
};

typedef std::pair<cand_pos_t, energy_t> cand_entry_t;
typedef std::vector<cand_entry_t> cand_list_t;

typedef quatret_PF cand_entry_td1;
typedef std::vector<cand_entry_td1> cand_list_td1;

class Spark_PF {
  public:

    Spark_PF(const std::string &seq, std::string restricted, sparse_tree *tree, bool pseudoknot, bool pk_only, bool garbage_collect,bool mark_candidates);
    ~Spark_PF();
    void print_candidates();
    void print_trace_arrows();
    energy_t fold();
    const std::string& trace_back();


  private:
    std::string seq_;
    cand_pos_t n_;

    short *S_;
    short *S1_;

    vrna_exp_param_t *exp_params_;
    sparse_tree* tree;

    std::string structure_;
    std::string restricted_;

    bool garbage_collect_;
    bool mark_candidates_;
    bool pseudoknot;
    bool pk_only;

    LocARNA::Matrix<energy_t> V_; // store V[i..i+MAXLOOP-1][1..n]
    std::vector<energy_t> W_;
    std::vector<energy_t> WM_;
    std::vector<energy_t> WM2_;

    std::vector<energy_t> dmli1_; // WM2 from 1 iteration ago
    std::vector<energy_t> dmli2_; // WM2 from 2 iterations ago

    // Pseudoknot portion
    LocARNA::Matrix<energy_t> VP_; // store VP[i..i+MAXLOOP-1][1..n]
    std::vector<energy_t> WVe_;
    std::vector<energy_t> WMB_;
    std::vector<energy_t> WMBP_;
    std::vector<energy_t> WMBA_;
    std::vector<energy_t> WI_;
    std::vector<energy_t> dwi1_; // WI from 1 iteration ago
    std::vector<energy_t> WIP_;
    std::vector<energy_t> dwip1_; // WIP from 1 iteration ago
    std::vector<energy_t> WV_;
    std::vector<energy_t> dwvp_; // WV from 1 iteration ago;

    std::vector<energy_t> WI_Bbp;  // WI from band borders on left
    std::vector<energy_t> WIP_Bbp; // WIP from band borders on left
    std::vector<energy_t> WIP_Bp;  // WIP from band borders on the right

    TraceArrows ta_;
    TraceArrows taVP_;

    std::vector<cand_list_td1> CL_;
    std::vector<cand_list_t> CLVP_;
    std::vector<cand_list_t> CLWMB_;
    std::vector<cand_list_t> CLBE_;
    std::vector<cand_list_t> CLBEO_;

    std::vector<pf_t> scale;
    std::vector<pf_t> expMLbase;
    std::vector<pf_t> expcp_pen;
    std::vector<pf_t> expPUP_pen;

    // compare candidate list entries by keys (left index i) in descending order
    struct Cand_comp {
        bool operator()(const cand_entry_t &x, cand_pos_t y) const { return x.first > y; }
        bool operator()(const cand_entry_td1 &x, cand_pos_t y) const { return x.first > y; }
    } cand_comp;

    // Folding
    energy_t compute_internal(cand_pos_t i, cand_pos_t j, cand_pos_t &best_k, cand_pos_t &best_l, energy_t &best_e);
    energy_t compute_VP(cand_pos_t i, cand_pos_t j, cand_pos_t b_ij, cand_pos_t bp_ij, cand_pos_t Bp_ij, cand_pos_t B_ij);
    energy_t compute_VP_internal(cand_pos_t i, cand_pos_t j, cand_pos_t b_ij, cand_pos_t bp_ij, cand_pos_t Bp_ij, cand_pos_t B_ij, cand_pos_t &best_k, cand_pos_t &best_l, energy_t &best_e);
    energy_t compute_WV(cand_pos_t j, cand_pos_t bound_left, cand_pos_t bound_right);
    energy_t compute_WVe(cand_pos_t i, cand_pos_t j);
    void compute_WMB_case1(cand_pos_t j, energy_t &m1, energy_t &BE_en, cand_pos_t &best_border);
    energy_t compute_WMBP(cand_pos_t i, cand_pos_t j);
    energy_t compute_WMBA(cand_pos_t j);
    energy_t compute_WMB(cand_pos_t i, cand_pos_t j);
    energy_t compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp);

    // Traceback
    const std::vector<energy_t> recompute_WM(cand_pos_t i, cand_pos_t max_j);
    const std::vector<energy_t> recompute_WM2(cand_pos_t i, cand_pos_t max_j);
    void recompute_WMBP(cand_pos_t i, cand_pos_t max_j);
    void recompute_WI(cand_pos_t i, cand_pos_t max_j);
    void recompute_WIP(cand_pos_t i, cand_pos_t max_j);
    void recompute_WVe(cand_pos_t i, cand_pos_t max_j);
    void recompute_WV(cand_pos_t i, cand_pos_t max_j);
    void trace_V(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_W(cand_pos_t i, cand_pos_t j);
    void trace_WM(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_WM2(cand_pos_t i, cand_pos_t j);
    void trace_WMB(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_VP(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_WI(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_WIP(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_WV(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_WVe(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_WMBP(cand_pos_t i, cand_pos_t j, energy_t e);
    void trace_BE(cand_pos_t i, cand_pos_t ip, energy_t e);
    void trace_WMBA(cand_pos_t i, cand_pos_t j, energy_t e);

    // Util
    pf_t to_Energy(pf_t energy, cand_pos_t length);
    pf_t to_PF(pf_t energy, cand_pos_t length);
    void rescale_pk_globals();
    void exp_params_rescale(double mfe);
    void rotate_arrays();
    int compute_exterior_cases(cand_pos_t l, cand_pos_t j);
    void find_mb_dangle(const energy_t WM2ij, const energy_t WM2i1j, const energy_t WM2ij1, const energy_t WM2i1j1, const cand_pos_t i, const cand_pos_t j, cand_pos_t &k, cand_pos_t &l);
    pf_t HairpinE(cand_pos_t i, cand_pos_t j);
    pf_t ILoopE(const pair_type &ptype_closing, const cand_pos_t &i, const cand_pos_t &j, const cand_pos_t &k, const cand_pos_t &l);
    
    pf_t exp_Extloop(cand_pos_t i, cand_pos_t j);
    pf_t exp_MLstem(cand_pos_t i, cand_pos_t j);
    pf_t exp_Mbloop(cand_pos_t i, cand_pos_t j);

    /**
     * @brief Test existence of candidate. Used primarily for determining whether (i,j) is candidate for W/WM splits
     *
     */
    inline bool is_candidate(const std::vector<cand_list_td1> &CL, const Spark_PF::Cand_comp &cand_comp, cand_pos_t i, cand_pos_t j) {
        const cand_list_td1 &list = CL[j];

        auto it = std::lower_bound(list.begin(), list.end(), i, cand_comp);

        return it != list.end() && it->first == i;
    }
    inline bool is_candidate(const std::vector<cand_list_t> &CL, const Spark_PF::Cand_comp &cand_comp, cand_pos_t i, cand_pos_t j) {
        const cand_list_t &list = CL[j];

        auto it = std::lower_bound(list.begin(), list.end(), i, cand_comp);

        return it != list.end() && it->first == i;
    }

};
#endif