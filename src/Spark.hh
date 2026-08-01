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


struct quatret {
    cand_pos_t first;
    energy_t second;
    energy_t third;
    energy_t fourth;
    quatret() {
        first = 1;
        second = 2;
        third = 3;
        fourth = 4;
    }
    quatret(cand_pos_t x, energy_t y, energy_t z, energy_t w) {
        first = x;
        second = y;
        third = z;
        fourth = w;
    }
};

typedef std::pair<cand_pos_t, energy_t> cand_entry_t;
typedef std::vector<cand_entry_t> cand_list_t;

typedef quatret cand_entry_td1;
typedef std::vector<cand_entry_td1> cand_list_td1;

class Spark {
  public:

    Spark(const std::string &seq, std::string restricted, sparse_tree *tree, int dangles, bool pseudoknot, bool pk_only, bool garbage_collect,bool mark_candidates);
    ~Spark();
    void print_candidates();
    void print_trace_arrows();
    energy_t fold();
    const std::string& trace_back();


  private:
    std::string seq_;
    cand_pos_t n_;

    short *S_;
    short *S1_;

    vrna_param_t *params_;
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
    void rotate_arrays();
    int compute_exterior_cases(cand_pos_t l, cand_pos_t j);
    void find_mb_dangle(const energy_t WM2ij, const energy_t WM2i1j, const energy_t WM2ij1, const energy_t WM2i1j1, const cand_pos_t i, const cand_pos_t j, cand_pos_t &k, cand_pos_t &l);
    energy_t HairpinE(cand_pos_t i, cand_pos_t j);
    energy_t ILoopE(const pair_type &ptype_closing, const cand_pos_t &i, const cand_pos_t &j, const cand_pos_t &k, const cand_pos_t &l);
    energy_t E_MbLoop(const std::vector<energy_t> &dmli1, const std::vector<energy_t> &dmli2, cand_pos_t i, cand_pos_t j);
    energy_t E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1, const cand_pos_t i,const cand_pos_t j, Dangle &d);
    energy_t E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,cand_pos_t i, cand_pos_t j, Dangle &d);

    /**
     * @brief Test existence of candidate. Used primarily for determining whether (i,j) is candidate for W/WM splits
     *
     */
    inline bool is_candidate(const std::vector<cand_list_td1> &CL, const Spark::Cand_comp &cand_comp, cand_pos_t i, cand_pos_t j) {
        const cand_list_td1 &list = CL[j];

        auto it = std::lower_bound(list.begin(), list.end(), i, cand_comp);

        return it != list.end() && it->first == i;
    }
    inline bool is_candidate(const std::vector<cand_list_t> &CL, const Spark::Cand_comp &cand_comp, cand_pos_t i, cand_pos_t j) {
        const cand_list_t &list = CL[j];

        auto it = std::lower_bound(list.begin(), list.end(), i, cand_comp);

        return it != list.end() && it->first == i;
    }

    /**
     * @brief Register a candidate
     */
    inline void register_candidate(std::vector<cand_list_td1> &CL, cand_pos_t const &i, cand_pos_t const &j, energy_t const &e, energy_t const &wmij, energy_t const &wij) {
        assert(i <= j + TURN + 1);
        CL[j].emplace_back(cand_entry_td1(i, e, wmij, wij));
    }
    inline void register_candidate(std::vector<cand_list_t> &CL, cand_pos_t const &i, cand_pos_t const &j, energy_t const &e) {
        assert(i <= j + TURN + 1);
        CL[j].emplace_back(cand_entry_t(i, e));
    }
 
    /**
     * @brief Sums the number of Candidates at each index over all indices
     *
     * @param CL_ Candidate list
     * @return total number of candidates
     */
    inline cand_pos_t num_of_candidates(const std::vector<cand_list_td1> &CL_) {
        cand_pos_t c = 0;
        for (const cand_list_td1 &x : CL_) {
            c += x.size();
        }
        return c;
    }
    inline cand_pos_t num_of_candidates(const std::vector<cand_list_t> &CL_) {
        cand_pos_t c = 0;
        for (const cand_list_t &x : CL_) {
            c += x.size();
        }
        return c;
    }
    /**
     * @brief Finds the size of allocated storage capacity across all indices
     *
     * @param CL_ Candidate List
     * @return the amount of allocated storage
     */
    inline cand_pos_t capacity_of_candidates(const std::vector<cand_list_td1> &CL_) {
        cand_pos_t c = 0;
        for (const cand_list_td1 &x : CL_) {
            c += x.capacity();
        }
        return c;
    }
    inline cand_pos_t capacity_of_candidates(const std::vector<cand_list_t> &CL_) {
        cand_pos_t c = 0;
        for (const cand_list_t &x : CL_) {
            c += x.capacity();
        }
        return c;
    }

};