/**
 * @mainpage
 *
 * Space-efficient sparse variant of an RNA (loop-based) free energy
 * minimization algorithm (RNA folding equivalent to the Zuker
 * algorithm).
 *
 * The results are equivalent to HFold.
 */
#include "PK_globals.hh"
#include "Spark.hh"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

Spark::Spark(const std::string &seq, std::string restricted, sparse_tree *tree, int dangles, bool pseudoknot, bool pk_only, bool garbage_collect,bool mark_candidates) : seq_(seq), n_(seq.length()), params_(vrna_params(NULL)), garbage_collect_(garbage_collect), mark_candidates_(mark_candidates), ta_(n_), taVP_(n_) {
    make_pair_matrix();

    S_ = encode_sequence(seq.c_str(), 0);
    S1_ = encode_sequence(seq.c_str(), 1);
    this->tree = tree;
    params_->model_details.dangles = dangles;
    this->pseudoknot= pseudoknot;
    this->pk_only = pk_only;

    V_.resize(MAXLOOP + 1, n_ + 1, INF);
    W_.resize(n_ + 1, 0);
    WM_.resize(n_ + 1, INF);
    WM2_.resize(n_ + 1, INF);
    dmli1_.resize(n_ + 1, INF);
    dmli2_.resize(n_ + 1, INF);

    // Pseudoknot portion

    VP_.resize(MAXLOOP + 1, n_ + 1, INF);
    WVe_.resize(n_ + 1, INF);
    WMB_.resize(n_ + 1, INF);
    WMBP_.resize(n_ + 1, INF);
    WMBA_.resize(n_ + 1, INF);
    WI_.resize(n_ + 1, 0);
    dwi1_.resize(n_ + 1, 0);
    WIP_.resize(n_ + 1, INF);
    dwip1_.resize(n_ + 1, INF);
    WV_.resize(n_ + 1, INF);
    dwvp_.resize(n_ + 1, INF);

    WI_Bbp.resize(n_ + 1, 0);
    WIP_Bbp.resize(n_ + 1, INF);
    WIP_Bp.resize(n_ + 1, INF);

    // init candidate lists
    CL_.resize(n_ + 1);
    CLWMB_.resize(n_ + 1);
    CLVP_.resize(n_ + 1);
    CLBE_.resize(n_ + 1);
    CLBEO_.resize(n_ + 1);

    resize(ta_, n_ + 1);
    resize(taVP_, n_ + 1);

    // resize(ta_dangle_,n_+1);

    restricted_ = restricted;

    int count = 0;
    for (cand_pos_t i = 1; i <= n_; ++i) {
        if (tree->tree[i].pair > i || (tree->tree[i].pair < i && tree->tree[i].pair > 0)) count = 4;
        if (tree->tree[i].pair < 0 && count > 0) {
            WI_Bbp[i] = (5 - count) * PUP_penalty;
            count--;
        }
    }
}

Spark::~Spark() {
    free(params_);
    free(S_);
    free(S1_);
}

void Spark::print_candidates(){
    std::cout << "Can num:\t" << num_of_candidates(CL_);
    std::cout << ", Can cap:\t" << capacity_of_candidates(CL_) << std::endl;

    std::cout << "\nPsuedoknotted\n" << std::endl;
    std::cout << std::endl;
    std::cout << "WMB Can num:\t" << num_of_candidates(CLWMB_);
    std::cout << ", WMB Can cap:\t" << capacity_of_candidates(CLWMB_) << std::endl;
    std::cout << "VP Can num:\t" << num_of_candidates(CLVP_);
    std::cout << ", VP Can cap:\t" << capacity_of_candidates(CLVP_) << std::endl;
    std::cout << "BE Can num:\t" << num_of_candidates(CLBE_);
    std::cout << ", BE Can cap:\t" << capacity_of_candidates(CLBE_) << std::endl;
}

void Spark::print_trace_arrows(){
    std::cout << "PKfree: ";
    std::cout << "TA cnt:\t" << sizeT(ta_);
    std::cout << ", TA max:\t" << maxT(ta_);
    std::cout << ", TA av:\t" << avoidedT(ta_);
    std::cout << ", TA rm:\t" << erasedT(ta_) << std::endl;
    std::cout << "PK: " << std::endl;
    std::cout << "TAs num:\t" << sizeT(taVP_);
    std::cout << ", TAs cap:\t" << capacityT(taVP_);
    std::cout << ", TA av:\t" << avoidedT(taVP_);
    std::cout << ", TA rm:\t" << erasedT(taVP_) << std::endl;

    // std::cout << "TAs num:\t" << sizeT(spark.ta_) << std::endl;
    //     std::cout << "TAs cap:\t" << capacityT(spark.ta_) << std::endl;
}

/**
 * @brief Rotate WM2 and WI arrays to store the previous and previous previous iterations
 * @param WM2 WM2 array
 * @param dmli1 WM2 from one iteration ago
 * @param dmli2 WM2 from two iterations ago
 * @param WI WI array
 * @param dwi1 WM2 from one iteration ago
 * @param WIP WIP array
 * @param dwip1 WIP from one iteration ago
 * @param WV WV array
 * @param dwvp WV from one iteration ago
 */
void Spark::rotate_arrays() {
    dmli2_.swap(dmli1_);
    dmli1_.swap(WM2_);
    dwi1_.swap(WI_);
    dwip1_.swap(WIP_);
    dwvp_.swap(WV_);
}

/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
 */
energy_t Spark::HairpinE(cand_pos_t i, cand_pos_t j) {

    const pair_type ptype_closing = pair[S_[i]][S_[j]];

    if (ptype_closing == 0) return INF;

    return E_Hairpin(j-i-1, ptype_closing, S1_[i + 1], S1_[j - 1], &seq_.c_str()[i - 1], const_cast<vrna_param_t *>(params_));
}

/**
 * @brief Returns the internal loop energy for a given i.j and k.l
 *
 */
energy_t Spark::ILoopE(const pair_type &ptype_closing, const cand_pos_t &i, const cand_pos_t &j, const cand_pos_t &k, const cand_pos_t &l) {
    assert(ptype_closing > 0);
    assert(1 <= i);
    assert(i < k);
    assert(k < l);
    assert(l < j);

    // note: enclosed bp type 'turned around' for lib call
    const pair_type ptype_enclosed = rtype[pair[S_[k]][S_[l]]];

    if (ptype_enclosed == 0) return INF;

    return E_IntLoop(k - i - 1, j - l - 1, ptype_closing, ptype_enclosed, S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1], const_cast<vrna_param_t *>(params_));
}


/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
*/
energy_t Spark::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1, const cand_pos_t i,const cand_pos_t j, Dangle &d){

	energy_t e = INF;

    auto consider = [&](energy_t v, bool valid, pair_type tt, base_type s5, base_type s3, Dangle &d, int d_type) {
        if (!valid || v == INF) return;
        if(v + E_ExtLoop(tt, s5, s3, params_) < e){
            e = v + E_ExtLoop(tt, s5, s3, params_);
            d = d_type;
        }
        e = std::min(e, v + E_ExtLoop(tt, s5, s3, params_));
    };
	base_type si1  = i > 1 ? S_[i-1] : -1;
    base_type sj1  = j < n_ ? S_[j+1] : -1;
    base_type si = S_[i];
    base_type sj = S_[j];

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(vij, ((tree->tree[i].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i].pair == j && tree->tree[j].pair == i)), pair[S_[i]][S_[j]], dangle2 ? si1 : -1, dangle2 ? sj1 : -1,d, dangle2 ? d : 0);
	if (dangle1) {
        consider(vi1j,j-i-1>TURN && (((tree->tree[i + 1].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i + 1].pair == j)) && tree->tree[i].pair < 0), pair[S_[i+1]][S_[j]], si, -1,d,1);
        consider(vij1,j-1-i>TURN && (((tree->tree[i].pair < -1 && tree->tree[j - 1].pair < -1) || (tree->tree[i].pair == j - 1)) && tree->tree[j].pair < 0), pair[S_[i]][S_[j-1]], -1, sj,d,2);
        consider(vi1j1,j-1-i-1>TURN && (((tree->tree[i + 1].pair < -1 && tree->tree[j - 1].pair < -1) || (tree->tree[i + 1].pair == j - 1)) &&tree-> tree[i].pair < 0 && tree->tree[j].pair < 0), pair[S_[i+1]][S_[j-1]], si, sj,d,3);
    }
	return e;
}
/**
* @brief Computes the multiloop V contribution. This gives back essentially VM(i,j).
* 
*/
energy_t Spark::E_MbLoop(const std::vector<energy_t> &dmli1, const std::vector<energy_t> &dmli2, cand_pos_t i, cand_pos_t j){
	energy_t e = INF;

    bool pairable = (tree->tree[i].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i].pair == j);
    pair_type tt = pair[S_[j]][S_[i]];
    base_type si1 = S_[i+1];
    base_type sj1 = S_[j-1];

	auto consider = [&](energy_t v, bool check, base_type s5, base_type s3, int ml_count) {
        if (check && v == INF) return;
        e = std::min(e, v + E_MLstem(tt, s5, s3, params_) + params_->MLclosing + ml_count * params_->MLbase);
    };

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(dmli1[j-1],pairable, dangle2 ? sj1 : -1, dangle2 ? si1 : -1, 0);
	if(dangle1){
		// ML pair 5 — closing (i,j) with mb part [i+2, j-1]
		consider(dmli2[j-1],pairable && tree->tree[i+1].pair < 0, -1, si1, 1);
        // ML pair 3 — closing (i,j) with mb part [i+1, j-2]
        consider(dmli1[j-2],pairable && tree->tree[j-1].pair < 0, sj1, -1, 1);
        // ML pair 53 — closing (i,j) with mb part [i+2, j-2]
        consider(dmli2[j-2],pairable && tree->tree[i+1].pair < 0 && tree->tree[j-1].pair < 0, sj1, si1, 2);
	}
	return e;
}

/**
 * @brief Gives the WM(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
*/
energy_t Spark::E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,cand_pos_t i, cand_pos_t j, Dangle &d){

	energy_t e = INF;

    auto consider = [&](energy_t v, bool valid, pair_type type, base_type s5, base_type s3, int ml_count, Dangle &d, int d_type) {
        if (!valid || v == INF) return;
        if(v + E_MLstem(type, s5, s3, params_) + ml_count * params_->MLbase < e){
            e = v + E_MLstem(type, s5, s3, params_) + ml_count * params_->MLbase;
            d = d_type;
        }
    };

	base_type si1  = i > 1 ? S_[i-1] : -1;
    base_type sj1  = j < n_ ? S_[j+1] : -1;
    base_type si = S_[i];
    base_type sj = S_[j];

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(vij, (tree->tree[i].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i].pair == j), pair[S_[i]][S_[j]], dangle2 ? si1 : -1, dangle2 ? sj1 : -1, 0,d,dangle2 ? d : 0);
	if (dangle1) {
		consider(vi1j,j-i-1>TURN && (((tree->tree[i + 1].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i + 1].pair == j)) && tree->tree[i].pair < 0), pair[S_[i+1]][S_[j]], si, -1, 1,d,1);
        consider(vij1,j-1-i>TURN && (((tree->tree[i].pair < -1 && tree->tree[j - 1].pair < -1) || (tree->tree[i].pair == j - 1)) && tree->tree[j].pair < 0), pair[S_[i]][S_[j-1]], -1, sj, 1,d,2);
        consider(vi1j1,j-1-i-1>TURN && (((tree->tree[i + 1].pair < -1 && tree->tree[j - 1].pair < -1) || (tree->tree[i + 1].pair == j - 1)) && tree->tree[i].pair < 0 && tree->tree[j].pair < 0), pair[S_[i+1]][S_[j-1]], si, sj, 2,d,3);
	}
    return e;
}

/**
 * In cases where the band border is not found, if specific cases are met, the value is Inf(n) not -1.
 * Mateo Jan 2025: Added to Fix WMBP problem
 */
int Spark::compute_exterior_cases(cand_pos_t l, cand_pos_t j) {

    // Case 1 -> l is not covered
    bool case1 = tree->tree[l].parent->index <= 0;
    // Case 2 -> l is paired
    bool case2 = tree->tree[l].pair > 0;
    // Case 3 -> l is part of a closed subregion
    bool case3 = 0;
    // Case 4 -> l.bp(l) i.e. l.j does not cross anything -- could I compare parents instead?
    bool case4 = j < tree->Bp(l, j);
    // By bitshifting each one, we have a more granular idea of what cases fail and is faster than branching
    return (case1 << 3) | (case2 << 2) | (case3 << 1) | case4;
}

/**
 * @brief Recompute row of WM. This is used in the traceback when we haved decided the current i.j pair closes a multiloop,
 * and the WM energies need to be recomputed fom the candidates.
 *
 * @param WM WM array
 * @param CL Candidate List
 * @param i Current i
 * @param max_j Current j
 */
const std::vector<energy_t> Spark::recompute_WM(cand_pos_t i, cand_pos_t max_j) {
    assert(i >= 1);
    assert(max_j <= spark.n_);
    std::vector<energy_t> temp = WM_;

    for (cand_pos_t j = i - 1; j <= std::min(i + TURN, max_j); j++) {
        temp[j] = INF;
    }

    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wm = INF;
        for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {
            cand_pos_t k = it->first;
            energy_t v_kj = it->third >> 2;

            bool can_pair = tree->up[k - 1] >= (k - i);
            if (can_pair) wm = std::min(wm, static_cast<energy_t>(params_->MLbase * (k - i)) + v_kj);
            wm = std::min(wm, temp[k - 1] + v_kj);
        }
        for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first > i + TURN + 1; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second + PSM_penalty + b_penalty;
            bool can_pair = tree->up[k - 1] >= (k - i);
            if (can_pair) wm = std::min(wm, static_cast<energy_t>(params_->MLbase * (k - i)) + wmb_kj);

            wm = std::min(wm, temp[k - 1] + wmb_kj);
        }
        if (tree->tree[j].pair < 0) wm = std::min(wm, temp[j - 1] + params_->MLbase);
        temp[j] = wm;
    }
    return temp;
}

/**
 * @brief Recompute row of WM2. This is used in the traceback when we haved decided the current i.j pair closes a multiloop,
 * and the WM2 energies need to be recomputed fom the candidates to get the corresponding energy for it.
 *
 * @param WM WM array
 * @param WM2 WM2 array
 * @param CL Candidate List
 * @param i Current i
 * @param max_j Current j
 */
const std::vector<energy_t> Spark::recompute_WM2(cand_pos_t i, cand_pos_t max_j) {
    assert(i >= 1);
    assert(max_j <= spark.n_);
    std::vector<energy_t> temp = WM2_;

    for (cand_pos_t j = i - 1; j <= std::min(i + 2 * TURN + 2, max_j); j++) {
        temp[j] = INF;
    }

    for (cand_pos_t j = i + 2 * TURN + 3; j <= max_j; j++) {
        energy_t wm2 = INF;
        for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first > i + TURN + 1; ++it) {

            cand_pos_t k = it->first;
            energy_t v_kj = it->third >> 2;

            wm2 = std::min(wm2, WM_[k - 1] + v_kj);
        }
        for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second + PSM_penalty + b_penalty;
            bool can_pair = tree->up[k - 1] >= (k - i);

            wm2 = std::min(wm2, WM_[k - 1] + wmb_kj);
            if (can_pair) wm2 = std::min(wm2, static_cast<energy_t>(params_->MLbase * (k - i)) + wmb_kj);
        }
        if (tree->tree[j].pair < 0) wm2 = std::min(wm2, temp[j - 1] + params_->MLbase);
        temp[j] = wm2;
    }
    return temp;
}

/**
 * @brief Recompute row of WMBP. This is used in the traceback when we haved decided the current ij pair closes a psuedoknot,
 * and the WMBP energies need to be recomputed such that the pseudoknot can be broken down. This happens infrequently enough such that it not too
 * costly
 *
 * @param CL V candidate list
 * @param CLVP VP candidate list
 * @param CLWMB WMB Candidate list
 * @param CLBE BE candidate list
 * @param WI_Bbp array which holds the values where the left index is either an outer closing base or an inner opening base
 * @param i Current i
 * @param max_j Current j
 */
void Spark::recompute_WMBP(cand_pos_t i, cand_pos_t max_j) {
    assert(i >= 1);
    assert(max_j <= spark.n_);

    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wmbp = INF;
        energy_t wmba = INF;

        energy_t vp_ij = INF;
        for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t v_kj = it->second;

            wmba = std::min(wmba, WMBP_[k - 1] + v_kj + PPS_penalty);
        }
        for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second;

            wmba = std::min(wmba, WMBP_[k - 1] + wmb_kj + PPS_penalty + PSM_penalty);
        }
        // m2
        if (tree->tree[j].pair < 0) {
            for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it && it->first >= i; ++it) {
                cand_pos_t k = it->first;
                if (k == i) vp_ij = it->second; // second?
                cand_pos_t bp_ik = tree->bp(i, k);
                cand_pos_t Bp_kj = tree->Bp(k, j);
                cand_pos_t b_ij = tree->b(i, j);
                int ext_case = compute_exterior_cases(k, j);
                if ((b_ij > 0 && k < b_ij) || (b_ij < 0 && ext_case == 0)) {
                    if (bp_ik >= 0 && k > bp_ik && Bp_kj > 0 && k < Bp_kj) {
                        energy_t BE_energy = INF;
                        cand_pos_t B_kj = tree->B(k, j);
                        cand_pos_t b_kj = (B_kj > 0) ? tree->tree[B_kj].pair : -2;
                        for (auto it2 = CLBE_[Bp_kj].begin(); CLBE_[Bp_kj].end() != it2; ++it2) {
                            cand_pos_t l = it2->first;
                            if (l == b_kj) {
                                BE_energy = it2->second;
                                break;
                            }
                        }
                        wmbp = std::min(wmbp, BE_energy + WMBA_[k - 1] + it->second + 2 * PB_penalty);
                    }
                }
            }
        }
        // m3
        if (tree->tree[j].pair < 0 && tree->tree[i].pair >= 0 && tree->tree[i].pair < j) {
            for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it && it->first >= i; ++it) {
                cand_pos_t k = it->first;

                cand_pos_t bp_ik = tree->bp(i, k);
                if (bp_ik >= 0 && k + TURN <= j) {
                    cand_pos_t Bp_ik = tree->tree[bp_ik].pair;
                    energy_t BE_energy = INF;
                    for (auto it2 = CLBE_[Bp_ik].begin(); CLBE_[Bp_ik].end() != it2; ++it2) {
                        cand_pos_t l = it2->first;
                        if (l == i) {
                            BE_energy = it2->second;
                            break;
                        }
                    }
                    wmbp = std::min(wmbp, BE_energy + WI_Bbp[k - 1] + it->second + 2 * PB_penalty);
                }
            }
        }
        wmbp = std::min(wmbp, vp_ij + PB_penalty);
        if (tree->tree[j].pair < 0) wmba = std::min(wmba, WMBA_[j - 1] + PUP_penalty);
        wmba = std::min(wmba, wmbp);
        WMBA_[j] = wmba;
        WMBP_[j] = wmbp;
    }
}

/**
 * @brief Recompute row of WI. This is used in the traceback when we haved decided the current ij pair is a crossing base pair, we need to calculate
 * the energy, and the WI energies need to be recomputed such that the VP energy can be broken down. This happens only when the calculations are
 * confirmed to be needed such that it not too costly
 *
 * @param CL V candidate list
 * @param CLWMB WMB candidate list
 * @param i Current i
 * @param max_j Current j
 */
void Spark::recompute_WI(cand_pos_t i, cand_pos_t max_j) {
    assert(i >= 1);
    assert(max_j <= n_);

    // Causes vector resize error if the ifs are not there because if i is close to n, it would go past the bounds
    for (cand_pos_t j = 0; j < 4; ++j) {
        if (i + j < n_) {
            WI_[i + j] = (j + 1) * PUP_penalty;
        }
    }

    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wi = INF;
        for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t v_kj = it->second + PPS_penalty;

            wi = std::min(wi, WI_[k - 1] + v_kj);
            bool can_pair = tree->up[k - 1] >= (k - i);
            if (can_pair) wi = std::min(wi, static_cast<energy_t>(PUP_penalty * (k - i)) + v_kj);
        }
        for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second + PSM_penalty + PPS_penalty;

            wi = std::min(wi, WI_[k - 1] + wmb_kj);
            bool can_pair = tree->up[k - 1] >= (k - i);
            if (can_pair) wi = std::min(wi, static_cast<energy_t>(PUP_penalty * (k - i)) + wmb_kj);
        }
        if (tree->tree[j].pair < 0) wi = std::min(wi, WI_[j - 1] + PUP_penalty);
        WI_[j] = wi;
    }
}

/**
 * @brief Recompute row of WIP. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band or part of a
 * BE, we need to calculate the energy, and the WIP energies need to be recomputed such that the VP or BE energy can be broken down. This happens only
 * when the calculations are confirmed to be needed such that it not too costly
 *
 * @param CL V candidate list
 * @param CLWMB WMB candidate list
 * @param i Current i
 * @param max_j Current j
 */
void Spark::recompute_WIP(cand_pos_t i, cand_pos_t max_j){
    assert(i >= 1);
    assert(max_j <= spark.n_);
    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wip = INF;
        for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t v_kj = it->second + bp_penalty;

            wip = std::min(wip, WIP_[k - 1] + v_kj);
            bool can_pair = tree->up[k - 1] >= (k - i);
            if (can_pair) wip = std::min(wip, static_cast<energy_t>(cp_penalty * (k - i)) + v_kj);
        }
        for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second + PSM_penalty + bp_penalty;

            wip = std::min(wip, WIP_[k - 1] + wmb_kj);
            bool can_pair = tree->up[k - 1] >= (k - i);
            if (can_pair) wip = std::min(wip, static_cast<energy_t>(cp_penalty * (k - i)) + wmb_kj);
        }
        if (tree->tree[j].pair < 0) wip = std::min(wip, WIP_[j - 1] + cp_penalty);
        WIP_[j] = wip;
    }
}
/**
 * @brief Recompute row of WVe. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band, we need to
 * calculate the energy, and the WVe energies need to be recomputed such that the VP energy can be broken down. This happens infrequently such that it
 * not too costly
 *
 * @param CLVP VP candidate list
 * @param i Current i
 * @param max_j Current j
 */
void Spark::recompute_WVe(cand_pos_t i, cand_pos_t max_j) {
    assert(i >= 1);
    assert(max_j <= spark.n_);
    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wve = INF;
        for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it && it->first >= i; ++it) {
            cand_pos_t k = it->first;
            energy_t vp_kj = it->second;

            bool can_pair = tree->up[k - 1] >= (k - i);
            if (can_pair) wve = std::min(wve, static_cast<energy_t>(cp_penalty * (k - i)) + vp_kj);
        }
        if (tree->tree[j].pair < 0) wve = std::min(wve, WVe_[j - 1] + cp_penalty);
        WVe_[j] = wve;
    }
}
/**
 * @brief Recompute row of WV. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band, we need to
 * calculate the energy, and the WV energies need to be recomputed such that the VP energy can be broken down. This happens infrequently such that it
 * not too costly
 *
 * @param CL V candidate list
 * @param CLWMB WMB candidate list
 * @param i Current i
 * @param CLVP VP candidate list
 * @param WIP WIP energies
 * @param i Current i
 * @param max_j Current j
 */
void Spark::recompute_WV(cand_pos_t i, cand_pos_t max_j) {
    assert(i >= 1);
    assert(max_j <= n_);

    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        cand_pos_t bound_right = std::max(tree->bp(i, j), tree->B(i, j));
        cand_pos_t bound_left = std::min((cand_pos_tu)tree->Bp(i, j), (cand_pos_tu)tree->b(i, j));

        energy_t wv = INF;
        if (!tree->weakly_closed(i, j)) {
            for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first > bound_right; ++it) {

                cand_pos_t k = it->first;
                // energy_t v_kj = it->third >> 2;
                energy_t v_kj = it->second;
                // Dangle d = it->third & 3;
                // cand_pos_t num = 0;
                // if(d == 1 || d==2) num = 1;
                // else if(d==3 && params->model_details.dangles == 1) num =2;
                // energy_t fix = num*cp_penalty - num*params->MLbase;

                wv = std::min(wv, WVe_[k - 1] + v_kj + bp_penalty);
                wv = std::min(wv, WV_[k - 1] + v_kj + bp_penalty);
            }

            for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it; ++it) {

                cand_pos_t k = it->first;
                energy_t vp_kj = it->second;
                if (k < bound_left) wv = std::min(wv, WIP_[k - 1] + vp_kj);
            }
            for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first > bound_right; ++it) {

                cand_pos_t k = it->first;
                energy_t wmb_kj = it->second;

                wv = std::min(wv, WVe_[k - 1] + wmb_kj + PSM_penalty + bp_penalty);
                wv = std::min(wv, WV_[k - 1] + wmb_kj + PSM_penalty + bp_penalty);
            }
            if (tree->tree[j].pair < 0) wv = std::min(wv, WV_[j - 1] + cp_penalty);
            WV_[j] = wv;
        }
    }
}

/**
 * @brief Sparse function to move through WMB case 1 without linear check
 *
 * @param j right closing base pair
 */
void Spark::compute_WMB_case1(cand_pos_t j, energy_t &m1, energy_t &BE_en, cand_pos_t &best_border) {
    // We are moving through the elements of BE/ the stucture in G
    // We are taking advantage of the fact that WMBA holds both the actual pseudoknotted base pair from WMB' and anything to the right of it
    // If there exists candidates in BEO, then they encompass the areas to be sparsely decomposed
    // We can calculate the BE from the candidate energy and the WMBA encompasses all values within and can be sparsely decomposed
    cand_pos_t bp_j = tree->tree[j].pair; // j' or opening base pair
    for (auto it = CLBEO_[bp_j].begin(); CLBEO_[bp_j].end() != it; ++it) {
        cand_pos_t candidate_index = it->first; // j or an inner closing base pair
        energy_t BE_energy = it->second;
        energy_t WMBA_energy = WMBA_[candidate_index - 1];
        if (BE_energy + WMBA_energy < m1) {
            m1 = BE_energy + WMBA_energy;
            BE_en = BE_energy;
            best_border = candidate_index;
        }
    }
}
/**
 * @brief Determines the type of dangle being used for a closing multiloop while in traceback.
 *
 * @param WM2ij The WM2 energy for the region [i,j]
 * @param WM2i1j The WM2 energy for the region [i+1,j]
 * @param WM2ij1 The WM2 energy for the region [i,j-1]
 * @param WM2i1j1 The WM2 energy for the region [i+1,j-1]
 */
void Spark::find_mb_dangle(const energy_t WM2ij, const energy_t WM2i1j, const energy_t WM2ij1, const energy_t WM2i1j1, const cand_pos_t i, const cand_pos_t j, cand_pos_t &k, cand_pos_t &l) {

    const pair_type tt = pair[S_[j]][S_[i]];
    const energy_t e1 = WM2ij + E_MLstem(tt, -1, -1, params_);
    const energy_t e2 = WM2i1j + E_MLstem(tt, -1, S_[i + 1], params_) + params_->MLbase;
    const energy_t e3 = WM2ij1 + E_MLstem(tt, S_[j - 1], -1, params_) + params_->MLbase;
    const energy_t e4 = WM2i1j1 + E_MLstem(tt, S_[j - 1], S_[i + 1], params_) + 2*params_->MLbase;
    energy_t e = e1;

    if (e2 < e && tree->tree[i+1].pair < 0) {
        e = e2;
        k = i + 2;
        l = j - 1;
    }
    if (e3 < e && tree->tree[j-1].pair < 0) {
        e = e3;
        k = i + 1;
        l = j - 2;
    }
    if (e4 < e && tree->tree[i+1].pair < 0 && tree->tree[j-1].pair < 0) {
        e = e4;
        k = i + 2;
        l = j - 2;
    }
}

/**
 * @brief Traceback from W entry.
 * pre: W contains values of row i in interval i..j
 *
 * @param seq Sequence
 * @param structure Final structure
 * @param W W array
 * @param i row index
 * @param j column index
 */
void Spark::trace_W(cand_pos_t i, cand_pos_t j) {
    if (debug) printf("W at %d and %d with %d\n", i, j, W_[j]);

    if (i + TURN + 1 > j) return;
    // case j unpaired
    if (W_[j] == W_[j - 1]) {
        trace_W(i, j-1);
        return;
    }

    cand_pos_t m = j + 1;
    energy_t w = INF;
    for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t wmb_kj = it->second;
        w = W_[m - 1] + wmb_kj;
        if (W_[j] == w + PS_penalty) {
            trace_W(i, m - 1);
            trace_WMB(m, j, wmb_kj);
            return;
        }
    }
    energy_t v = INF;
    w = INF;
    Dangle dangle = 3;
    energy_t vk = INF;
    for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t v_kj = it->fourth >> 2;
        const Dangle d = it->fourth & 3;
        w = W_[m - 1] + v_kj;
        if (W_[j] == w) {
            v = it->second;
            dangle = d;
            vk = v_kj;
            break;
        }
    }
    cand_pos_t k = m;
    cand_pos_t l = j;
    pair_type ptype = 0;
    switch (dangle) {
    case 0:
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_ExtLoop(ptype, -1, -1, params_);
        break;
    case 1:
        k = m + 1;
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_ExtLoop(ptype, S_[m], -1, params_);
        break;
    case 2:
        l = j - 1;
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_ExtLoop(ptype, -1, S_[j], params_);
        break;
    case 3:
        if (params_->model_details.dangles == 1) {
            k = m + 1;
            l = j - 1;
            ptype = pair[S_[k]][S_[l]];
            v = vk - E_ExtLoop(ptype, S_[m], S_[j], params_);
        }
        break;
    }
    assert(i <= m && m < j);
    assert(v < INF);
    // don't recompute W, since i is not changed
    trace_W(i, m - 1);
    trace_V(k, l, v);
}

/**
 * @brief Traceback from V entry
 *
 * @param structure Final Structure
 * @param mark_candidates Whether Candidates should be [ ]
 * @param i row index
 * @param j column index
 */
void Spark::trace_V(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("V at %d and %d with %d\n", i, j, e);

    assert(i + TURN + 1 <= j);
    assert(j <= spark.n_);

    if (mark_candidates_ && is_candidate(CL_, cand_comp, i, j)) {
        structure_[i] = '{';
        structure_[j] = '}';
    } else {
        structure_[i] = '(';
        structure_[j] = ')';
    }
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    if (exists_trace_arrow_from(ta_, i, j)) {

        const TraceArrow &arrow = trace_arrow_from(ta_, i, j);
        const cand_pos_t k = arrow.k(i);
        const cand_pos_t l = arrow.l(j);
        assert(i < k);
        assert(l < j);
        trace_V(k, l, arrow.target_energy());
        return;

    } else {

        // try to trace back to a candidate: (still) interior loop case
        cand_pos_t l_min = std::max(i, j - 31);
        for (cand_pos_t l = j - 1; l >= l_min; --l) {
            // Break if it's an assured dangle case
            for (auto it = CL_[l].begin(); CL_[l].end() != it && it->first > i; ++it) {
                const cand_pos_t k = it->first;
                if (k - i > 31) continue;
                if (e == it->second + E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1], const_cast<vrna_param_t *>(params_))) {
                    trace_V(k, l, it->second);
                    return;
                }
            }
        }
    }
    // is this a hairpin?
    if (e == HairpinE(i, j)) {
        return;
    }

    // if we are still here, trace to wm2 (split case);
    // in this case, we know the 'trace arrow'; the next row has to be recomputed
    std::vector<energy_t> temp;
    if (params_->model_details.dangles == 1) {
        temp = recompute_WM(i + 2, j - 1);
        WM_ = temp;
        dmli1_ = recompute_WM2(i + 2, j - 1);
    }
    WM_  = recompute_WM(i + 1, j - 1);
    WM2_  = recompute_WM2(i + 1, j - 1);

    // Dangle for Multiloop
    cand_pos_t k = i + 1;
    cand_pos_t l = j - 1;
    if (params_->model_details.dangles == 1) {
        find_mb_dangle(WM2_[j - 1], dmli1_[j - 1], WM2_[j - 2], dmli1_[j - 2], i, j, k, l); // Check whether k and l should be used
        if (k > i + 1) {
            WM_.swap(temp);
            WM2_.swap(dmli1_);
        }
    }

    trace_WM2(k,l);
}

/**
 * @brief Traceback from WM
 *
 * @param WM WM array at [i,j]
 * @param WM2 WM2 array at [i,j]
 * @param i row index
 * @param j column index
 * @param e energy in WM(i,j)
 */
void Spark::trace_WM(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("WM at %d and %d with %d\n", i, j, e);

    if (i + TURN + 1 > j) {
        return;
    }

    if (e == WM_[j - 1] + params_->MLbase) {
        trace_WM(i, j - 1, WM_[j - 1]);
        return;
    }
    cand_pos_t m = j + 1;
    for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t wmb_kj = it->second + PSM_penalty;
        energy_t wmb_up = static_cast<energy_t>((m - i) * params_->MLbase) + wmb_kj;
        energy_t wmb_wm = WM_[m - 1] + wmb_kj;
        if (e == wmb_up + PSM_penalty) {
            trace_WMB(m, j, wmb_kj);
            return;
        } else if (e == wmb_wm + PSM_penalty) {
            trace_WM(i, m - 1, WM_[m - 1]);
            trace_WMB(m, j, wmb_kj);
            return;
        }
    }

    energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 3;
    for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t v_kj = it->third >> 2;
        const Dangle d = it->third & 3;
        if (e == WM_[m - 1] + v_kj) {
            dangle = d;
            vk = v_kj;
            v = it->second;
            // no recomp, same i
            break;
        } else if (e == static_cast<energy_t>((m - i) * params_->MLbase) + v_kj) {
            dangle = d;
            vk = v_kj;
            v = it->second;
            break;
        }
    }
    cand_pos_t k = m;
    cand_pos_t l = j;
    pair_type ptype = 0;
    switch (dangle) {
    case 0:
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_MLstem(ptype, -1, -1, params_);
        break;
    case 1:
        k = m + 1;
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_MLstem(ptype, S_[m], -1, params_) - params_->MLbase;
        break;
    case 2:
        l = j - 1;
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_MLstem(ptype, -1, S_[j], params_) - params_->MLbase;
        break;
    case 3:
        if (params_->model_details.dangles == 1) {
            k = m + 1;
            l = j - 1;
            ptype = pair[S_[k]][S_[l]];
            v = vk - E_MLstem(ptype, S_[m], S_[j], params_) - 2 * params_->MLbase;
        }
        break;
    }

    if (e == WM_[m - 1] + vk) {
        // no recomp, same i
        trace_WM(i, m - 1, WM_[m - 1]);
        trace_V(k, l, v);
        return;
    } else if (e == static_cast<energy_t>((m - i) * params_->MLbase) + vk) {
        trace_V(k, l, v);
        return;
    }
    assert(false);
}

/**
 * @brief Traceback from WM2
 *
 * @param WM WM array at [i,j]
 * @param WM2 Wm2 array at [i,j]
 * @param i row index
 * @param j column index
 */
void Spark::trace_WM2(cand_pos_t i, cand_pos_t j) {
    if (debug) printf("WM2 at %d and %d with %d\n", i, j, WM2_[j]);

    if (i + 2 * TURN + 3 > j) {
        return;
    }
    const energy_t e = WM2_[j];

    // case j unpaired
    if (e == WM2_[j - 1] + params_->MLbase) {
        // same i, no recomputation
        trace_WM2(i, j - 1);
        return;
    }

    cand_pos_t m = j + 1;
    for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t wmb_kj = it->second + PSM_penalty;
        energy_t wmb_up = static_cast<energy_t>((m - i) * params_->MLbase) + wmb_kj;
        energy_t wmb_wm = WM_[m - 1] + wmb_kj;
        if (e == wmb_up) {
            trace_WMB(m, j, wmb_kj);
            return;
        } else if (e == wmb_wm) {
            trace_WM(i, m - 1, WM_[m - 1]);
            trace_WMB(m, j, wmb_kj);
            return;
        }
    }

    energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 4;
    for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i + TURN + 1; ++it) {
        m = it->first;

        const energy_t v_kj = it->third >> 2;
        const Dangle d = it->third & 3;
        if (e == WM_[m - 1] + v_kj) {
            vk = v_kj;
            dangle = d;
            v = it->second;
            break;
        }
    }
    cand_pos_t k = m;
    cand_pos_t l = j;
    pair_type ptype = 0;
    switch (dangle) {
    case 0:
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_MLstem(ptype, -1, -1, params_);
        break;
    case 1:
        k = m + 1;
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_MLstem(ptype, S_[m], -1, params_) - params_->MLbase;
        break;
    case 2:
        l = j - 1;
        ptype = pair[S_[k]][S_[l]];
        v = vk - E_MLstem(ptype, -1, S_[j], params_) - params_->MLbase;
        break;
    case 3:
        if (params_->model_details.dangles == 1) {
            k = m + 1;
            l = j - 1;
            ptype = pair[S_[k]][S_[l]];
            v = vk - E_MLstem(ptype, S_[m], S_[j], params_) - 2 * params_->MLbase;
        }
        break;
    }

    if (e == WM_[m - 1] + vk) {
        trace_WM(i, m - 1, WM_[m - 1]);
        trace_V(k, l, v);
        return;
    }
    assert(false);
}

/**
 * @brief Traceback from WMB
 *
 * @param CLBEO BE candidate list outer->in
 * @param i row index
 * @param j column index
 */
void Spark::trace_WMB(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("WMB at i is %d and j is %d and e is %d\n", i, j, e);
    assert(i + TURN + 1 <= j);
    assert(j <= spark.n_);

    recompute_WMBP(i, j);

    cand_pos_t bp_j = tree->tree[j].pair;

    if (tree->tree[j].pair >= 0 && j > tree->tree[j].pair && tree->tree[j].pair > i) {
        energy_t en = INF;
        energy_t BE_energy = INF;
        cand_pos_t best_border = j;
        compute_WMB_case1(j, en, BE_energy, best_border);
        if (e == en + PB_penalty) {
            trace_BE(bp_j, tree->tree[best_border].pair, BE_energy);
            trace_WMBA(i, best_border - 1, en - BE_energy);
        }

        return;
    }
    trace_WMBP(i, j, WMBP_[j]);
    return;
}

/**
 * @brief Traceback from VP entry
 *
 * @param structure Final Structure
 * @param taVP VP trace arrows
 * @param CLVP VP candidate list
 * @param i row index
 * @param j column index
 */
void Spark::trace_VP(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("VP at %d and %d with %d\n", i, j, e);
    structure_[i] = '[';
    structure_[j] = ']';
    if (e == 0) return;

    const pair_type ptype_closing = pair[S_[i]][S_[j]];

    cand_pos_t B_ij = tree->B(i, j);
    cand_pos_t Bp_ij = tree->Bp(i, j);
    cand_pos_t b_ij = tree->b(i, j);
    cand_pos_t bp_ij = tree->bp(i, j);
    if (tree->tree[i].parent->index > 0 && tree->tree[j].parent->index < tree->tree[i].parent->index && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0) {
        recompute_WI(i + 1, Bp_ij - 1);
        recompute_WI(B_ij + 1, j - 1);
        if (e == WI_[Bp_ij - 1] + WI_[j - 1]) {
            trace_WI(i + 1, Bp_ij - 1, WI_[Bp_ij - 1]);
            trace_WI(B_ij + 1, j - 1, WI_[j - 1]);
            return;
        }
    }
    if (tree->tree[i].parent->index < tree->tree[j].parent->index && tree->tree[j].parent->index > 0 && b_ij >= 0 && bp_ij >= 0 && Bp_ij < 0) {
        recompute_WI(i + 1, b_ij - 1);
        recompute_WI(bp_ij + 1, j - 1);
        if (e == (WI_[b_ij - 1] + WI_[j - 1])) {
            trace_WI(i + 1, b_ij - 1, WI_[b_ij - 1]);
            trace_WI(bp_ij + 1, j - 1, WI_[j - 1]);
            return;
        }
    }
    if (tree->tree[i].parent->index > 0 && tree->tree[j].parent->index > 0 && Bp_ij >= 0 && B_ij >= 0 && b_ij >= 0 && bp_ij >= 0) {
        recompute_WI(i + 1, Bp_ij - 1);
        recompute_WI(B_ij + 1, b_ij - 1);
        recompute_WI(bp_ij + 1, j - 1);

        if (e == WI_[Bp_ij - 1] + WI_[b_ij - 1] + WI_[j - 1]) {
            trace_WI(i + 1, Bp_ij - 1, WI_[Bp_ij + 1]);
            trace_WI(B_ij + 1, b_ij - 1, WI_[b_ij - 1]);
            trace_WI(bp_ij + 1, j - 1, WI_[j - 1]);
            return;
        }
    }
    if (exists_trace_arrow_from(taVP_, i, j)) {

        const TraceArrow &arrow = trace_arrow_from(taVP_, i, j);
        const size_t k = arrow.k(i);
        const size_t l = arrow.l(j);
        assert(i < k);
        assert(l < j);
        trace_VP(k, l, arrow.target_energy());
        return;

    } else {

        // try to trace back to a candidate: (still) interior loop case
        cand_pos_t l_min = std::max(i, j - 31);
        for (cand_pos_t l = j - 1; l >= l_min; l--) {
            // Break if it's an assured dangle case
            for (auto it = CLVP_[l].begin(); CLVP_[l].end() != it && it->first > i; ++it) {
                const cand_pos_t k = it->first;

                if (k - i > 31) continue;
                energy_t temp = lrint(((j - l == 1 && k - i == 1) ? e_stP_penalty : e_intP_penalty)
                                      * E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1], const_cast<vrna_param_t *>(params_)));
                if (e == it->second + temp) {
                    trace_VP(k, l, it->second);
                    return;
                }
            }
        }
    }
    // 	// If not other cases, must be WV multiloop
    recompute_WVe(i + 1, j - 1);
    recompute_WIP(i + 1, j - 1);
    recompute_WV(i + 1, j - 1);

    trace_WV(i + 1, j - 1, WV_[j - 1]);
}

/**
 * @brief Traceback from WVe entry
 *
 * @param CLVP VP candidate list
 * @param i row index
 * @param j column index
 */
void Spark::trace_WVe(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("WVe at %d and %d with %d\n", i, j, e);

    if (i + TURN + 1 >= j) return;
    if (WVe_[j] == WVe_[j - 1] + cp_penalty) {
        trace_WVe(i, j - 1, WVe_[j - 1]);
        return;
    }
    cand_pos_t bound_left = j;
    if (tree->b(i, j) > 0) bound_left = tree->b(i, j);
    if (tree->Bp(i, j) > 0) bound_left = std::min((cand_pos_tu)bound_left, (cand_pos_tu)tree->Bp(i, j));
    for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (k > bound_left) continue;
        if (e == static_cast<energy_t>(cp_penalty * (k - i)) + it->second) {
            trace_VP(k, j, it->second);
            return;
        }
    }
}

/**
 * @brief Traceback from WV entry
 *
 * @param CL V candidate list
 * @param CLVP VP candidate list
 * @param CLWMB WMB candidate list
 * @param WIP WIP values in region [i,j]
 * @param i row index
 * @param j column index
 */
void Spark::trace_WV(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("WV at %d and %d with %d\n", i, j, e);

    if (WV_[j] == WV_[j - 1] + cp_penalty) {
        trace_WV(i, j - 1, WV_[j - 1]);
        return;
    }
    cand_pos_t bound_left = std::min(tree->b(i, j), tree->Bp(i, j));
    cand_pos_t bound_right = std::min(tree->b(i, j), tree->Bp(i, j));

    for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first > bound_right; ++it) {
        cand_pos_t k = it->first;
        // energy_t wm_v = it->third >> 2;
        energy_t v = it->second;
        // Dangle d = it->third & 3;
        // cand_pos_t num = 0;
        // if(d == 1 || d==2) num = 1;
        // else if(d==3 && params->model_details.dangles == 1) num =2;
        // energy_t fix = num*cp_penalty - num*params->MLbase - b_penalty;

        if (e == WV_[k - 1] + v + bp_penalty) {
            cand_pos_t m = k;
            cand_pos_t l = j;
            // pair_type ptype = 0;
            energy_t v = it->second;
            // switch(d){
            // 	case 0:
            // 		ptype= pair[S[m]][S[l]];
            // 		v = wm_v - E_MLstem(ptype,-1,-1,params);
            // 		break;
            // 	case 1:
            // 		m=k+1;
            // 		ptype= pair[S[m]][S[l]];
            // 		v = wm_v - E_MLstem(ptype,S[k],-1,params) - params->MLbase;
            // 		break;
            // 	case 2:
            // 		l=j-1;
            // 		ptype= pair[S[m]][S[l]];
            // 		v = wm_v - E_MLstem(ptype,-1,S[j],params) - params->MLbase;
            // 		break;
            // 	case 3:
            // 		if(params->model_details.dangles == 1){
            // 			m=k+1;
            // 			l=j-1;
            // 			ptype= pair[S[m]][S[l]];
            // 			v = wm_v - E_MLstem(ptype,S[k],S[j],params) - 2*params->MLbase;
            // 		}
            // 		break;
            // }

            trace_WV(i, k - 1, WV_[k - 1]);
            trace_V( m, l, v);
            return;
        }
        if (e == WVe_[k - 1] + v + bp_penalty) {
            cand_pos_t m = k;
            cand_pos_t l = j;
            // pair_type ptype = 0;
            energy_t v = it->second;

            // switch(d){
            // 	case 0:
            // 		ptype= pair[S[m]][S[l]];
            // 		v = wm_v - E_MLstem(ptype,-1,-1,params);
            // 		break;
            // 	case 1:
            // 		m=k+1;
            // 		ptype= pair[S[m]][S[l]];
            // 		v = wm_v - E_MLstem(ptype,S[k],-1,params) - params->MLbase;
            // 		break;
            // 	case 2:
            // 		l=j-1;
            // 		ptype= pair[S[m]][S[l]];
            // 		v = wm_v - E_MLstem(ptype,-1,S[j],params) - params->MLbase;
            // 		break;
            // 	case 3:
            // 		if(params->model_details.dangles == 1){
            // 			m=k+1;
            // 			l=j-1;
            // 			ptype= pair[S[m]][S[l]];
            // 			v = wm_v - E_MLstem(ptype,S[k],S[j],params) - 2*params->MLbase;
            // 		}
            // 		break;
            // }
            trace_WVe(i, k - 1,WVe_[k - 1]);
            trace_V(m, l, v);
            return;
        }
    }

    for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first > bound_right; ++it) {
        cand_pos_t k = it->first;
        if (e == WV_[k - 1] + it->second + PSM_penalty + bp_penalty) {
            trace_WV(i, k - 1, WV_[k - 1]);
            trace_WMB(k, j, it->second);
            return;
        }
        if (e == WVe_[k - 1] + it->second + PSM_penalty + bp_penalty) {
            trace_WVe(i, k - 1, WV_[k - 1]);
            trace_WMB(k, j, it->second);
            return;
        }
    }
    for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (k > bound_left) continue;
        if (e == WIP_[k - 1] + it->second) {
            trace_WIP(i, k - 1, WIP_[k - 1]);
            trace_VP(k, j, it->second);
            return;
        }
    }
}

/**
 * @brief Traceback from WI entry
 *
 * @param CL V candidate list
 * @param CLWMB WMB candidate list
 * @param WI WI values in region [i,j]
 * @param i row index
 * @param j column index
 */
void Spark::trace_WI(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("WI at %d and %d with %d\n", i, j, e);

    if (i + TURN + 1 >= j) return;

    // How to do one base backwards?

    if (e == WI_[j - 1] + PUP_penalty) {
        trace_WI(i, j - 1, WI_[j - 1]);
        return;
    }

    for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == WI_[k - 1] + it->second + PPS_penalty) {
            trace_WI(i, k - 1, WI_[k - 1]);
            trace_V(k, j, it->second);
            return;
        }
    }

    for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == WI_[k - 1] + it->second + PPS_penalty + PSM_penalty) {
            trace_WI(i, k - 1, WI_[k - 1]);
            trace_WMB(k, j, it->second);
            return;
        }
    }
    assert(false);
}

/**
 * @brief Traceback from WIP entry
 *
 * @param CL V candidate list
 * @param CLVP VP candidate list
 * @param WIP WIP values in region [i,j]
 * @param i row index
 * @param j column index
 */
void Spark::trace_WIP(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("WIP at %d and %d with %d\n", i, j, e);

    if (i + TURN + 1 >= j) return;

    // 	// How to do one base backwards?

    if (e == WIP_[j - 1] + cp_penalty) {
        trace_WIP(i, j - 1, WIP_[j - 1]);
        return;
    }
    for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == WIP_[k - 1] + it->second + bp_penalty) {
            trace_WIP(i, k - 1, WIP_[k - 1]);
            trace_V(k, j, it->second);
            return;
        }
        if (e == static_cast<energy_t>((k - i) * cp_penalty) + it->second + bp_penalty) {
            trace_V(k, j, it->second);
            return;
        }
    }

    for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == WIP_[k - 1] + it->second + bp_penalty + PSM_penalty) {
            // Why do I pick two different variables for WIP
            trace_WIP(i, k - 1, WIP_Bbp[k - 1]);
            trace_WMB(k, j, it->second);
            return;
        }
        if (e == static_cast<energy_t>((k - i) * cp_penalty) + it->second + bp_penalty + PSM_penalty) {
            trace_WMB(k, j, it->second);
            return;
        }
    }
}

/**
 * @brief Traceback from WMBP entry
 *
 * @param CL V candidate list
 * @param CLVP VP candidate list
 * @param CLWMB WMB candidate list
 * @param CLBE BE candidate list
 * @param WI_Bbp WI values in regions where the left index is a outer closing base or an inner opening base
 * @param i row index
 * @param j column index
 */
void Spark::trace_WMBP(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("WMBP at %d and %d with %d\n", i, j, e);

    energy_t VP_ij = INF;

    if (tree->tree[j].pair < 0) {
        cand_pos_t b_ij = tree->b(i, j);
        for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it && it->first >= i; ++it) {
            cand_pos_t k = it->first;
            if (k == i) VP_ij = it->second; // second?

            int ext_case = compute_exterior_cases(k, j);
            if ((b_ij > 0 && k < b_ij) || (b_ij < 0 && ext_case == 0)) {
                cand_pos_t bp_ik = tree->bp(i, k);
                cand_pos_t Bp_kj = tree->Bp(k, j);
                if (bp_ik >= 0 && k > bp_ik && Bp_kj > 0 && k < Bp_kj) {
                    energy_t BE_energy = INF;
                    cand_pos_t B_kj = tree->B(k, j);
                    cand_pos_t bp_kj = tree->tree[Bp_kj].pair;
                    cand_pos_t b_kj = (B_kj > 0) ? tree->tree[B_kj].pair : -2;
                    for (auto it2 = CLBE_[Bp_kj].begin(); CLBE_[Bp_kj].end() != it2; ++it2) {
                        cand_pos_t l = it2->first;
                        if (l == b_kj) {
                            BE_energy = it2->second;
                            break;
                        }
                    }
                    if (e == WMBA_[k - 1] + it->second + 2 * PB_penalty + BE_energy) {
                        trace_BE(b_kj, bp_kj, BE_energy);
                        trace_WMBA(i, k - 1, WMBP_[k - 1]);
                        trace_VP(k, j, it->second);
                        return;
                    }
                }
            }
        }
    }

    if (tree->tree[j].pair < 0 && tree->tree[i].pair >= 0 && tree->tree[i].pair < j) {
        for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it && it->first >= i; ++it) {
            cand_pos_t k = it->first;
            cand_pos_t bp_ik = tree->bp(i, k);
            if (bp_ik >= 0 && k + TURN <= j) {
                cand_pos_t Bp_ik = tree->tree[bp_ik].pair;
                energy_t BE_energy = INF;
                for (auto it2 = CLBE_[Bp_ik].begin(); CLBE_[Bp_ik].end() != it2; ++it2) {
                    cand_pos_t l = it2->first;
                    if (l == i) {
                        BE_energy = it2->second;
                        break;
                    }
                }

                if (e == WI_Bbp[k - 1] + it->second + 2 * PB_penalty + BE_energy) {
                    recompute_WI(bp_ik + 1, k - 1);
                    trace_BE(i, bp_ik, BE_energy);
                    trace_WI(bp_ik + 1, k - 1, WI_[k - 1]);
                    trace_VP(k, j, it->second);
                    return;
                }
            }
        }
    }

    if (e == VP_ij + PB_penalty) {
        trace_VP(i, j, VP_ij);
        return;
    }
}

/**
 * @brief Traceback from WMBA entry
 *
 * @param CL V candidate list
 * @param CLWMB WMB candidate list
 * @param i row index
 * @param j column index
 */
void Spark::trace_WMBA(cand_pos_t i, cand_pos_t j, energy_t e) {
    if (debug) printf("WMBA at %d and %d with %d\n", i, j, e);

    if (WMBA_[j] == WMBA_[j - 1] + PUP_penalty) {
        trace_WMBA(i, j - 1, WMBA_[j - 1]);
        return;
    }

    if (WMBA_[j] == WMBP_[j]) {
        trace_WMBP(i, j, WMBP_[j]);
        return;
    }

    for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == WMBA_[k - 1] + it->second + PPS_penalty) {
            trace_WMBA(i, k - 1, WMBP_[k - 1]);
            trace_V(k, j, it->second);
            return;
        }
    }

    for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == WMBA_[k - 1] + it->second + PPS_penalty + PSM_penalty) {
            trace_WMBA(i, k - 1, WMBP_[k - 1]);
            trace_WMB(k, j, it->second);
            return;
        }
    }
}

/**
 * @brief Traceback from BE entry
 *
 * @param CLBE BE candidate list
 * @param WIP WIP values in region [i,j]
 * @param i row index
 * @param j column index
 */
void Spark::trace_BE(cand_pos_t i, cand_pos_t ip, energy_t e) {
    cand_pos_t j = tree->tree[i].pair;
    cand_pos_t jp = tree->tree[ip].pair;
    cand_pos_t lp = jp;
    // currently this is i lp ip  l j and it should be i lp jp j l ip

    if (debug) printf("BE at [%d,%d] U [%d,%d] with %d\n", i, ip, jp, j, e);

    structure_[i] = '(';
    structure_[j] = ')';
    const pair_type ptype_closing_ij = pair[S_[i]][S_[j]];
    if (i == ip) return;

    energy_t BE_energy = INF;

    for (auto it = CLBE_[jp].begin(); CLBE_[jp].end() != it && it->first > i; ++it) {
        lp = it->first;
        BE_energy = it->second;
    }
    cand_pos_t l = tree->tree[lp].pair;
    // i lp ip       jp  l   j

    if (e == lrint(e_stP_penalty * ILoopE(ptype_closing_ij,i,j,lp,l)) + BE_energy) {
        trace_BE(lp, ip, BE_energy);
        return;
    }
    if (e == lrint(e_intP_penalty * ILoopE(ptype_closing_ij,i,j,lp,l)) + BE_energy) {
        trace_BE(lp, ip, BE_energy);
        return;
    }

    if (e == WIP_Bbp[lp - 1] + BE_energy + WIP_Bp[l + 1] + ap_penalty + 2 * bp_penalty) {
        recompute_WIP(i + 1, lp - 1);
        recompute_WIP(l + 1, j - 1);
        trace_WIP(i + 1, lp - 1, WIP_Bbp[lp - 1]);
        trace_BE(lp, ip, BE_energy);
        trace_WIP(l + 1, j - 1, WIP_Bp[j - 1]);
        return;
    }
    if (e == cp_penalty * ((lp - i - 1)) + BE_energy + WIP_Bp[l + 1] + ap_penalty + 2 * bp_penalty) {
        recompute_WIP(l + 1, j - 1);
        trace_BE(lp, ip, BE_energy);
        trace_WIP(l + 1, j - 1, WIP_Bbp[j - 1]);
        return;
    }
    if (e == WIP_Bbp[lp - 1] + BE_energy + cp_penalty * ((j - l - 1)) + ap_penalty + 2 * bp_penalty) {
        recompute_WIP(i + 1, lp - 1);
        trace_WIP(i + 1, lp - 1, WIP_Bbp[lp - 1]);
        trace_BE(lp, ip, BE_energy);
        return;
    }
}
/**
 * @brief Trace back
 * pre: row 1 of matrix W is computed
 * @return mfe structure (reference)
 */
const std::string& Spark::trace_back() {

    structure_.resize(n_ + 1, '.');

    /* Traceback */
    trace_W(1, n_);
    structure_ = structure_.substr(1, n_);

    return structure_;
}

/**
 * @brief Computes the values for the WVe matrix
 * @param i start
 * @param j end
 * @param spark datastructure
 * @param tree tree for boundary determination
 */
energy_t Spark::compute_WVe(cand_pos_t i, cand_pos_t j) {
    energy_t wve = INF;
    cand_pos_t bound = std::min((cand_pos_tu)tree->Bp(i, j), (cand_pos_tu)tree->b(i, j));
    for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it; ++it) {
        cand_pos_t k = it->first;
        bool can_pair = tree->up[k - 1] >= (k - i);
        if (can_pair && k < bound) wve = std::min(wve, static_cast<energy_t>(cp_penalty * (k - i)) + it->second);
    }
    if (tree->tree[j].pair < 0) wve = std::min(wve, WVe_[j - 1] + cp_penalty);
    return wve;
}

/**
 * @brief Computes the values for the WV matrix
 * @param j end
 * @param bound_left boundary for non-closed region
 * @param bound_right boundary for closed region
 * @param spark datastructure
 * @param tree tree for boundary determination
 */
energy_t Spark::compute_WV(cand_pos_t j, cand_pos_t bound_left, cand_pos_t bound_right) {
    energy_t m1 = INF, m2 = INF, m3 = INF, m5 = INF, m6 = INF, wv = INF;

    for (auto it = CL_[j].begin(); CL_[j].end() != it && it->first > bound_right; ++it) {
        cand_pos_t k = it->first;
        // energy_t val = it->third >> 2;
        energy_t val = it->second;
        // Dangle d = it->third & 3;
        // cand_pos_t num = 0;
        // if(d == 1 || d==2) num = 1;
        // else if(d==3 && params->model_details.dangles == 1) num =2;
        // energy_t fix = num*cp_penalty - num*params->MLbase-b_penalty;
        // m1 = std::min(m1, WVe[k-1] + val + fix + bp_penalty);
        // m5 = std::min(m5, WV[k-1] + val+ fix + bp_penalty);
        m1 = std::min(m1, WVe_[k - 1] + val + bp_penalty);
        m5 = std::min(m5, WV_[k - 1] + val + bp_penalty);
    }
    for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it && it->first > bound_right; ++it) {
        cand_pos_t k = it->first;
        energy_t val = it->second;
        m2 = std::min(m1, WVe_[k - 1] + val + PSM_penalty + bp_penalty);
        m6 = std::min(m6, WV_[k - 1] + val + PSM_penalty + bp_penalty);
    }

    for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it; ++it) {
        cand_pos_t k = it->first;
        energy_t val = it->second;
        if (k < bound_left) m3 = std::min(m1, dwip1_[k - 1] + val);
    }
    wv = std::min({m1, m2, m3, m5, m6});

    if (tree->tree[j].pair < 0) wv = std::min(wv, WV_[j - 1] + cp_penalty);

    return wv;
}

/**
 * @brief computes entries for BE
 * @param i outer opening base pair
 * @param j inner closing base pair
 * @param ip outer closing base pair
 * @param jp inner opening base pair
 * @param spark datastructure
 * @param tree tree for boundary estimation
 */
energy_t Spark::compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp) {
    // We are checking for the closest pair that we have already calculated to i/ip from j/jp
    // If there is nothing, then i is the closest encompassing pair to jp
    // If it is not, then we get the energy for everything from jp to lp so that we calculate less

    // (.....(..(....)..).....)
    // i     lp      j  l     ip
    energy_t BE_energy = INF;
    cand_pos_t lp = jp;
    if (!CLBE_[j].empty()) {
        auto const [k, vbe] = CLBE_[j].back();
        BE_energy = vbe;
        lp = k;
    }
    cand_pos_t l = tree->tree[lp].pair; // right closing base for lp

    const pair_type ptype_closing_iip = pair[S_[i]][S_[ip]];

    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF, val = INF;
    // 1
    if (i + 1 == lp && ip - 1 == l) {
        m1 = lrint(e_stP_penalty * ILoopE(ptype_closing_iip, i, ip, lp, l)) + BE_energy;
        val = std::min(val, m1);
    }

    bool empty_region_ilp = (tree->up[lp - 1] >= lp - i - 1);    // empty between i+1 and lp-1
    bool empty_region_lip = (tree->up[ip - 1] >= ip - l - 1);    // empty between l+1 and ip-1
    bool weakly_closed_ilp = tree->weakly_closed(i + 1, lp - 1); // weakly closed between i+1 and lp-1
    bool weakly_closed_lip = tree->weakly_closed(l + 1, ip - 1); // weakly closed between l+1 and ip-1

    // 2
    if (empty_region_ilp && empty_region_lip) {
        m2 = lrint(e_intP_penalty * ILoopE(ptype_closing_iip, i, ip, lp, l)) + BE_energy;
        val = std::min(val, m2);
    }

    // 3
    if (weakly_closed_ilp && weakly_closed_lip) {
        m3 = dwip1_[lp - 1] + BE_energy + WIP_Bp[l + 1] + ap_penalty + 2 * bp_penalty;
        val = std::min(val, m3);
    }

    // 4
    if (weakly_closed_ilp && empty_region_lip) {
        m4 = dwip1_[lp - 1] + BE_energy + cp_penalty * (ip - l - 1) + ap_penalty + 2 * bp_penalty;
        val = std::min(val, m4);
    }

    // 5
    if (empty_region_ilp && weakly_closed_lip) {

        m5 = ap_penalty + 2 * bp_penalty + (cp_penalty * (lp - i - 1)) + BE_energy + WIP_Bp[l + 1];
        val = std::min(val, m5);
    }

    return val;
}
/**
 * @brief Computes entries for WMBP
 * @param i start
 * @param j end
 * @param spark datastructure
 * @param tree tree for boundary estimation
 */
energy_t Spark::compute_WMBP(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF, m2 = INF, m3 = INF, wmbp = INF;
    // 1) WMBP(i,j) = BE(bpg(Bp(l,j)),Bp(l,j),bpg(B(l,j)),B(l,j)) + WMBP(i,l) + VP(l+1,j)
    if (tree->tree[j].pair < 0) {
        energy_t tmp = INF;
        cand_pos_t b_ij = tree->b(i, j);
        for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it; ++it) {
            cand_pos_t k = it->first;
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(k, j);
            if ((b_ij > 0 && k < b_ij) || (b_ij < 0 && ext_case == 0)) {
                cand_pos_t bp_ik = tree->bp(i, k);
                cand_pos_t Bp_kj = tree->Bp(k, j);
                if (bp_ik >= 0 && k > bp_ik && Bp_kj > 0 && k < Bp_kj) { // if(sparse_tree->b(i,j)>=0 && l <sparse_tree->b(i,j)){//
                    cand_pos_t B_kj = tree->B(k, j);
                    if (i <= tree->tree[k].parent->index && tree->tree[k].parent->index < j && k + 3 <= j) {
                        energy_t BE_energy = INF;
                        cand_pos_t b_kj = tree->tree[B_kj].pair;
                        for (auto it2 = CLBE_[Bp_kj].begin(); CLBE_[Bp_kj].end() != it2; ++it2) {
                            cand_pos_t l = it2->first;
                            if (l == b_kj) {
                                BE_energy = it2->second;
                                break;
                            }
                        }
                        energy_t WMBA_energy = WMBA_[k - 1];
                        energy_t VP_energy = it->second;
                        energy_t sum = BE_energy + WMBA_energy + VP_energy;

                        tmp = std::min(tmp, sum);
                    }
                }
            }

            m1 = 2 * PB_penalty + tmp;
        }
    }

    // 2) WMBP(i,j) = VP(i,j) + P_b
    cand_pos_t i_mod = i % (MAXLOOP + 1);
    m2 = VP_(i_mod, j) + PB_penalty;

    // check later if <0 or <-1

    // WMBP(i,j) = BE(i,,,) _ WI(bp(i,k),k-1) + VP(k,j)
    if (tree->tree[j].pair < 0 && tree->tree[i].pair >= 0 && tree->tree[i].pair < j) {
        energy_t tmp = INF;
        for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it; ++it) {
            cand_pos_t k = it->first;
            cand_pos_t bp_ik = tree->bp(i, k);
            if (bp_ik >= 0 && k + TURN <= j) {
                energy_t BE_energy = INF;
                cand_pos_t Bp_ik = tree->tree[bp_ik].pair;
                if (!CLBE_[Bp_ik].empty()) {
                    auto const [l, vbe] = CLBE_[Bp_ik].back();
                    if (i == l) BE_energy = vbe;
                }

                energy_t WI_energy = (k - 1 - (bp_ik + 1)) > 4 ? WI_Bbp[k - 1] : PUP_penalty * (k - 1 - (bp_ik + 1) + 1);
                energy_t VP_energy = it->second;
                energy_t sum = BE_energy + WI_energy + VP_energy;

                tmp = std::min(tmp, sum);
            }
        }

        m3 = 2 * PB_penalty + tmp;
    }
    // get the min for WMB
    wmbp = std::min({m1, m2, m3});

    return (wmbp);
}
/**
 * @brief Computes entries for WMBA
 * @param j end
 * @param spark datastructure
 * @param tree tree for boundary estimation
 */
energy_t Spark::compute_WMBA(cand_pos_t j) {

    // WMBA criteria
    energy_t wmba = INF;
    if (tree->tree[j].parent->index > 0) {

        for (auto it = CL_[j].begin(); CL_[j].end() != it; ++it) {
            cand_pos_t k = it->first;
            if (tree->tree[k].pair < 0 && tree->tree[k].parent->index > -1 && tree->tree[j].parent->index > -1
                && tree->tree[j].parent->index == tree->tree[k].parent->index) {
                wmba = std::min(wmba, WMBA_[k - 1] + it->second + PPS_penalty);
            }
        }

        for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it; ++it) {
            cand_pos_t k = it->first;
            if (tree->tree[k].pair < 0 && tree->tree[k].parent->index > -1 && tree->tree[j].parent->index > -1
                && tree->tree[j].parent->index == tree->tree[k].parent->index) {
                wmba = std::min(wmba, WMBA_[k - 1] + it->second + PPS_penalty + PSM_penalty);
            }
        }
        if (tree->tree[j].pair < 0) wmba = std::min(wmba, WMBA_[j - 1] + PUP_penalty);
    } else {
        wmba = INF;
    }
    wmba = std::min(wmba, WMBP_[j]);
    return wmba;
}
/**
 * @brief Computes entries for WMB
 * @param i start
 * @param j end
 * @param spark datastructure
 * @param tree tree for boundary estimation
 */
energy_t Spark::compute_WMB(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF, m2 = INF, wmb = INF;

    if (tree->tree[j].pair >= 0 && j > tree->tree[j].pair && tree->tree[j].pair > i) {
        cand_pos_t best_border = j - 1;
        energy_t BE_energy = INF;
        compute_WMB_case1(j, m1, BE_energy, best_border);
        m1 = PB_penalty + m1;
    }
    // check the WMBP_ij value
    m2 = WMBP_[j];

    wmb = std::min(m1, m2);
    return wmb;
}

/**
 * @brief Computes VP internal loop value for case 5 of VP
 *
 * @param VP VP array
 * @param i row index
 * @param j column index
 */
energy_t Spark::compute_VP_internal(cand_pos_t i, cand_pos_t j, cand_pos_t b_ij, cand_pos_t bp_ij, cand_pos_t Bp_ij, cand_pos_t B_ij, cand_pos_t &best_k, cand_pos_t &best_l, energy_t &best_e){

    energy_t m5 = INF;
    // By doing uint, we make sure it can't be negative
    cand_pos_t min_borders = std::min((cand_pos_tu)Bp_ij, (cand_pos_tu)b_ij);
    cand_pos_t edge_i = std::min(i + MAXLOOP + 1, j - TURN - 1);
    min_borders = std::min({min_borders, edge_i});
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    for (cand_pos_t k = i + 1; k <= min_borders; k++) {
        cand_pos_t k_mod = k % (MAXLOOP + 1);

        energy_t cank = ((tree->up[k - 1] >= (k - i - 1)) - 1);
        cand_pos_t max_borders = std::max(bp_ij, B_ij) + 1;
        cand_pos_t edge_j = k + j - i - MAXLOOP - 2;
        max_borders = std::max({max_borders, edge_j});

        for (cand_pos_t l = j - 1; l >= max_borders; --l) {
            assert(k - i + j - l - 2 <= MAXLOOP);

            energy_t canl = (((tree->up[j - 1] >= (j - l - 1)) - 1) | cank);
            energy_t v_iloop_kl = INF & canl;
            v_iloop_kl = v_iloop_kl + VP_(k_mod, l) + lrint(e_intP_penalty
                                 * E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1], const_cast<vrna_param_t *>(params_)));
            if (v_iloop_kl < m5) {
                m5 = v_iloop_kl;
                best_l = l;
                best_k = k;
                best_e = VP_(k_mod, l);
            }
        }
    }
    return m5;
}

/**
 * @brief Computes the energy values for VP
 * @param i row index
 * @param j column index
 * @param b_ij left outer boundary
 * @param bp_ij left inner boundary
 * @param BP_ij right inner boundary
 * @param B_ij right outer boundary
 * @param spark datastructure
 * @param tree tree for boundary determination
 */
energy_t Spark::compute_VP(cand_pos_t i, cand_pos_t j, cand_pos_t b_ij, cand_pos_t bp_ij, cand_pos_t Bp_ij, cand_pos_t B_ij) {
    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF, m6 = INF, vp = INF;
    const pair_type ptype_closing = pair[S_[i]][S_[j]];

    if (tree->tree[i].parent->index > 0 && tree->tree[j].parent->index < tree->tree[i].parent->index && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0) {
        energy_t WI_ipus1_BPminus = dwi1_[Bp_ij - 1];
        energy_t WI_Bplus_jminus = (j - 1 - (B_ij + 1)) > 4 ? WI_Bbp[j - 1] : PUP_penalty * (j - 1 - (B_ij + 1) + 1);

        m1 = WI_ipus1_BPminus + WI_Bplus_jminus;
    }
    if (tree->tree[i].parent->index < tree->tree[j].parent->index && tree->tree[j].parent->index > 0 && b_ij >= 0 && bp_ij >= 0 && Bp_ij < 0) {
        energy_t WI_i_plus_b_minus = dwi1_[b_ij - 1];
        energy_t WI_bp_plus_j_minus = (j - 1 - (bp_ij + 1)) > 4 ? WI_Bbp[j - 1] : PUP_penalty * (j - 1 - (bp_ij + 1) + 1);

        m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
    }

    if (tree->tree[i].parent->index > 0 && tree->tree[j].parent->index > 0 && Bp_ij >= 0 && B_ij >= 0 && b_ij >= 0 && bp_ij >= 0) {
        energy_t WI_i_plus_Bp_minus = dwi1_[Bp_ij - 1];
        energy_t WI_B_plus_b_minus = (b_ij - 1 - (B_ij + 1)) > 4 ? WI_Bbp[b_ij - 1] : PUP_penalty * (b_ij - 1 - (B_ij + 1) + 1);
        energy_t WI_bp_plus_j_minus = (j - 1 - (bp_ij + 1)) > 4 ? WI_Bbp[j - 1] : PUP_penalty * (j - 1 - (bp_ij + 1) + 1);

        m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
    }
    if (tree->tree[i + 1].pair < -1 && tree->tree[j - 1].pair < -1) {
        cand_pos_t ip1_mod = (i + 1) % (MAXLOOP + 1);

        m4 = lrint(e_stP_penalty * ILoopE(ptype_closing, i, j, i + 1, j - 1)) + VP_(ip1_mod, j - 1);
    }

    cand_pos_t best_k = 0;
    cand_pos_t best_l = 0;
    energy_t best_e = 0;

    m5 = compute_VP_internal(i, j, b_ij, bp_ij, Bp_ij, B_ij, best_k, best_l, best_e);

    // case 6 and 7
    m6 = dwvp_[j - 1] + ap_penalty + 2 * bp_penalty;

    energy_t vp_h = std::min({m1, m2, m3});
    energy_t vp_iloop = std::min(m4, m5);
    if (m4 < m5) {
        best_k = i + 1;
        best_l = j - 1;
        cand_pos_t ip1_mod = (i + 1) % (MAXLOOP + 1);
        best_e = VP_(ip1_mod, j - 1);
    }
    energy_t vp_split = m6;
    vp = std::min({vp_h, vp_iloop, vp_split});

    if (vp_iloop < std::min(vp_h, vp_split)) {
        if (is_candidate(CLVP_, cand_comp, best_k, best_l)) {
            avoid_trace_arrow(taVP_);
        } else {
            register_trace_arrow(taVP_, i, j, best_k, best_l, best_e);
        }
    }

    return vp;
}

/**
 * @brief Computes the internal loop value for V
 *
 * @param V V array
 * @param i row index
 * @param j column index
 */
energy_t Spark::compute_internal(cand_pos_t i, cand_pos_t j, cand_pos_t &best_k, cand_pos_t &best_l, energy_t &best_e){
    energy_t v_iloop = INF;
    cand_pos_t max_k = std::min(j - TURN - 2, i + MAXLOOP + 1);
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    for (cand_pos_t k = i + 1; k <= max_k; k++) {
        cand_pos_t k_mod = k % (MAXLOOP + 1);

        energy_t cank = ((tree->up[k - 1] >= (k - i - 1)) - 1);
        cand_pos_t min_l = std::max(k + TURN + 1 + MAXLOOP + 2, k + j - i) - MAXLOOP - 2;
        // cand_pos_t ind = k_mod*V.ydim_;
        for (cand_pos_t l = j - 1; l >= min_l; --l) {
            assert(k - i + j - l - 2 <= MAXLOOP);
            energy_t canl = (((tree->up[j - 1] >= (j - l - 1)) - 1) | cank);
            energy_t v_iloop_kl = INF & canl;

            v_iloop_kl = v_iloop_kl + V_(k_mod,l)
                         + E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1], const_cast<vrna_param_t *>(params_));
            if (v_iloop_kl < v_iloop) {
                v_iloop = v_iloop_kl;
                best_l = l;
                best_k = k;
                best_e = V_(k_mod, l);
            }
        }
    }
    return v_iloop;
}

/**
 * @brief Determines the MFE energy for a given sequence
 * @param spark datastructure
 * @param tree tree for boundaary determination
 * @param n length of sequence
 * @param garbage_collect whether to garbage collect for trace_arrows
 */
energy_t Spark::fold(){
    Dangle d = 3;
    if (params_->model_details.dangles == 0 || params_->model_details.dangles == 1) d = 0;

    for (cand_pos_t i = n_; i > 0; --i) {
        if (pseudoknot) {
            for (cand_pos_t j = i; j <= n_ && tree->tree[j].pair < 0; ++j) {
                WI_[j] = (j - i + 1) * PUP_penalty;
            }
        }

        for (cand_pos_t j = i + TURN + 1; j <= n_; j++) {

            bool evaluate = tree->weakly_closed(i, j);
            // ------------------------------
            // W: split case
            bool pairedkj = 0;
            energy_t w_split = INF;
            energy_t wm_split = INF;
            energy_t wm2_split = INF;
            energy_t wi_split = INF;
            energy_t wip_split = INF;
            for (auto it = CL_[j].begin(); CL_[j].end() != it; ++it) {
                cand_pos_t k = it->first;

                const energy_t v_kj = it->third >> 2;
                const energy_t v_kjw = it->fourth >> 2;
                bool can_pair = tree->up[k - 1] >= (k - i);
                // WM Portion
                wm_split = std::min(wm_split, WM_[k - 1] + v_kj);
                if (can_pair) wm_split = std::min(wm_split, static_cast<energy_t>((k - i) * params_->MLbase) + v_kj);
                // WM2 Portion
                wm2_split = std::min(wm2_split, WM_[k - 1] + v_kj);
                // W Portion
                w_split = std::min(w_split, W_[k - 1] + v_kjw);

                // WI portion
                energy_t v_kjj = it->second + PPS_penalty;
                wi_split = std::min(wi_split, WI_[k - 1] + v_kjj);
                // WIP portion
                v_kjj = it->second + bp_penalty;
                wip_split = std::min(wip_split, WIP_[k - 1] + v_kjj);
                if (can_pair) wip_split = std::min(wip_split, static_cast<energy_t>((k - i) * cp_penalty) + v_kjj);
            }

            if (tree->weakly_closed(i, j)) {
                for (auto it = CLWMB_[j].begin(); CLWMB_[j].end() != it; ++it) {

                    if (pairedkj) break; // Not needed I believe as there shouldn't be any candidates there if paired anyways
                    // Maybe this would just avoid this loop however

                    cand_pos_t k = it->first;
                    bool can_pair = tree->up[k - 1] >= (k - i);

                    // For W
                    energy_t wmb_kj = it->second + PS_penalty;
                    w_split = std::min(w_split, W_[k - 1] + wmb_kj);
                    // For WM -> I believe this would add a PSM penalty for every pseudoknot which would be bad
                    wmb_kj = it->second + PSM_penalty + b_penalty;
                    wm_split = std::min(wm_split, WM_[k - 1] + wmb_kj);
                    if (can_pair) wm_split = std::min(wm_split, static_cast<energy_t>((k - i) * params_->MLbase) + wmb_kj);
                    wm2_split = std::min(wm2_split, WM_[k - 1] + wmb_kj);
                    if (can_pair) wm2_split = std::min(wm2_split, static_cast<energy_t>((k - i) * params_->MLbase) + wmb_kj);
                    // For WI
                    wmb_kj = it->second + PSM_penalty + PPS_penalty;
                    wi_split = std::min(wi_split, WI_[k - 1] + wmb_kj);

                    // For WIP
                    wmb_kj = it->second + PSM_penalty + bp_penalty;
                    wip_split = std::min(wip_split, WIP_[k - 1] + wmb_kj);
                    if (can_pair) wip_split = std::min(wip_split, static_cast<energy_t>((k - i) * cp_penalty) + wmb_kj);
                }
            }

            if (tree->tree[j].pair < 0) w_split = std::min(w_split, W_[j - 1]);
            if (tree->tree[j].pair < 0) wm2_split = std::min(wm2_split, WM2_[j - 1] + params_->MLbase);
            if (tree->tree[j].pair < 0) wm_split = std::min(wm_split, WM_[j - 1] + params_->MLbase);
            if (tree->tree[j].pair < 0) wi_split = std::min(wi_split, WI_[j - 1] + PUP_penalty);
            if (tree->tree[j].pair < 0) wip_split = std::min(wip_split, WIP_[j - 1] + cp_penalty);

            energy_t w = w_split;   // entry of W w/o contribution of V
            energy_t wm = wm_split; // entry of WM w/o contribution of V

            size_t i_mod = i % (MAXLOOP + 1);

            const pair_type ptype_closing = pair[S_[i]][S_[j]];
            const bool restricted = tree->tree[i].pair == -1 || tree->tree[j].pair == -1;

            const bool unpaired = (tree->tree[i].pair < -1 && tree->tree[j].pair < -1);
            const bool paired = (tree->tree[i].pair == j && tree->tree[j].pair == i);
            const bool pkonly = (!pk_only || paired);
            energy_t v = INF;
            // ----------------------------------------
            // cases with base pair (i,j)
            // if(ptype_closing>0 && !restricted && evaluate) { // if i,j form a canonical base pair
            if (ptype_closing > 0 && !restricted && evaluate && pkonly) {
                bool canH = (paired || unpaired);
                if (tree->up[j - 1] < (j - i - 1)) canH = false;

                energy_t v_h = canH ? HairpinE(i, j) : INF;
                // info of best interior loop decomposition (if better than hairpin)
                cand_pos_t best_l = 0;
                cand_pos_t best_k = 0;
                energy_t best_e = INF;

                energy_t v_iloop = INF;

                // constraints for interior loops
                // i<k; l<j
                // k-i+j-l-2<=MAXLOOP  ==> k <= MAXLOOP+i+1
                //            ==> l >= k+j-i-MAXLOOP-2
                // l-k>=TURN+1         ==> k <= j-TURN-2
                //            ==> l >= k+TURN+1
                // j-i>=TURN+3
                //
                if ((tree->tree[i].pair < -1 && tree->tree[j].pair < -1) || tree->tree[i].pair == j) {
                    v_iloop = compute_internal(i, j, best_k, best_l, best_e);
                }
                const energy_t v_split = E_MbLoop(dmli1_, dmli2_, i, j);

                v = std::min(v_h, std::min(v_iloop, v_split));
                // register required trace arrows from (i,j)
                if (v_iloop < std::min(v_h, v_split)) {
                    if (is_candidate(CL_, cand_comp, best_k, best_l)) {
                        avoid_trace_arrow(ta_);
                    } else {
                        register_trace_arrow(ta_, i, j, best_k, best_l, best_e);
                    }
                }

                V_(i_mod, j) = v;
            } else {
                V_(i_mod, j) = INF;
            } // end if (i,j form a canonical base pair)

            cand_pos_t ip1_mod = (i + 1) % (MAXLOOP + 1);
            energy_t vi1j = V_(ip1_mod, j);
            energy_t vij1 = V_(i_mod, j - 1);
            energy_t vi1j1 = V_(ip1_mod, j - 1);

            // Checking the dangle positions for W
            energy_t w_v = E_ext_Stem(v, vi1j, vij1, vi1j1, i, j, d);
            // Checking the dangle positions for W
            const energy_t wm_v = E_MLStem(v, vi1j, vij1, vi1j1, i, j, d);

            cand_pos_t k = i;
            cand_pos_t l = j;
            if (params_->model_details.dangles == 1) {
                if (d > 0) {
                    switch (d) {
                    case 1:
                        k = i + 1;
                        break;
                    case 2:
                        l = j - 1;
                        break;
                    case 3:
                        k = i + 1;
                        l = j - 1;
                        break;
                    }
                    if (exists_trace_arrow_from(ta_, k, l) && (wm_v < wm_split || w_v < w_split)) inc_source_ref_count(ta_, k, l);
                }
            }
            energy_t wi_v = INF;
            energy_t wip_v = INF;
            energy_t wi_wmb = INF;
            energy_t wip_wmb = INF;
            energy_t w_wmb = INF, wm_wmb = INF;
            if (pseudoknot) {
                cand_pos_t Bp_ij = tree->Bp(i, j);
                cand_pos_t B_ij = tree->B(i, j);
                cand_pos_t b_ij = tree->b(i, j);
                cand_pos_t bp_ij = tree->bp(i, j);

                // Start of VP ---- Will have to change the bounds to 1 to n instead of 0 to n-1
                bool weakly_closed_ij = tree->weakly_closed(i, j);
                if (weakly_closed_ij || tree->tree[i].pair >= -1 || tree->tree[j].pair >= -1 || ptype_closing == 0) {

                    VP_(i_mod, j) = INF;

                } else {
                    const energy_t vp = compute_VP(i, j, b_ij, bp_ij, Bp_ij, B_ij);

                    VP_(i_mod, j) = vp;
                }

                // -------------------------------------------End of VP----------------------------------------------------------------

                // Start of WMBP
                if ((tree->tree[i].pair >= -1 && tree->tree[i].pair > j) || (tree->tree[j].pair >= -1 && tree->tree[j].pair < i)
                    || (tree->tree[i].pair >= -1 && tree->tree[i].pair < i) || (tree->tree[j].pair >= -1 && j < tree->tree[j].pair)) {
                    WMB_[j] = INF;
                    WMBP_[j] = INF;
                    WMBA_[j] = INF;

                } else {

                    const energy_t wmbp = compute_WMBP(i, j);
                    WMBP_[j] = wmbp;

                    const energy_t wmba = compute_WMBA(j);
                    WMBA_[j] = wmba;

                    const energy_t wmb = compute_WMB(i, j);
                    WMB_[j] = wmb;
                }

                // -------------------------------------------------------End of WMB------------------------------------------------------

                // Start of WI -- the conditions on calculating WI is the same as WIP, so we combine them

                if (!weakly_closed_ij) {
                    WI_[j] = INF;
                    WIP_[j] = INF;
                } else {

                    wi_v = V_(i_mod, j) + PPS_penalty;
                    wip_v = V_(i_mod, j) + bp_penalty;

                    wi_wmb = WMB_[j] + PSM_penalty + PPS_penalty;
                    wip_wmb = WMB_[j] + PSM_penalty + bp_penalty;

                    WI_[j] = std::min({wi_v, wi_wmb, wi_split});
                    WIP_[j] = std::min({wip_v, wip_wmb, wip_split});

                    if ((tree->tree[i - 1].pair > (i - 1) && tree->tree[i - 1].pair > j) || tree->tree[i - 1].pair < (i - 1)) {
                        WI_Bbp[j] = WI_[j];
                        WIP_Bbp[j] = WIP_[j];
                    }
                    if (j + 1 < n_ && tree->tree[j + 1].pair < i) {
                        WIP_Bp[i] = WIP_[j];
                    }
                }

                // ------------------------------------------------End of Wi/Wip--------------------------------------------------

                // start of WV and WVe
                if (!weakly_closed_ij) {
                    cand_pos_t bound_right = std::max(bp_ij, B_ij);
                    // Since bound_left is a min, by doing an uint, I make it so that it's a large number if negative (i.e. not possible)
                    cand_pos_t bound_left = std::min((cand_pos_tu)tree->Bp(i, j), (cand_pos_tu)tree->b(i, j));

                    const energy_t wve = compute_WVe(i, j);
                    const energy_t wv = compute_WV(j, bound_left, bound_right);

                    WVe_[j] = wve;

                    WV_[j] = wv;
                } else {
                    WV_[j] = INF;
                    WVe_[j] = INF;
                }

                // ------------------------------------------------End of WV------------------------------------------------------

                /*
                The order for these should be          i        lp               jp         j           l      ip
                                                     ( (        (    (   (       (          )      )  ) )      ) )
                i and ip are the outer base pair
                lp and l are the closest encompassing base pair to i/ip
                                                                                                //    jp>i   j>jp  ip>i ip>j
                jp and j are some inner base pair;j has the be the closing due to the j=i+4 setup we have
                */
                // Start of BE
                cand_pos_t ip = tree->tree[i].pair; // i's pair ip should be right side so ip = )
                cand_pos_t jp = tree->tree[j].pair; // j's pair jp should be left side so jp = (

                // base case: i.j and ip.jp must be in G

                if (jp > i && j > jp && ip > j && ip > i) { // Don't need to check if they are pairs separately because it is checked this
                    if (tree->tree[jp + 1].pair == j - 1) {
                        // BE_avoided++;
                    } else {
                        energy_t BE = compute_BE(i, j, ip, jp);
                        register_candidate(CLBE_, i, j, BE);
                        register_candidate(CLBEO_, j, i, BE);
                    }

                } else if (i == jp && ip == j) {
                    if (tree->tree[jp + 1].pair == j - 1) { // huh?
                        // BE_avoided++;
                    } else {
                        register_candidate(CLBE_, i, j, 0);
                        register_candidate(CLBEO_, j, i, 0);
                    }
                }

                // // ------------------------------------------------End of BE---------------------------------------------------------
            }

            energy_t vp_min1 = INF;
            energy_t vp_min2 = INF;
            for (auto it = CLVP_[j].begin(); CLVP_[j].end() != it; ++it) {
                const cand_pos_t k = it->first;
                bool can_pair = tree->up[k - 1] >= (k - i);
                energy_t WIk = WI_[k - 1];
                energy_t WIPk = WIP_[k - 1];
                vp_min1 = std::min(vp_min1, WIk + it->second);
                vp_min2 = std::min(vp_min2, WIPk + it->second);
                if (can_pair) vp_min2 = std::min(vp_min2, static_cast<energy_t>((k - i) * cp_penalty) + it->second);
            }
            if ((VP_(i_mod, j) < INFover2)) {
                if ((VP_(i_mod, j) < vp_min1 || VP_(i_mod, j) < vp_min2)) {
                    register_candidate(CLVP_, i, j, VP_(i_mod, j));
                    inc_source_ref_count(taVP_, i, j);
                }
            }

            // Things that needed to happen later like W's wmb
            w_wmb = tree->weakly_closed(i, j) ? WMB_[j] + PS_penalty : INF;
            wm_wmb = tree->weakly_closed(i, j) ? WMB_[j] + PSM_penalty + b_penalty : INF;
            w = std::min({w_v, w_split, w_wmb});
            wm = std::min({wm_v, wm_split, wm_wmb});

            // Some case of WI candidate splits will show INf as a vkj value as we are unable to do dangle versions yet but there are dangle
            // candidates
            if (w_v < w_split || wm_v < wm_split || wi_v < wi_split || wip_v < wip_split || paired) {
                // cand_pos_t k_mod = k%(MAXLOOP+1);
                // Encode the dangles into the energies
                energy_t w_enc = (static_cast<cand_pos_tu>(w_v) << 2) | d;
                energy_t wm_enc = (static_cast<cand_pos_tu>(wm_v) << 2) | d;
                register_candidate(CL_, i, j, V_(i_mod, j), wm_enc, w_enc);
                // always keep arrows starting from candidates
                inc_source_ref_count(ta_, i, j);
            }
            if ((WMB_[j] < INFover2) && (w_wmb < w_split || wm_wmb < wm_split || wi_wmb < wi_split || wip_wmb < wip_split)) {

                register_candidate(CLWMB_, i, j, WMB_[j]);
            }

            W_[j] = w;
            WM_[j] = wm;
            WM2_[j] = std::min(wm2_split, WMB_[j] + PSM_penalty + b_penalty);

        } // end loop j
        rotate_arrays();

        // Clean up trace arrows in i+MAXLOOP+1
        if (garbage_collect_ && i + MAXLOOP + 1 <= n_) {
            gc_row(ta_, i + MAXLOOP + 1);
            gc_row(taVP_, i + MAXLOOP + 1);
        }
        // Reallocate candidate lists in i
        for (auto &x : CL_) {
            if (x.capacity() > 1.5 * x.size()) {
                cand_list_td1 vec(x.size());
                copy(x.begin(), x.end(), vec.begin());
                vec.swap(x);
            }
        }
        for (auto &x : CLVP_) {
            if (x.capacity() > 1.5 * x.size()) {
                cand_list_t vec(x.size());
                copy(x.begin(), x.end(), vec.begin());
                vec.swap(x);
            }
        }
        for (auto &x : CLWMB_) {
            if (x.capacity() > 1.5 * x.size()) {
                cand_list_t vec(x.size());
                copy(x.begin(), x.end(), vec.begin());
                vec.swap(x);
            }
        }
        for (auto &x : CLBE_) {
            if (x.capacity() > 1.5 * x.size()) {
                cand_list_t vec(x.size());
                copy(x.begin(), x.end(), vec.begin());
                vec.swap(x);
            }
        }

        compactify(ta_);
        compactify(taVP_);
    }
    return W_[n_];
}