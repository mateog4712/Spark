#include "part_func.hh"
#include "h_externs.hh"
#include "pf_globals.hh"
#include "ViennaRNA/utils.hh"


Spark_PF::Spark_PF(const std::string &seq, std::string restricted, sparse_tree *tree, bool pseudoknot, bool pk_only, bool garbage_collect,bool mark_candidates) : seq_(seq), n_(seq.length()), exp_params_(vrna_exp_params(NULL)), garbage_collect_(garbage_collect), mark_candidates_(mark_candidates), ta_(n_), taVP_(n_) {
    make_pair_matrix();

    S_ = encode_sequence(seq.c_str(), 0);
    S1_ = encode_sequence(seq.c_str(), 1);
    this->tree = tree;
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

void Spark_PF::rescale_pk_globals() {
    double kT = exp_params_->model_details.betaScale * (exp_params_->model_details.temperature + K0) * GASCONST; /* kT in cal/mol  */
    double TT = (exp_params_->model_details.temperature + K0) / (Tmeasure);
    int pf_smooth = exp_params_->model_details.pf_smooth;
    // ShapeData->rescale_calculate(kT,TT,pf_smooth);
    // for(int i = 0; i<30;++i){
    //     std::cout << i << " " << ShapeData->get_calculated(i) << " " << ShapeData->get_expcalculated(i) << std::endl;
    // }

    expPS_penalty = RESCALE_BF(PS_penalty, PS_penalty * 3, TT, kT);
    expPSM_penalty = RESCALE_BF(PSM_penalty, PSM_penalty * 3, TT, kT);
    expPSP_penalty = RESCALE_BF(PSP_penalty, PSP_penalty * 3, TT, kT);
    expPB_penalty = RESCALE_BF(PB_penalty, PB_penalty * 3, TT, kT);
    expPUP_penalty = RESCALE_BF(PUP_penalty, PUP_penalty * 3, TT, kT);
    expPPS_penalty = RESCALE_BF(PPS_penalty, PPS_penalty * 3, TT, kT);

    expa_penalty = RESCALE_BF(a_penalty, ML_closingdH, TT, kT);
    expb_penalty = RESCALE_BF(b_penalty, ML_interndH, TT, kT);
    expc_penalty = RESCALE_BF(c_penalty, ML_BASEdH, TT, kT);

    expap_penalty = RESCALE_BF(ap_penalty, ap_penalty * 3, TT, kT);
    expbp_penalty = RESCALE_BF(bp_penalty, bp_penalty * 3, TT, kT);
    expcp_penalty = RESCALE_BF(cp_penalty, cp_penalty * 3, TT, kT);
}

pf_t Spark_PF::exp_Extloop(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[i]][S_[j]];

    if (exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i - 1] : -1;
        base_type sj1 = j < n_ ? S_[j + 1] : -1;
        return exp_E_ExtLoop(tt, si1, sj1, exp_params_);
    } else {
        return exp_E_ExtLoop(tt, -1, -1, exp_params_);
    }
}

pf_t Spark_PF::exp_MLstem(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[i]][S_[j]];
    if (exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i - 1] : -1;
        base_type sj1 = j < n_ ? S_[j + 1] : -1;
        return exp_E_MLstem(tt, si1, sj1, exp_params_);
    } else {
        return exp_E_MLstem(tt, -1, -1, exp_params_);
    }
}

pf_t Spark_PF::exp_Mbloop(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[j]][S_[i]];
    if (exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i + 1] : -1;
        base_type sj1 = j < n_ ? S_[j - 1] : -1;
        return exp_E_MLstem(tt, sj1, si1, exp_params_);
    } else {
        return exp_E_MLstem(tt, -1, -1, exp_params_);
    }
}

pf_t Spark_PF::HairpinE(cand_pos_t i, cand_pos_t j) {
    const int ptype_closing = pair[S_[i]][S_[j]];
    if (ptype_closing == 0) return 0;
    pf_t e_h = static_cast<pf_t>(exp_E_Hairpin(j - i - 1, ptype_closing, S1_[i + 1], S1_[j - 1], &seq_.c_str()[i - 1], exp_params_));
    e_h *= scale[j - i + 1];
    return e_h;
}