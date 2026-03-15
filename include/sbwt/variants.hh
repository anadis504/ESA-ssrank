#pragma once

#include <cstring>
#include <filesystem>
#include <string>

#include "MEF.hpp"
#include "SBWT.hh"
#include "SubsetBlockedSplitCorrectionSets.hh"
#include "SubsetConcatRank.hh"
#include "SubsetConcatSplitLenthsRank.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetNewConcatRank.hh"
#include "SubsetNewPlainConcatRank.hh"
#include "SubsetNewSplitRank.hh"
#include "SubsetSplitCorrectionSets.hh"
#include "SubsetSplitRank.hh"
#include "SubsetSplitSmallerSizeRank.hh"
#include "SubsetWT.hh"
#include "cxxopts.hpp"
#include "globals.hh"
#include "stdlib_printing.hh"
#include "SubsetBlockedSplit.hh"
#include "SubsetSplitRankPred8.hh"

namespace sbwt {

// matrices
typedef SBWT<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    plain_matrix_sbwt_t;
typedef SBWT<
    SubsetMatrixRank<sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type>>
    rrr_matrix_sbwt_t;
typedef SBWT<SubsetMatrixRank<mod_ef_vector<>, mod_ef_vector<>::rank_1_type>>
    mef_matrix_sbwt_t;  // Currently does not support extracting all k-mers
                        // because mod_ef_vector does not support access.

// splits
typedef SBWT<SubsetSplitRank<sdsl::bit_vector, sdsl::rank_support_v5<>,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>>
    plain_split_sbwt_t;

typedef SBWT<
    SubsetSplitRank<sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type,
                    sdsl::bit_vector, sdsl::rank_support_v5<>>>
    rrr_split_sbwt_t;

typedef SBWT<SubsetSplitRank<mod_ef_vector<>, mod_ef_vector<>::rank_1_type,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>>
    mef_split_sbwt_t;  // Currently does not support extracting all k-mers
                       // because mod_ef_vector does not support access.

// concats
typedef SBWT<SubsetConcatRank<
    sdsl::bit_vector, sdsl::bit_vector::select_0_type,
    sdsl::wt_blcd<sdsl::bit_vector, sdsl::rank_support_v5<>,
                  sdsl::select_support_scan<1>, sdsl::select_support_scan<0>>>>
    plain_concat_sbwt_t;

typedef SBWT<
    SubsetConcatRank<sd_vector<>, sd_vector<>::select_0_type,
                     sdsl::wt_blcd<rrr_vector<63>, rrr_vector<>::rank_1_type,
                                   rrr_vector<>::select_1_type,
                                   rrr_vector<>::select_0_type>>>
    mef_concat_sbwt_t;  // Currently does not support extracting all k-mers
                        // because mod_ef_vector does not support access.

// wavelet trees
typedef SBWT<SubsetWT<
    sdsl::wt_blcd<sdsl::bit_vector, sdsl::rank_support_v5<>,
                  sdsl::select_support_scan<1>, sdsl::select_support_scan<0>>>>
    plain_sswt_sbwt_t;

typedef SBWT<SubsetWT<
    sdsl::wt_blcd<sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type,
                  rrr_vector<>::select_1_type, rrr_vector<>::select_0_type>>>
    rrr_sswt_sbwt_t;

typedef SBWT<SubsetNewConcatRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    new_concat_sbwt_t;
typedef SBWT<
    SubsetSplitCorrectionSetsRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    split_correction_sets_sbwt_t;
typedef SBWT<
    SubsetSplitSmallerSizeRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    split_smaller_size_sbwt_t;
typedef SBWT<
    SubsetConcatSplitLengthsRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    concat_split_lengths_sbwt_t;
typedef SBWT<
    SubsetNewPlainConcatRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    new_plain_concat_sbwt_t;
typedef SBWT<SubsetNewSplitRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    new_split_sbwt_t;
typedef SBWT<
    SubsetBlockedSplitCorrectionSetsRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    blocked_split_correction_sets_sbwt_t;
typedef SBWT<
    SubsetBlockedSplitRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    blocked_split_sbwt_t;
typedef SBWT<SubsetSplitRankPred8<mod_ef_vector<>, mod_ef_vector<>::rank_1_type,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>>
    pred8_split_sbwt_t;

}  // namespace sbwt
