#ifndef TNPTREE_H
#define TNPTREE_H

#include "../git/foliage/include/foliage.h"

#include "TTree.h"

#include <vector>

#define B_VAL_TNP(ACTION, ...)                                              \
    ACTION(float,           tag_pt,                     ## __VA_ARGS__)     \
    ACTION(float,           tag_eta,                    ## __VA_ARGS__)     \
    ACTION(float,           tag_phi,                    ## __VA_ARGS__)     \
    ACTION(float,           probe_pt,                   ## __VA_ARGS__)     \
    ACTION(float,           probe_eta,                  ## __VA_ARGS__)     \
    ACTION(float,           probe_phi,                  ## __VA_ARGS__)     \
    ACTION(float,           mass,                       ## __VA_ARGS__)     \
    ACTION(float,           dr2_l1,                     ## __VA_ARGS__)     \
    ACTION(float,           dr2_hlt,                    ## __VA_ARGS__)     \
    ACTION(int32_t,         pass_l1,                    ## __VA_ARGS__)     \
    ACTION(int32_t,         pass_hlt,                   ## __VA_ARGS__)     \
    ACTION(int32_t,         pass_veto_id,               ## __VA_ARGS__)     \
    ACTION(int32_t,         pass_loose_id,              ## __VA_ARGS__)     \
    ACTION(int32_t,         pass_medium_id,             ## __VA_ARGS__)     \
    ACTION(int32_t,         pass_tight_id,              ## __VA_ARGS__)     \
    ACTION(float,           weight,                     ## __VA_ARGS__)     \
    ACTION(int32_t,         pass_v1,                    ## __VA_ARGS__)     \
    ACTION(int32_t,         pass_v2,                    ## __VA_ARGS__)     \

class tnptree {
  public:
    tnptree(TTree* t) {
        B_VAL_TNP(SETMONE)

        B_VAL_TNP(BRANCHVAL, t)
    }

    ~tnptree() = default;

    B_VAL_TNP(DECLVAL)
};

#endif /* TNPTREE_H */
