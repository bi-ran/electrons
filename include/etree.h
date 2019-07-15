#ifndef ETREE_H
#define ETREE_H

#include "../git/foliage/include/foliage.h"
#include "../git/foliage/include/electrons.h"
#include "../git/foliage/include/triggers.h"
#include "../git/foliage/include/event.h"

#include "TTree.h"

#include <vector>

#define B_VEC_TRG(ACTION, ...)                                              \
    ACTION(sv<int32_t>,     accepts,                    ## __VA_ARGS__)     \

#define B_VEC_ELE_EXT(ACTION, ...)                                          \
    ACTION(sv<int32_t>,     gen_index,                  ## __VA_ARGS__)     \

class etree {
  public:
    etree(TTree* t, bool gen, bool hlt)
            : _gen(gen),
              _hlt(hlt) {
        B_VAL_ELE_RECO(SETMONE)
        B_VEC_ELE_RECO(ALLOCOBJ)

        B_VAL_EVT_RECO(SETMONE)

        if (_gen) {
            B_VAL_ELE_GEN(SETMONE)
            B_VEC_ELE_GEN(ALLOCOBJ)
            B_VEC_ELE_EXT(ALLOCOBJ)

            B_VAL_EVT_GEN(SETMONE)
        }

        if (_hlt) {
            B_VEC_TRG(ALLOCOBJ)
        }

        branch(t);
    }

    etree(TTree* t, bool gen)
        : etree(t, gen, false) { }

    etree(bool gen, bool hlt, TTree* t)
            : _gen(gen),
              _hlt(hlt) {
        B_VAL_ELE_RECO(SETZERO)
        B_VEC_ELE_RECO(SETZERO)

        B_VAL_EVT_RECO(SETZERO)

        if (_gen) {
            B_VAL_ELE_GEN(SETZERO)
            B_VEC_ELE_GEN(SETZERO)
            B_VEC_ELE_EXT(SETZERO)

            B_VAL_EVT_GEN(SETZERO)
        }

        if (_hlt) {
            B_VEC_TRG(SETZERO)
        }

        read(t);
    }

    etree(bool gen, TTree* t)
        : etree(gen, false, t) { }

    ~etree() = default;

    void clear() {
        B_VEC_ELE_RECO(CLEAROBJ)

        if (_gen) {
            B_VEC_ELE_GEN(CLEAROBJ)
            B_VEC_ELE_EXT(CLEAROBJ)
        }

        if (_hlt) {
            B_VEC_TRG(CLEAROBJ)
        }
    }

    void copy(electrons* t) {
        B_VAL_ELE_RECO(COPYVAL, t)
        B_VEC_ELE_RECO(COPYOBJ, t)

        if (_gen) {
            B_VAL_ELE_GEN(COPYVAL, t)
            B_VEC_ELE_GEN(COPYOBJ, t)
        }
    }

    void copy(event* t) {
        B_VAL_EVT_RECO(COPYVAL, t)

        if (_gen) {
            B_VAL_EVT_GEN(COPYVAL, t)
        }
    }

    void copy(triggers* t) {
        if (_hlt) {
            B_VEC_TRG(COPYPTR, t, t->size())
        }
    }

    B_VAL_ELE_RECO(DECLVAL)
    B_VEC_ELE_RECO(DECLPTR)
    B_VAL_ELE_GEN(DECLVAL)
    B_VEC_ELE_GEN(DECLPTR)
    B_VEC_ELE_EXT(DECLPTR)
    B_VAL_EVT_RECO(DECLVAL)
    B_VAL_EVT_GEN(DECLVAL)
    B_VEC_TRG(DECLPTR)

  private:
    void branch(TTree* t) {
        B_VAL_ELE_RECO(BRANCHVAL, t)
        B_VEC_ELE_RECO(BRANCHPTR, t)

        B_VAL_EVT_RECO(BRANCHVAL, t)

        if (_gen) {
            B_VAL_ELE_GEN(BRANCHVAL, t)
            B_VEC_ELE_GEN(BRANCHPTR, t)
            B_VEC_ELE_EXT(BRANCHPTR, t)

            B_VAL_EVT_GEN(BRANCHVAL, t)
        }

        if (_hlt) {
            B_VEC_TRG(BRANCHPTR, t)
        }
    }

    void read(TTree* t) {
        B_VAL_ELE_RECO(SETVALADDR, t)
        B_VEC_ELE_RECO(SETVALADDR, t)

        B_VAL_EVT_RECO(SETVALADDR, t)

        if (_gen) {
            B_VAL_ELE_GEN(SETVALADDR, t)
            B_VEC_ELE_GEN(SETVALADDR, t)
            B_VEC_ELE_EXT(SETVALADDR, t)

            B_VAL_EVT_GEN(SETVALADDR, t)
        }

        if (_hlt) {
            B_VEC_TRG(SETVALADDR, t)
        }
    }

    bool _gen;
    bool _hlt;
};

#endif /* ETREE_H */
