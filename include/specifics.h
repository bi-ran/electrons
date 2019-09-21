#ifndef SPECIFICS_H
#define SPECIFICS_H

template <typename T>
bool within_hem_failure_region(T* t, int64_t i) {
    return ((*t->eleSCEta)[i] < -1.3
        && (*t->eleSCPhi)[i] < -0.87
        && (*t->eleSCPhi)[i] > -1.57);
}

template <typename T>
bool passes_basic_electron_selections(T* t, int64_t i) {
    return (*t->eleConvVeto)[i] && (*t->eleMissHits)[i] <= 1
        && (*t->eleIP3D)[i] < 0.03;
}

enum ip { incl, cent, peri, nip };
enum det { barrel, endcap, ndet };
enum wp { veto, loose, medium, tight, claustro, nwp };
enum var { hoe, see, deta, dphi, eop, nele,
        /* hoe, see, */ iso = 2, npho};

constexpr float ecuts[ip::nip][det::ndet][wp::nwp][var::nele] = {
    {    /* ip::incl */
        {    /* det::barrel */
            { 0.0607, 0.0103, 0.0048, 0.0625, 0.1327 }, /* wp::veto */
            { 0.0271, 0.0102, 0.0032, 0.0394, 0.0530 }, /* wp::loose */
            { 0.0246, 0.0097, 0.0024, 0.0292, 0.0447 }, /* wp::medium */
            { 0.0205, 0.0093, 0.0023, 0.0279, 0.0392 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }, { /* det::endcap */
            { 0.0452, 0.0306, 0.0066, 0.0881, 0.9066 }, /* wp::veto */
            { 0.0375, 0.0295, 0.0057, 0.0382, 0.0236 }, /* wp::loose */
            { 0.0113, 0.0294, 0.0056, 0.0283, 0.0234 }, /* wp::medium */
            { 0.0014, 0.0283, 0.0047, 0.0267, 0.0154 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }
    }, { /* ip::cent */
        {    /* det::barrel */
            { 0.2733, 0.0147, 0.0041, 0.0853, 0.0367 }, /* wp::veto */
            { 0.1616, 0.0135, 0.0038, 0.0376, 0.0177 }, /* wp::loose */
            { 0.1589, 0.0116, 0.0037, 0.0224, 0.0173 }, /* wp::medium */
            { 0.1459, 0.0104, 0.0029, 0.0206, 0.0105 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }, { /* det::endcap */
            { 0.1898, 0.0480, 0.0097, 0.2348, 0.0300 }, /* wp::veto */
            { 0.1317, 0.0466, 0.0063, 0.1186, 0.0201 }, /* wp::loose */
            { 0.1092, 0.0418, 0.0062, 0.0373, 0.0133 }, /* wp::medium */
            { 0.0925, 0.0358, 0.0051, 0.0266, 0.0065 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }
    }, { /* ip::peri */
        {    /* det::barrel */
            { 0.1814, 0.0113, 0.0037, 0.1280, 0.1065 }, /* wp::veto */
            { 0.1268, 0.0107, 0.0035, 0.0327, 0.0774 }, /* wp::loose */
            { 0.0311, 0.0101, 0.0033, 0.0210, 0.0701 }, /* wp::medium */
            { 0.0067, 0.0099, 0.0026, 0.0170, 0.0077 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }, { /* det::endcap */
            { 0.1138, 0.0376, 0.0074, 0.2085, 0.0237 }, /* wp::veto */
            { 0.0977, 0.0339, 0.0067, 0.0838, 0.0193 }, /* wp::loose */
            { 0.0810, 0.0316, 0.0051, 0.0384, 0.0192 }, /* wp::medium */
            { 0.0655, 0.0288, 0.0044, 0.0266, 0.0123 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }
    }
};

constexpr float acuts[det::ndet][2] = {
    { -1, 1.4442 }, /* det::barrel */
    { 1.556, 2.1 }  /* det::endcap */
};

template <det T, typename U>
bool within_acceptance(U* t, int64_t i) {
    auto abs_eta = std::abs((*t->eleSCEta)[i]);
    return abs_eta > acuts[T][0] && abs_eta < acuts[T][1];
}

template <ip T, det U, wp V, typename W>
bool passes_electron_id(W* t, int64_t i) {
    if (!passes_basic_electron_selections(t, i)) { return false; }
    if (!within_acceptance<U>(t, i)) { return false; }

    return (*t->eleHoverEBc)[i] < ecuts[T][U][V][var::hoe]
        && (*t->eleSigmaIEtaIEta_2012)[i] < ecuts[T][U][V][var::see]
        && std::abs((*t->eledEtaSeedAtVtx)[i]) < ecuts[T][U][V][var::deta]
        && std::abs((*t->eledPhiAtVtx)[i]) < ecuts[T][U][V][var::dphi]
        && std::abs((*t->eleEoverPInv)[i]) < ecuts[T][U][V][var::eop];
}

template <ip T, wp U, typename V>
bool passes_electron_id(V* t, int64_t i) {
    if (!passes_basic_electron_selections(t, i)) { return false; }
    auto d = within_acceptance<det::barrel>(t, i) ? det::barrel
        : within_acceptance<det::endcap>(t, i) ? det::endcap : det::ndet;
    if (d == det::ndet) { return false; }

    return (*t->eleHoverEBc)[i] < ecuts[T][d][U][var::hoe]
        && (*t->eleSigmaIEtaIEta_2012)[i] < ecuts[T][d][U][var::see]
        && std::abs((*t->eledEtaSeedAtVtx)[i]) < ecuts[T][d][U][var::deta]
        && std::abs((*t->eledPhiAtVtx)[i]) < ecuts[T][d][U][var::dphi]
        && std::abs((*t->eleEoverPInv)[i]) < ecuts[T][d][U][var::eop];
}

template <det T, wp U, typename V>
bool passes_electron_id(V* t, int64_t i, bool heavyion) {
    if (!passes_basic_electron_selections(t, i)) { return false; }
    if (!within_acceptance<T>(t, i)) { return false; }

    auto iptype = heavyion ? (t->hiBin < 60 ? ip::cent : ip::peri) : ip::incl;
    return (*t->eleHoverEBc)[i] < ecuts[iptype][T][U][var::hoe]
        && (*t->eleSigmaIEtaIEta_2012)[i] < ecuts[iptype][T][U][var::see]
        && std::abs((*t->eledEtaSeedAtVtx)[i]) < ecuts[iptype][T][U][var::deta]
        && std::abs((*t->eledPhiAtVtx)[i]) < ecuts[iptype][T][U][var::dphi]
        && std::abs((*t->eleEoverPInv)[i]) < ecuts[iptype][T][U][var::eop];
}

#endif /* SPECIFICS_H */
