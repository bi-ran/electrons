#include "../include/lambdas.h"
#include "../include/specifics.h"
#include "../include/tnptree.h"

#include "../git/config/include/configurer.h"

#include "../git/foliage/include/electrons.h"
#include "../git/foliage/include/hltobjs.h"
#include "../git/foliage/include/l1objs.h"
#include "../git/foliage/include/triggers.h"

#include "../git/tricks-and-treats/include/maglev.h"
#include "../git/tricks-and-treats/include/overflow_angles.h"
#include "../git/tricks-and-treats/include/train.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

using namespace std::literals::string_literals;

auto oadphi = [](float phi1, float phi2) {
    return revert_radian(convert_radian(phi1) - convert_radian(phi2));
};

std::vector<int64_t> clean(std::vector<double> const& values, int64_t size) {
    std::vector<double> handled;
    std::vector<int64_t> indices;

    for (int64_t i = 0; i < static_cast<int64_t>(values.size()); ++i) {
        if (std::find(std::begin(handled), std::end(handled), values[i])
            != std::end(handled)) { continue; }
        if (std::count(std::begin(values), std::end(values), values[i])
            == size) { handled.push_back(values[i]); indices.push_back(i); }
    }

    return indices;
}

float nearest_neighbour(float pt, float eta, float phi,
                        std::vector<float> const& other_pt,
                        std::vector<float> const& other_eta,
                        std::vector<float> const& other_phi) {
    float mindr2 = 999.f;

    auto count = static_cast<int64_t>(other_eta.size());
    for (int64_t i = 0; i < count; ++i) {
        if (other_pt[i] < pt) { continue; }

        float deta = eta - other_eta[i];
        float dphi = oadphi(phi, other_phi[i]);
        float dr2 = deta * deta + dphi * dphi;

        if (dr2 < mindr2) { mindr2 = dr2; }
    }

    return mindr2;
}

int flatten(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto paths = conf->get<std::vector<std::string>>("paths");
    auto tree = conf->get<std::string>("tree");

    auto tag_pt_min = conf->get<float>("tag_pt_min");

    auto l1pt = conf->get<float>("l1pt");
    auto l1dr = conf->get<float>("l1dr");
    auto hltpath = conf->get<std::string>("hltpath");
    auto hltsteps = conf->get<uint32_t>("hltsteps");
    auto hltpt = conf->get<float>("hltpt");
    auto hltdr = conf->get<float>("hltdr");

    auto l1dr2 = l1dr * l1dr;
    auto hltdr2 = hltdr * hltdr;

    auto forest = new train(files);
    auto chain_eg = forest->attach("ggHiNtuplizerGED/EventTree", true);
    auto chain_l1 = forest->attach("l1object/L1UpgradeFlatTree", true);
    auto chain_hlt = forest->attach(("hltobject/"s + hltpath).data(), true);
    auto chain_trg = forest->attach((tree + "/HltTree").data(), true);

    (*forest)();

    auto tree_egm = harvest<electrons>(chain_eg);
    auto tree_trg = harvest<triggers>(chain_trg, paths);
    auto tree_l1 = harvest<l1objs>(chain_l1);
    auto tree_hlt = harvest<hltobjs>(chain_hlt);

    TTree::SetMaxTreeSize(1000000000000LL);

    TFile* fout = new TFile(output, "recreate");
    TTree* tout = new TTree("tnp", "electrons");
    auto tree_tnp = new tnptree(tout);

    int64_t nentries = forest->count();
    if (max_entries) nentries = std::min(nentries, max_entries);
    for (int64_t i = 0; i < nentries; ++i) {
        if (i % 10000 == 0)
            printf("entry: %li/%li\n", i, nentries);

        tree_trg->reset();
        forest->get(i);

        if (tree_trg->accept() != 1)
            continue;

        if (tree_egm->nEle < 2) { continue; }

        std::vector<float> l1mindr2;
        std::vector<float> hltmindr2;
        std::vector<int32_t> veto_id;
        std::vector<int32_t> loose_id;
        std::vector<int32_t> medium_id;
        std::vector<int32_t> tight_id;

        int64_t tag = -1;

        for (int64_t j = 0; j < tree_egm->nEle; ++j) {
            float eta = (*tree_egm->eleEta)[j];
            float phi = (*tree_egm->elePhi)[j];

            l1mindr2.push_back(nearest_neighbour(l1pt, eta, phi,
                *tree_l1->egEt, *tree_l1->egEta, *tree_l1->egPhi));

            /* select hlt objects passing final filter */
            auto indices = clean(*tree_hlt->pt, hltsteps);

            std::vector<float> ptfinal;
            std::vector<float> etafinal;
            std::vector<float> phifinal;

            for (auto index : indices) {
                ptfinal.push_back((*tree_hlt->pt)[index]);
                etafinal.push_back((*tree_hlt->eta)[index]);
                phifinal.push_back((*tree_hlt->phi)[index]);
            }

            hltmindr2.push_back(nearest_neighbour(hltpt, eta, phi,
                ptfinal, etafinal, phifinal));

            /* evaluate id */
            veto_id.push_back(passes_electron_id<
                ip::incl, wp::veto, electrons>(tree_egm, j));
            loose_id.push_back(passes_electron_id<
                ip::incl, wp::loose, electrons>(tree_egm, j));
            medium_id.push_back(passes_electron_id<
                ip::incl, wp::medium, electrons>(tree_egm, j));
            tight_id.push_back(passes_electron_id<
                ip::incl, wp::tight, electrons>(tree_egm, j));

            if (tag < 0 && l1mindr2.back() < l1dr2 && hltmindr2.back() < hltdr2
                    && tight_id.back() && (*tree_egm->elePt)[j] > tag_pt_min) {
                tag = j; }
        }

        if (tag < 0) { continue; }

        for (int64_t j = 0; j < tree_egm->nEle; ++j) {
            if (j == tag) { continue; }

            if ((*tree_egm->eleCharge)[tag] == (*tree_egm->eleCharge)[j])
                continue;

            tree_tnp->tag_pt = (*tree_egm->elePt)[tag];
            tree_tnp->tag_eta = (*tree_egm->eleEta)[tag];
            tree_tnp->tag_phi = (*tree_egm->elePhi)[tag];
            tree_tnp->probe_pt = (*tree_egm->elePt)[j];
            tree_tnp->probe_eta = (*tree_egm->eleEta)[j];
            tree_tnp->probe_abseta = std::abs(tree_tnp->probe_eta);
            tree_tnp->probe_phi = (*tree_egm->elePhi)[j];
            tree_tnp->dr2_l1 = l1mindr2[j];
            tree_tnp->dr2_hlt = hltmindr2[j];
            tree_tnp->pass_l1 = l1mindr2[j] < l1dr2;
            tree_tnp->pass_hlt = hltmindr2[j] < hltdr2;
            tree_tnp->pass_veto_id = veto_id[j];
            tree_tnp->pass_loose_id = loose_id[j];
            tree_tnp->pass_medium_id = medium_id[j];
            tree_tnp->pass_tight_id = tight_id[j];

            tree_tnp->mass = std::sqrt(
                ml_invariant_mass<coords::collider>(
                    (*tree_egm->elePt)[tag],
                    (*tree_egm->eleEta)[tag],
                    (*tree_egm->elePhi)[tag],
                    0.000511f,
                    (*tree_egm->elePt)[j],
                    (*tree_egm->eleEta)[j],
                    (*tree_egm->elePhi)[j],
                    0.000511f));

            tree_tnp->weight = 1.f;

            /* special variables for trigger versions */
            tree_tnp->pass_v1 = tree_tnp->pass_hlt && tree_trg->accept(
                "HLT_HIEle20_WPLoose_Gsf_v1"s) == 1;
            tree_tnp->pass_v2 = tree_tnp->pass_hlt && tree_trg->accept(
                "HLT_HIEle20_WPLoose_Gsf_v2"s) == 1;

            tout->Fill();
        }
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return flatten(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
