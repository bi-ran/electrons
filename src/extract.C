#include "TAxis.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

#include <string>
#include <vector>

#include "../include/etree.h"
#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/foliage/include/foliage.h"
#include "../git/foliage/include/eggen.h"
#include "../git/foliage/include/electrons.h"
#include "../git/foliage/include/triggers.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"
#include "../git/tricks-and-treats/include/train.h"

auto oadphi = [](float phi1, float phi2) {
    return revert_radian(convert_radian(phi1) - convert_radian(phi2));
};

int extract(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto paths = conf->get<std::vector<std::string>>("paths");
    auto skim = conf->get<std::vector<std::string>>("skim");
    auto weights = conf->get<std::string>("weights");

    auto heavyion = conf->get<bool>("heavyion");
    auto mc_branches = conf->get<bool>("mc_branches");
    auto hlt_branches = conf->get<bool>("hlt_branches");

    auto apply_weight = !weights.empty();

    TH2F* hweights = nullptr;
    TAxis* xaxis = nullptr;
    TAxis* yaxis = nullptr;
    if (apply_weight) {
        TFile* fw = new TFile(weights.data(), "read");
        hweights = (TH2F*)fw->Get("hweights");
        xaxis = hweights->GetXaxis();
        yaxis = hweights->GetYaxis();
    }

    auto forest = new train(files);
    auto chain_eg = forest->attach("ggHiNtuplizerGED/EventTree", true);
    auto chain_hlt = forest->attach("hltanalysis/HltTree", hlt_branches);
    auto chain_evt = forest->attach("hiEvtAnalyzer/HiTree", true);

    (*forest)();

    auto tegg = harvest<eggen>(chain_eg, mc_branches);
    auto tegm = harvest<electrons>(chain_eg);
    auto thlt = harvest<triggers>(chain_hlt, paths);
    auto tevt = harvest<event>(chain_evt, mc_branches);

    TTree::SetMaxTreeSize(1000000000000LL);

    TFile* fout = new TFile(output, "recreate");
    TTree* tout = new TTree("e", "electrons");
    auto te = new etree(tout, mc_branches, hlt_branches);

    int64_t nentries = forest->count();
    if (max_entries) nentries = std::min(nentries, max_entries);
    for (int64_t i = 0; i < nentries; ++i) {
        te->clear();

        if (i % 10000 == 0)
            printf("entry: %li/%li\n", i, nentries);

        forest->get(i);

        if (tegm->nEle < 1) { continue; }

        if (!skim.empty()) {
            bool pass_skim = false;
            for (auto const& path : skim)
                if (thlt->accept(path) == 1)
                    pass_skim = true;

            if (!pass_skim) { continue; }
        }

        te->copy(tegg);
        te->copy(tegm);
        te->copy(thlt);
        te->copy(tevt);

        if (!heavyion) {
            te->hiBin = 0;
            te->hiHF = 0;
            te->Ncoll = 1;
        }

        if (mc_branches) {
            constexpr float max_dr2 = 0.15 * 0.15;

            for (int32_t j = 0; j < te->nEle; ++j) {
                auto weight = !apply_weight ? 1. : hweights->GetBinContent(
                    xaxis->FindBin((*te->elePt)[j]),
                    yaxis->FindBin((*te->eleEta)[j]));
                te->ele_weight->push_back(weight);

                float maxpt = -1;
                int32_t match = -1;

                float ele_eta = (*te->eleEta)[j];
                float ele_phi = (*te->elePhi)[j];

                for (int32_t k = 0; k < te->nMC; ++k) {
                    if (std::abs((*te->mcPID)[k]) != 11) { continue; }
                    if ((*te->mcStatus)[k] != 1) { continue; }
                    if ((*te->mcPt)[k] < maxpt) { continue; }

                    float deta = ele_eta - (*te->mcEta)[k];
                    float dphi = oadphi(ele_phi, (*te->mcPhi)[k]);
                    float dr2 = dphi * dphi + deta * deta;

                    if (dr2 < max_dr2) {
                        maxpt = (*te->mcPt)[k];
                        match = k;
                    }
                }

                te->gen_index->push_back(match);
            }
        }

        tout->Fill();
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();
    
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return extract(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
