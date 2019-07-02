#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include <string>
#include <vector>

#include "../include/etree.h"

#include "../git/config/include/configurer.h"

#include "../git/foliage/include/foliage.h"
#include "../git/foliage/include/electrons.h"
#include "../git/foliage/include/triggers.h"

#include "../git/tricks-and-treats/include/train.h"
#include "../git/tricks-and-treats/include/maglev.h"

int extract(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto mc_branches = conf->get<bool>("mc_branches");
    auto hlt_branches = conf->get<bool>("hlt_branches");
    auto paths = conf->get<std::vector<std::string>>("paths");
    auto skim = conf->get<std::vector<std::string>>("skim");

    auto forest = new train(files);
    auto chain_eg = forest->attach("ggHiNtuplizerGED/EventTree", true);
    auto chain_hlt = forest->attach("hltanalysisReco/HltTree", hlt_branches);

    (*forest)();

    auto tree_eg = new electrons(chain_eg, mc_branches);
    auto tree_hlt = new triggers(chain_hlt, paths);

    TTree::SetMaxTreeSize(1000000000000LL);

    TFile* fout = new TFile(output, "recreate");
    TTree* tout = new TTree("e", "electrons");
    auto tree_e = new etree(tout, mc_branches, hlt_branches);

    int64_t nentries = forest->count();
    if (max_entries) nentries = std::min(nentries, max_entries);
    for (int64_t i = 0; i < nentries; ++i) {
        tree_e->clear();

        if (i % 10000 == 0)
            printf("entry: %li/%li\n", i, nentries);

        forest->get(i);

        if (tree_eg->nEle < 1) { continue; }

        if (!skim.empty()) {
            bool pass_skim = false;
            for (auto const& path : skim)
                if (tree_hlt->accept(path) == 1)
                    pass_skim = true;

            if (!pass_skim) { continue; }
        }

        tree_e->copy(tree_eg);
        tree_e->copy(tree_hlt);

        /* extra variables */
        if (mc_branches) {
            constexpr float max_dr2 = 0.15 * 0.15;
            for (int32_t j = 0; j < tree_e->nEle; ++j) {
                float maxpt = -1;
                int32_t match = -1;

                float ele_eta = (*tree_e->eleEta)[j];
                float ele_phi = (*tree_e->elePhi)[j];

                for (int32_t k = 0; k < tree_e->nMC; ++k) {
                    if (std::abs((*tree_e->mcPID)[k]) != 11) { continue; }
                    if ((*tree_e->mcStatus)[k] != 1) { continue; }
                    if ((*tree_e->mcPt)[k] < maxpt) { continue; }

                    float deta = ele_eta - (*tree_e->mcEta)[k];
                    float dphi = ml_dphi(ele_phi, (*tree_e->mcPhi)[k]);
                    float dr2 = dphi * dphi + deta * deta;

                    if (dr2 < max_dr2) {
                        maxpt = (*tree_e->mcPt)[k];
                        match = k;
                    }
                }

                tree_e->gen_index->push_back(match);
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

    return 0;
}
