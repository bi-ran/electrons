#include "../include/etree.h"
#include "../include/lambdas.h"
#include "../include/specifics.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/maglev.h"
#include "../git/tricks-and-treats/include/trunk.h"

#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

double f_double_sided_crystal_ball(double* x, double* params) {
    double x0 = x[0];
    /* gaussian */
    double n = params[0];
    double mean = params[1];
    double sigma = params[2];
    /* transition */
    double a1 = params[3];
    double n1 = params[4];
    double a2 = params[5];
    double n2 = params[6];

    double u = (x0 - mean) / sigma;
    double A1 = TMath::Power(n1 / TMath::Abs(a1), n1) * TMath::Exp(-a1 * a1 / 2);
    double A2 = TMath::Power(n2 / TMath::Abs(a2), n2) * TMath::Exp(-a2 * a2 / 2);
    double B1 = n1 / TMath::Abs(a1) - TMath::Abs(a1);
    double B2 = n2 / TMath::Abs(a2) - TMath::Abs(a2);

    if (u < -a1) return n * A1 * TMath::Power(B1 - u, -n1);
    if (u < a2) return n * TMath::Exp(-u * u / 2);
    return n * A2 * TMath::Power(B2 + u, -n2);
}

static float transverse_momentum(bool ecal, float var, float eta) {
    return ecal ? var / std::cosh(eta) : var;
}

static std::string index_to_string(int64_t i, int64_t j) {
    return std::to_string(i) + "_"s + std::to_string(j);
}

int64_t dielectrons(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto ecal = conf->get<bool>("ecal");
    auto tag = conf->get<std::string>("tag");
    auto dcent = conf->get<std::vector<float>>("cent");

    auto mc_branches = conf->get<bool>("mc_branches");

    std::vector<std::vector<float>> scale_factors;
    for (auto const& type : { "b"s, "e"s })
        scale_factors.push_back(
            conf->get<std::vector<float>>(type + "_scales"));

    std::vector<std::vector<float>> smear_factors;
    for (auto const& type : { "b"s, "e"s })
        smear_factors.push_back(
            conf->get<std::vector<float>>(type + "_smears"));

    std::vector<std::vector<float>> ref_smear_factors;
    for (auto const& type : { "bb"s, "be"s, "ee"s })
        ref_smear_factors.push_back(
            conf->get<std::vector<float>>("ref_"s + type + "_smears"));

    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("e");
    auto e = new etree(mc_branches, false, t);

    auto pt = ecal ? e->eleEcalE : e->elePt;

    auto cents = new interval(dcent);
    auto imass = new interval("mass (GeV/c^{2})"s, 30, 60., 120.);
    std::vector<int64_t> shape = { 3, cents->size(), 2 };

    auto fm = std::bind(&interval::book<TH1F>, imass, _1, _2, _3);
    auto minv = new history<TH1F>("mass"s, "counts"s, fm, shape);

    TRandom3* gen = new TRandom3(144);

    int64_t nentries = t->GetEntries();
    for (int64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        for (int64_t j = 0; j < e->nEle; ++j) {
            if ((*e->elePt)[j] < 20)
                continue;

            int64_t is_1_barrel = passes_electron_id<
                det::barrel, wp::loose, etree>(e, j, true);
            int64_t is_1_endcap = passes_electron_id<
                det::endcap, wp::loose, etree>(e, j, true);

            if (!is_1_barrel && !is_1_endcap)
                continue;

            /* double electron invariant mass */
            for (int64_t k = j + 1; k < e->nEle; ++k) {
                if ((*e->elePt)[k] < 20)
                    continue;

                int64_t is_2_barrel = passes_electron_id<
                    det::barrel, wp::loose, etree>(e, k, true);
                int64_t is_2_endcap = passes_electron_id<
                    det::endcap, wp::loose, etree>(e, k, true);

                if (!is_2_barrel && !is_2_endcap)
                    continue;

                int64_t type_x = is_1_endcap + is_2_endcap;
                auto cent_x = cents->index_for(e->hiBin);
                int64_t charge_x = std::abs(
                    (*e->eleCharge)[j] + (*e->eleCharge)[k]) / 2;

                auto scf1 = scale_factors[is_1_endcap][cent_x];
                auto smf1 = smear_factors[is_1_endcap][cent_x] / 91.1876;
                auto sf1 = scf1 * gen->Gaus(1., smf1);

                auto scf2 = scale_factors[is_2_endcap][cent_x];
                auto smf2 = smear_factors[is_2_endcap][cent_x] / 91.1876;
                auto sf2 = scf2 * gen->Gaus(1., smf2);

                auto mass = std::sqrt(ml_invariant_mass<coords::collider>(
                    transverse_momentum(ecal, (*pt)[j], (*e->eleEta)[j]) * sf1,
                    (*e->eleEta)[j],
                    (*e->elePhi)[j],
                    0.000511f,
                    transverse_momentum(ecal, (*pt)[k], (*e->eleEta)[k]) * sf2,
                    (*e->eleEta)[k],
                    (*e->elePhi)[k],
                    0.000511f));

                float weight = mc_branches ? e->Ncoll / 1000. : 1.;
                (*minv)[x{type_x, cent_x, charge_x}]->Fill(mass, weight);
            }
        }
    }

    TF1* fits[3][cents->size()] = { 0 };

    for (int64_t i = 0; i < 3; ++i) {
        for (int64_t j = 0; j < cents->size(); ++j) {
            auto index_string = index_to_string(i, j);

            fits[i][j] = new TF1(("f_"s + index_string).data(),
                f_double_sided_crystal_ball, 60, 120, 7);

            auto parameters = conf->get<std::vector<float>>(
                "pars_"s + index_string);
            fits[i][j]->SetParameter(0, parameters[0]);
            fits[i][j]->SetParameter(1, parameters[1]);
            fits[i][j]->SetParameter(2, parameters[2]);
            fits[i][j]->SetParameter(3, parameters[3]);
            fits[i][j]->SetParameter(4, parameters[4]);
            fits[i][j]->SetParameter(5, parameters[5]);
            fits[i][j]->SetParameter(6, parameters[6]);

            (*minv)[x{i, j, 0}]->Fit(("f_"s + index_string).data(),
                "LM", "", parameters[7], parameters[8]);

            conf->set<float>("mean_"s + index_string,
                fits[i][j]->GetParameter(1));
            conf->set<float>("sigma_"s + index_string,
                fits[i][j]->GetParameter(2));
        }
    }

    std::vector<std::string> types = { "bb"s, "be"s, "ee"s };

    for (int64_t i = 0; i < 3; ++i) {
        printf("std::vector<float> %s_scales =", types[i].data());
        for (int64_t j = 0; j < cents->size(); ++j) {
            auto index_string = index_to_string(i, j);
            auto mean = conf->get<float>("mean_"s + index_string);
            printf(" %.5f", 1.f + (91.1876 - mean) / mean);
        }
        printf("\n");
    }

    for (int64_t i = 0; i < 3; ++i) {
        printf("std::vector<float> %s_smears =", types[i].data());
        for (int64_t j = 0; j < cents->size(); ++j) {
            auto ref = ref_smear_factors[i][j];
            auto sigma = conf->get<float>("sigma_"s + index_to_string(i, j));
            auto relative = ref * ref - sigma * sigma;
            auto factor = std::sqrt(std::abs(relative));
            printf(" %.5f", relative < 0 ? -factor : factor);
        }
        printf("\n");
    }

    /* lambda to display mean, sigma from fit */
    auto info_text = [&](int64_t index) {
        int64_t i = (index - 1) / 3;
        int64_t j = (index - 1) % 3;

        auto index_string = index_to_string(i, j);

        TLatex* info = new TLatex();
        info->SetTextFont(43);
        info->SetTextSize(11);

        char buffer[128];

        sprintf(buffer, "%.0f - %.0f%%", dcent[j] / 2, dcent[j + 1] / 2);
        info->DrawLatexNDC(0.675, 0.84, buffer);
        auto mean = conf->get<float>("mean_"s + index_string);
        sprintf(buffer, "mean: %.2f", mean);
        info->DrawLatexNDC(0.675, 0.78, buffer);
        auto sigma = conf->get<float>("sigma_"s + index_string);
        sprintf(buffer, "sigma: %.2f", sigma);
        info->DrawLatexNDC(0.675, 0.75, buffer);
    };

    auto hb = new pencil();
    hb->category("type", "bb", "be", "ee");
    hb->category("sign", "os", "ss");

    hb->alias("bb", "EB #otimes EB");
    hb->alias("be", "EB #otimes EE");
    hb->alias("ee", "EE #otimes EE");
    hb->alias("os", "opp. sign");
    hb->alias("ss", "same sign");

    auto c1 = new paper("mass_"s + tag, hb);
    apply_style(c1, "PbPb #sqrt{s} = 5.02 TeV"s);
    c1->legend(std::bind(coordinates, 0.135, 0.4, 0.75, 0.04));
    c1->accessory(info_text);
    c1->divide(cents->size(), 3);

    for (int64_t i = 0; i < 3; ++i) {
        for (int64_t j = 0; j < cents->size(); ++j) {
            c1->add((*minv)[x{i, j, 0}], types[i], "os");
            c1->stack((*minv)[x{i, j, 1}], types[i], "ss");
        }
    }

    hb->set_binary("sign");
    hb->sketch();
    c1->draw("pdf");

    in(output, [&]() { minv->save(tag); });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return dielectrons(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
