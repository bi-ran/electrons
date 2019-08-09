#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;

static std::string index_to_string(int64_t i, int64_t j) {
    return std::to_string(i) + "_"s + std::to_string(j);
}

int closure(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto data = conf->get<std::string>("data");
    auto mc = conf->get<std::string>("mc");
    auto cent = conf->get<std::vector<float>>("cent");

    TH1::SetDefaultSumw2();

    TFile* fd = new TFile(data.data(), "read");
    TFile* fm = new TFile(mc.data(), "read");

    std::vector<std::string> types = { "bb"s, "be"s, "ee"s };

    auto cents = std::make_shared<interval>(cent);

    auto hdata = new history(fd, "data_scaled_mass");
    auto hmc = new history(fm, "mc_scaled_smeared_mass");

    auto normalise = [](TH1* h) { h->Scale(1. / h->Integral()); };
    hdata->apply(normalise); hmc->apply(normalise);

    auto hratio = new history(*hdata, "ratio");

    auto hb = new pencil();
    hb->category("type", "bb", "be", "ee");
    hb->category("sample", "data", "mc", "ratio");

    hb->alias("bb", "EB #otimes EB");
    hb->alias("be", "EB #otimes EE");
    hb->alias("ee", "EE #otimes EE");
    hb->alias("mc", "MC");
    hb->alias("ratio", "data / MC");

    hb->ditto("ratio", "data");

    auto ncents = cents->size();

    /* lambda to customise y-axis range per canvas */
    auto formatter = [&](TH1* obj, int64_t index) {
        obj->SetStats(0);
        obj->SetMarkerSize(0.84);
        obj->GetXaxis()->CenterTitle();
        obj->GetYaxis()->CenterTitle();

        if (index > ncents) {
            obj->GetYaxis()->SetTitle("ratio");
            obj->SetAxisRange(0.0, 2.0, "Y");
        }
    };

    /* lambda to display mean, sigma from fit */
    auto info_text = [&](int64_t index) {
        int64_t j = (index - 1) % 3;

        char buffer[128];
        sprintf(buffer, "%.0f - %.0f%%", cent[j] / 2, cent[j + 1] / 2);

        TLatex* info = new TLatex();
        info->SetTextFont(43);
        info->SetTextSize(11);
        info->DrawLatexNDC(0.675, 0.84, buffer);
    };

    /* lambda for line at unity */
    auto line_at_unity = [&](int64_t index) {
        if (index <= ncents) { return; }

        TLine* l1 = new TLine(60, 1., 120, 1.);
        l1->SetLineStyle(7);
        l1->Draw();
    };

    auto c1 = std::array<paper*, 3>();
    for (int64_t i = 0; i < 3; ++i) {
        c1[i] = new paper("closure_"s + types[i], hb);
        apply_default_style(c1[i],"PbPb #sqrt{s} = 5.02 TeV"s, 0., 0.3);
        c1[i]->legend(std::bind(coordinates, 0.135, 0.4, 0.75, 0.04));
        c1[i]->jewellery(formatter);
        c1[i]->accessory(info_text);
        c1[i]->accessory(line_at_unity);

        c1[i]->add(ncents * 2);
        c1[i]->divide(ncents, 2);

        for (int64_t j = 0; j < ncents; ++j) {
            auto index = hratio->index_for(x{i, j, 0});
            (*hratio)[index]->Divide((*hmc)[index]);

            for (auto h : { hmc, hdata, hratio })
                (*h)[index]->GetFunction(("f_"s + index_to_string(i, j)).data())
                    ->SetBit(TF1::kNotDraw);

            c1[i]->stack(j + 1, (*hmc)[index], types[i], "mc");
            c1[i]->stack(j + 1, (*hdata)[index], types[i], "data");
            c1[i]->stack(ncents + j + 1, (*hratio)[index], types[i], "ratio");
        }
    }

    hb->set_binary("sample");
    hb->sketch();

    for (auto c : c1) { c->draw("pdf"); }

    TFile* fout = new TFile(output, "recreate");

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("usage: %s [config] [output]\n", argv[0]);
        return 1;
    }

    return closure(argv[1], argv[2]);
}
