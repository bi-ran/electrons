#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLine.h"

using namespace std::placeholders;
using namespace std::literals::string_literals;

int64_t scale_factors(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto dir = conf->get<std::string>("dir");
    auto input_mc = conf->get<std::vector<std::string>>("input_mc");
    auto input_data = conf->get<std::vector<std::string>>("input_data");
    auto categories = conf->get<std::vector<std::string>>("categories");
    auto var = conf->get<std::string>("var");

    auto hb = new pencil();
    hb->category("sample", "data", "mc");
    hb->category("type", "barrel", "endcap", "incl");

    hb->alias("mc", "MC");
    hb->alias("incl", "");

    auto frame_formatter = [&](TH1* obj) {
        obj->SetTitle("");
        obj->GetXaxis()->SetTitleFont(43);
        obj->GetXaxis()->SetTitleSize(15);
        obj->GetYaxis()->SetTitleFont(43);
        obj->GetYaxis()->SetTitleSize(15);
        obj->GetXaxis()->SetLabelFont(43);
        obj->GetXaxis()->SetLabelSize(12);
        obj->GetYaxis()->SetLabelFont(43);
        obj->GetYaxis()->SetLabelSize(12);
        obj->GetXaxis()->SetTitle(var.data());
        obj->GetYaxis()->SetTitle("efficiency");
        obj->GetXaxis()->CenterTitle();
        obj->GetYaxis()->CenterTitle();
        obj->SetAxisRange(0.5, 1.2, "Y");
    };

    auto graph_formatter = [](TGraph* obj) {
        obj->SetMarkerSize(0.84);
    };

    auto c1 = new paper("scale_factors_"s + output, hb);
    apply_default_style(c1, "pp #sqrt{s} = 5.02 TeV"s, 0., 1.);
    c1->legend(std::bind(coordinates, 0.54, 0.9, 0.84, 0.04));
    c1->format(frame_formatter);
    c1->format(graph_formatter);

    std::vector<double> low_edges = { 0 };
    std::vector<double> high_edges = { 0 };

    auto panels = static_cast<int64_t>(input_data.size());
    for (int64_t i = 0; i < panels; ++i) {
        TFile* fm = new TFile((dir + "/"s + input_mc[i]).data(), "read");
        TFile* fd = new TFile((dir + "/"s + input_data[i]).data(), "read");

        TCanvas* cm = (TCanvas*)fm->Get(
            ("mc/efficiency/fit_eff_plots/probe_"s + var + "_PLOT"s).data());
        TCanvas* cd = (TCanvas*)fd->Get(
            ("data/efficiency/fit_eff_plots/probe_"s + var + "_PLOT"s).data());

        auto hframe = (TH1F*)cd->GetPrimitive("frame");

        low_edges.push_back(hframe->GetBinLowEdge(1));
        high_edges.push_back(hframe->GetBinLowEdge(hframe->GetNbinsX() + 1));

        auto gmc = (TGraphAsymmErrors*)cm->GetPrimitive("hxy_fit_eff");
        auto gdata = (TGraphAsymmErrors*)cd->GetPrimitive("hxy_fit_eff");

        c1->add(hframe);
        c1->stack(gmc, "mc", categories[i]);
        c1->stack(gdata, "data", categories[i]);
    }

    auto line_at_unity = [&](int64_t index) {
        TLine* l1 = new TLine(low_edges[index], 1., high_edges[index], 1.);
        l1->SetLineStyle(7);
        l1->Draw();
    };

    c1->accessory(line_at_unity);

    hb->ditto("incl", "barrel");
    hb->sketch();

    c1->draw("pdf");

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("usage: %s [config] [output]\n", argv[0]);
        return 1;
    }

    return scale_factors(argv[1], argv[2]);
}
