#include "../include/lambdas.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLine.h"

using namespace std::literals::string_literals;

int64_t scale_factors(std::string const& inputm, std::string const& inputd,
                      std::string const& var, std::string const& output) {
    TFile* fm = new TFile(inputm.data(), "read");
    TFile* fd = new TFile(inputd.data(), "read");

    TCanvas* cm = (TCanvas*)fm->Get(
        ("mc/efficiency/fit_eff_plots/probe_"s + var + "_PLOT"s).data());
    TCanvas* cd = (TCanvas*)fd->Get(
        ("data/efficiency/fit_eff_plots/probe_"s + var + "_PLOT"s).data());

    auto hframe = (TH1F*)cd->GetPrimitive("frame");

    auto gmc = (TGraphAsymmErrors*)cm->GetPrimitive("hxy_fit_eff");
    auto gdata = (TGraphAsymmErrors*)cd->GetPrimitive("hxy_fit_eff");

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

    auto line_at_unity = [&](int64_t) {
        double low_edge = hframe->GetBinLowEdge(1);
        double high_edge = hframe->GetBinLowEdge(hframe->GetNbinsX() + 1);

        TLine* l1 = new TLine(low_edge, 1., high_edge, 1.);
        l1->SetLineStyle(7);
        l1->Draw();
    };

    auto c1 = new paper("scale_factors_"s + output, hb);
    apply_default_style(c1, "pp #sqrt{s} = 5.02 TeV"s, 0., 1.);
    c1->legend(std::bind(coordinates, 0.54, 0.9, 0.84, 0.04));
    c1->format(frame_formatter);
    c1->format(graph_formatter);
    c1->accessory(line_at_unity);

    c1->add(hframe);
    c1->stack(gmc, "mc", "incl");
    c1->stack(gdata, "data", "incl");

    hb->ditto("incl", "barrel");
    hb->sketch();

    c1->draw("pdf");

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("usage: %s [mc] [data] [var] [output]\n", argv[0]);
        return 1;
    }

    return scale_factors(argv[1], argv[2], argv[3], argv[4]);
}
