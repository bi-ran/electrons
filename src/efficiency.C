#include "../git/config/include/configurer.h"

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooCmdArg.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooMappedCategory.h"
#include "RooMinimizer.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooProfileLL.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TTree.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;

int efficiency(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto pdfs = conf->get<std::vector<std::string>>("pdfs");
    auto target = conf->get<std::string>("target");
    auto abscissa = conf->get<std::string>("abscissa");
    auto limits = conf->get<std::vector<float>>("limits");
    auto variables = conf->get<std::vector<std::string>>("variables");

    auto bins = std::vector<std::vector<float>>();
    for (auto const& var : variables)
        bins.push_back(conf->get<std::vector<float>>(var + "_bins"));

    auto ncpu = conf->get<int32_t>("ncpu");

    /* roofit variables */
    RooArgSet args;

    RooRealVar x(abscissa.data(), abscissa.data(), limits[0], limits[1]);

    /* target category */
    RooCategory target_category(target.data(), target.data());
    target_category.defineType("pass", 1);
    target_category.defineType("fail", 0);

    args.add(x);
    args.add(target_category);

    /* remap target category */
    RooMappedCategory target_category_map("target_category", "target_category",
        target_category, "fail");
    target_category_map.map("pass", "pass");

    /* input data */
    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("tnp");

    TFile* fout = new TFile(output, "recreate");

    auto data = new RooDataSet("data", "data", t, args, "", nullptr);
    data->addColumn(target_category);
    data->addColumn(target_category_map);

    /* fitting */
    RooRealVar efficiency("efficiency", "efficiency", 0, 1);

    auto* w = new RooWorkspace("w");
    w->import(*data);

    /* construct pdfs */
    for (auto const& pdf : pdfs)
        w->factory(pdf.data());

    w->factory("total[1,0,1e10]");
    w->factory("signal_fraction[0.9]");
    w->factory("signal_efficiency[0.9,0,1]");
    w->factory("background_efficiency[0.9,0,1]");

    w->factory("expr::signal_passing("
               "'signal_efficiency*signal_fraction*total',"
               "signal_efficiency, signal_fraction, total)");
    w->factory("expr::signal_failing("
               "'(1-signal_efficiency)*signal_fraction*total',"
               "signal_efficiency, signal_fraction, total)");
    w->factory("expr::background_passing("
               "'background_efficiency*(1-signal_fraction)*total',"
               "background_efficiency, signal_fraction, total)");
    w->factory("expr::background_failing("
               "'(1-background_efficiency)*(1-signal_fraction)*total',"
               "background_efficiency, signal_fraction, total)");

    w->factory("SUM::passing_pdf(signal_passing*signal_pdf,"
               "background_passing*background_pdf_passing)");
    w->factory("SUM::failing_pdf(signal_failing*signal_pdf,"
               "background_failing*background_pdf_failing)");

    w->factory("SIMUL::simultaneous_pdf(target_category,"
               "pass=passing_pdf, fail=failing_pdf)");

    /* set initial values */
    w->var("signal_efficiency")->setConstant(false);

    auto total_passing = w->data("data")->sumEntries(
        "target_category==target_category::pass");
    if (total_passing == 0) {
        w->var("signal_efficiency")->setVal(0.);
        w->var("signal_efficiency")->setAsymError(0., 1.);
    }

    auto total_failing = w->data("data")->sumEntries(
        "target_category==target_category::fail");
    if (total_failing == 0) {
        w->var("signal_efficiency")->setVal(1.);
        w->var("signal_efficiency")->setAsymError(-1., 0.);
    }

    auto total = total_passing + total_failing;
    w->var("total")->setVal(total);
    w->var("total")->setMax(2. * total + 10.);

    RooAbsReal* nll = w->pdf("simultaneous_pdf")->createNLL(
        *data, RooFit::Extended(true), RooFit::NumCPU(ncpu));

    RooMinimizer minimizer(*nll);
    RooMinuit minuit(*nll);

    minuit.setStrategy(1);
    minuit.setProfile(true);

    RooProfileLL prof("profile", "", *nll, *w->var("signal_efficiency"));

    /* release parameters */

    minimizer.minimize("Minuit2", "Scan");
    minuit.migrad();
    minuit.hesse();

    auto result = w->pdf("simultaneous_pdf")->fitTo(
        *data,
        RooFit::Save(true),
        RooFit::Extended(true),
        RooFit::NumCPU(ncpu),
        RooFit::Strategy(2),
        RooFit::Minos(*w->var("signal_efficiency")),
        RooFit::PrintLevel(1),
        RooFit::PrintEvalErrors(1),
        RooFit::Warnings(true));

    /* save results */
    result->Write("results");

    auto data_passing = data->reduce(RooFit::Cut(
        "target_category==target_category::pass"));
    auto data_failing = data->reduce(RooFit::Cut(
        "target_category==target_category::fail"));
    auto data_pdf = w->pdf("simultaneous_pdf");

    TCanvas* c1 = new TCanvas("c1", "", 1200, 800);
    c1->Divide(3, 2);

    c1->cd(1);
    auto p_pass = x.frame(RooFit::Title("Passing probes"));
    data_passing->plotOn(p_pass);
    data_pdf->plotOn(p_pass,
        RooFit::Slice(target_category, "pass"),
        RooFit::ProjWData(*data_passing),
        RooFit::LineColor(kGreen),
        RooFit::Components("background_passing"),
        RooFit::LineStyle(kDashed));
    data_pdf->plotOn(p_pass,
        RooFit::Slice(target_category, "pass"),
        RooFit::ProjWData(*data_passing),
        RooFit::LineColor(kGreen));
    p_pass->Draw();

    c1->cd(2);
    auto p_fail = x.frame(RooFit::Title("Failing probes"));
    data_failing->plotOn(p_fail);
    data_pdf->plotOn(p_fail,
        RooFit::Slice(target_category, "fail"),
        RooFit::ProjWData(*data_failing),
        RooFit::LineColor(kRed),
        RooFit::Components("background_failing"),
        RooFit::LineStyle(kDashed));
    data_pdf->plotOn(p_fail,
        RooFit::Slice(target_category, "fail"),
        RooFit::ProjWData(*data_failing),
        RooFit::LineColor(kRed));
    p_fail->Draw();

    c1->cd(3);
    auto p_all = x.frame(RooFit::Title("All probes"));
    data->plotOn(p_all);
    data_pdf->plotOn(p_all,
        RooFit::ProjWData(*data),
        RooFit::LineColor(kBlue),
        RooFit::Components("background_passing, background_failing"),
        RooFit::LineStyle(kDashed));
    data->plotOn(p_all,
        RooFit::ProjWData(*data),
        RooFit::LineColor(kBlue));
    p_all->Draw();

    auto text_lambda = [](std::string const& text) {
        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(11);
        l->DrawLatexNDC(0.45, 0.8, ("#chi^{2}/ndf: "s + text).data());
    };

    int32_t ndof = result->floatParsFinal().getSize();

    c1->cd(4);
    auto p_pass_pull = x.frame();
    p_pass_pull->addPlotable(p_pass->pullHist());
    auto chi2_pass = p_pass->chiSquare(ndof);
    p_pass_pull->Draw();
    text_lambda(std::to_string(chi2_pass));

    c1->cd(5);
    auto p_fail_pull = x.frame();
    p_fail_pull->addPlotable(p_fail->pullHist());
    auto chi2_fail = p_fail->chiSquare(ndof);
    p_fail_pull->Draw();
    text_lambda(std::to_string(chi2_fail));

    c1->cd(6);
    auto p_all_pull = x.frame();
    p_all_pull->addPlotable(p_all->pullHist());
    auto chi2_all = p_all->chiSquare(ndof);
    p_all_pull->Draw();
    text_lambda(std::to_string(chi2_all));

    c1->SaveAs("tnp.pdf");

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return efficiency(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
