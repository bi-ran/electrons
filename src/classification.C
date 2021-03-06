#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/zip.h"

#include "TCut.h"
#include "TFile.h"
#include "TH1.h"
#include "TMarker.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Tools.h"
#include "TMVA/efficiencies.h"
#include "TMVA/variables.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>

using namespace std::literals::string_literals;
using namespace std::placeholders;

int outliers(configurer* conf) {
    auto signal = conf->get<std::string>("signal");
    auto variables = conf->get<std::vector<std::string>>("variables");
    auto type = conf->get<std::vector<uint32_t>>("type");
    auto target = conf->get<float>("target");

    auto count = static_cast<int64_t>(variables.size());
    std::vector<double> lower(count, std::numeric_limits<float>::lowest());
    std::vector<double> upper(count, std::numeric_limits<float>::max());

    TFile* f = new TFile(signal.data(), "read");
    TTree* t = (TTree*)f->Get("e");

    int64_t elements = t->Draw("elePt", "", "goff");
    t->SetEstimate(elements);

    zip([&](std::string const& variable, double& lower_, double& upper_,
            uint32_t type_) {
        t->Draw(variable.data(), "", "goff");
        auto data = t->GetV1();

        std::vector<double> values(data, data + elements);
        std::sort(std::begin(values), std::end(values));

        switch (type_) {
            case 0:
                upper_ = values[target * elements];
                break;
            case 1:
                lower_ = values[(1. - target) * elements];
                break;
            case 2:
                lower_ = values[((1. - target) / 2.) * elements];
                upper_ = values[((1. + target) / 2.) * elements];
                break;
        }
    }, variables, lower, upper, type);

    conf->set("lower", std::move(lower));
    conf->set("upper", std::move(upper));

    return 0;
}

int classify(configurer* conf, std::string const& output,
             std::string const& id, float efficiency) {
    auto signal = conf->get<std::string>("signal");
    auto background = conf->get<std::string>("background");
    auto base = conf->get<std::string>("base");
    auto variables = conf->get<std::vector<std::string>>("variables");
    auto lower = conf->get<std::vector<double>>("lower");
    auto upper = conf->get<std::vector<double>>("upper");
    auto type = conf->get<std::vector<uint32_t>>("type");
    auto weight = conf->get<std::string>("weight");

    auto nsig_train = conf->get<int32_t>("nsig_train");
    auto nbkg_train = conf->get<int32_t>("nbkg_train");

    auto tag = conf->get<std::string>("tag");

    /* 0: barrel, 1: endcap */
    auto options = conf->get<int32_t>("options");

    /* consistency checks */
    auto size = variables.size();
    if (lower.size() != size || upper.size() != size || type.size() != size) {
        printf("  inconsistent sizes\n"); exit(1); }

    zip([&](std::string const& variable, double lower_, double upper_) {
        if (lower_ >= upper_) {
            printf("  limits on %s: %f >= %f\n",
                variable.data(), lower_, upper_);
            exit(1);
        }
    }, variables, lower, upper);

    /* setup */
    TMVA::Tools::Instance();

    TFile* fsig = new TFile(signal.data(), "read");
    TTree* tsig = (TTree*)fsig->Get("e");
    TFile* fbkg = new TFile(background.data(), "read");
    TTree* tbkg = (TTree*)fbkg->Get("e");

    TFile* fout = TFile::Open(output.data(), "recreate");

    TMVA::Factory* factory = new TMVA::Factory(
        "TMVAClassification",
        fout,
        "!V:!Silent:Color:DrawProgressBar:"
        "Transformations=I:"
        "AnalysisType=Classification");

    TMVA::DataLoader* loader = new TMVA::DataLoader((tag + "_"s + id).data());
    for (auto const& var : variables)
        loader->AddVariable(var.data(), 'F');

    loader->AddSignalTree(tsig, 1.);
    loader->AddBackgroundTree(tbkg, 1.);

    if (!weight.empty()) {
        loader->SetSignalWeightExpression(weight.data());
        loader->SetBackgroundWeightExpression(weight.data());
    }

    /* selections */
    TCut minpt = "elePt > 20";
    TCut basic = base.data();

    TCut barrel = "abs(eleSCEta) < 1.442";
    TCut endcap = "abs(eleSCEta) > 1.556 && abs(eleSCEta) < 2.4";
    TCut region = options ? endcap : barrel;

    TCut match = "gen_index != -1 && abs(mcMomPID[gen_index]) < 25";
    TCut unmatch = "gen_index == -1";
    TCut nonprompt = "gen_index != -1 && abs(mcMomPID[gen_index]) > 24";

    TCut sig_sel = minpt && basic && match && region;
    TCut bkg_sel = minpt && basic && (unmatch || nonprompt) && region;

    auto train_options = "nTrain_Signal="s + std::to_string(nsig_train)
        + ":nTrain_Background="s + std::to_string(nbkg_train)
        + ":SplitMode=Random:NormMode=NumEvents:!V"s;

    loader->PrepareTrainingAndTestTree(
        sig_sel, bkg_sel, train_options.data());

    /* settings, variable limits */
    static std::string const bounds[3] = { "FMin"s, "FMax"s, "NotEnforced"s };

    auto settings = "!H:!V:FitMethod=GA:EffSel:Steps=30:Cycles=3:PopSize=400:"s
        + "SC_steps=10:SC_rate=5:SC_factor=0.95:"s;
    for (int64_t i = 0; i < static_cast<int64_t>(size); ++i) {
        settings += "VarProp["s + std::to_string(i) + "]="s
            + bounds[type[i]] + ":"s;
        settings += "CutRangeMin["s + std::to_string(i) + "]="s
            + std::to_string(lower[i]) + ":"s;
        settings += "CutRangeMax["s + std::to_string(i) + "]="s
            + std::to_string(upper[i]) + ":"s;
    }

    /* work */
    factory->BookMethod(loader, TMVA::Types::kCuts, "CutsGA", settings.data());

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    fout->Close();

    /* extract and save cuts */
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;

    auto method = dynamic_cast<TMVA::MethodCuts*>(
        factory->GetMethod((tag + "_"s + id).data(), "CutsGA"));
    method->GetCuts(efficiency, lower_bounds, upper_bounds);

    conf->set(id + "_lower"s, std::vector<double>(lower_bounds));
    conf->set(id + "_upper"s, std::vector<double>(upper_bounds));

    /* update limits */
    conf->unset<std::vector<double>>("lower");
    conf->set("lower", std::move(lower_bounds));
    conf->unset<std::vector<double>>("upper");
    conf->set("upper", std::move(upper_bounds));

    delete factory;
    delete loader;

    return 0;
}

std::pair<float, float> evaluate(configurer* conf, std::string const& id) {
    auto signal = conf->get<std::string>("signal");
    auto background = conf->get<std::string>("background");
    auto base = conf->get<std::string>("base");
    auto variables = conf->get<std::vector<std::string>>("variables");
    auto type = conf->get<std::vector<uint32_t>>("type");
    auto lower = conf->get<std::vector<double>>(id + "_lower");
    auto upper = conf->get<std::vector<double>>(id + "_upper");
    auto weight = conf->get<std::string>("weight");

    /* 0: barrel, 1: endcap */
    auto options = conf->get<int32_t>("options");

    /* evaluate signal, background efficiencies */
    TFile* fsig = new TFile(signal.data(), "read");
    TTree* tsig = (TTree*)fsig->Get("e");
    TFile* fbkg = new TFile(background.data(), "read");
    TTree* tbkg = (TTree*)fbkg->Get("e");

    TCut minpt = "elePt > 20";
    TCut basic = base.data();

    TCut barrel = "abs(eleSCEta) < 1.442";
    TCut endcap = "abs(eleSCEta) > 1.556 && abs(eleSCEta) < 2.4";
    TCut region = options ? endcap : barrel;

    TCut match = "gen_index != -1 && abs(mcMomPID[gen_index]) < 25";
    TCut unmatch = "gen_index == -1";
    TCut nonprompt = "gen_index != -1 && abs(mcMomPID[gen_index]) > 24";

    TCut sig_sel = minpt && basic && match && region;
    TCut bkg_sel = minpt && basic && (unmatch || nonprompt) && region;

    if (weight.empty()) { weight = "1"s; }
    TCut w(weight.data());

    std::string id_string;
    zip([&](std::string const& variable, double lower_, double upper_,
            uint32_t type_) {
        id_string += "&&"s;
        if (type_ != 0) {
            id_string += std::to_string(lower_) + "<"s + variable; }
        if (type_ == 2) {
            id_string += "&&"s; }
        if (type_ != 1) {
            id_string += std::to_string(upper_) + ">"s + variable; }
    }, variables, lower, upper, type);
    id_string.erase(0, 2);

    TCut id_sel(id_string.data());

    int64_t totalsig = tsig->Draw("elePt", sig_sel * w, "goff");
    int64_t selsig = tsig->Draw("elePt", (sig_sel && id_sel) * w, "goff");

    int64_t totalbkg = tbkg->Draw("elePt", bkg_sel * w, "goff");
    int64_t selbkg = tbkg->Draw("elePt", (bkg_sel && id_sel) * w, "goff");

    float sigeff = static_cast<float>(selsig) / totalsig;
    float bkgrej = 1.f - static_cast<float>(selbkg) / totalbkg;

    printf("\n  [ %.4f / %.4f ]\n\n", sigeff, bkgrej);

    return std::make_pair(sigeff, bkgrej);
}

void draw(configurer* conf, std::string const& output) {
    auto ids = conf->get<std::vector<std::string>>("ids");
    auto effs = conf->get<std::vector<float>>("effs");
    auto cols = conf->get<std::vector<std::string>>("cols");
    auto tag = conf->get<std::string>("tag");

    /* extract roc curve from output file */
    TFile* f = new TFile(output.data(), "read");
    auto roc = (TH1F*)f->Get(
        (tag + "_"s + ids[0] + "/Method_CutsGA/CutsGA/MVA_CutsGA_rejBvsS").data());

    /* lambda to format histogram approriately */
    auto formatter = [](TH1* obj) {
        obj->SetStats(0);
        obj->SetMarkerSize(0.4);
        obj->SetMarkerStyle(20);
        obj->SetAxisRange(0., 1., "Y");
        obj->SetTitle(";signal efficiency;background rejection");
        obj->GetXaxis()->CenterTitle();
        obj->GetYaxis()->CenterTitle();
    };

    auto label = output;
    auto ext = label.find(".root");
    if (ext != std::string::npos)
        label.erase(std::begin(label) + ext, std::end(label));

    auto c1 = new paper("working-points-"s + label);
    apply_style(c1, "pp #sqrt{s} = 5.02 TeV"s, formatter);
    c1->set(paper::key);

    c1->add(roc);
    c1->adjust(roc, "l", "");

    /* lambda to draw marker at working points */
    auto mark = [&](int64_t, std::pair<float, float> const& wp, int32_t col) {
        TMarker* marker = new TMarker(wp.first, wp.second, 21);
        marker->SetMarkerSize(0.8);
        marker->SetMarkerColor(col);
        marker->Draw("same");
    };

    std::unordered_map<std::string, int32_t> cmap;
    zip([&](std::string const& id, std::string const& col) {
        cmap[id] = TColor::GetColor(col.data()); }, ids, cols);

    /* print selections and evaluate efficiencies */
    auto variables = conf->get<std::vector<std::string>>("variables");
    auto type = conf->get<std::vector<uint32_t>>("type");

    for (auto const& id : ids) {
        auto lower = conf->get<std::vector<double>>(id + "_lower"s);
        auto upper = conf->get<std::vector<double>>(id + "_upper"s);

        printf("  %s:\n", id.data());
        zip([&](std::string const& variable, double lower_, double upper_,
                uint32_t type_) {
            printf("%24s: [ ", variable.data());
            if (type_ == 0) { printf("    -inf, "); }
            else { printf("%8.5f, ", lower_); }
            if (type_ == 1) { printf("+inf     ]\n"); }
            else { printf("%8.5f ]\n", upper_); }
        }, variables, lower, upper, type);

        printf("\n");

        auto wp = evaluate(conf, id);
        c1->accessory(std::bind(mark, _1, wp, cmap[id]));
    }

    c1->draw("pdf");
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("usage: %s [config] [output]\n", argv[0]);
        return 1;
    }

    auto conf = new configurer(argv[1]);

    auto stage = conf->get<uint64_t>("stage");
    auto ids = conf->get<std::vector<std::string>>("ids");
    auto effs = conf->get<std::vector<float>>("effs");
    auto tag = conf->get<std::string>("tag");

    auto base_stub = tag + "_"s + ids[0];
    auto base_tag = tag + "_"s + ids[0] + "_"s + argv[2];

    switch (stage) {
        case 0: goto _stage0;
        case 1: goto _stage1;
    }

_stage0:
    outliers(conf);

    zip([&](std::string const& id, float eff) {
        auto full_tag = tag + "_"s + id + "_"s + argv[2];
        classify(conf, full_tag, id, eff);
    }, ids, effs);

    draw(conf, base_tag);

_stage1:
    TMVA::efficiencies(base_stub.data(), base_tag.data(), 2, true);
    TMVA::variables(base_stub.data(), base_tag.data(),
                    "InputVariables_Id", "comparison", false, true);

    return 0;
}
