#include "../git/config/include/configurer.h"

#include "TCut.h"
#include "TFile.h"
#include "TH1.h"
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

using namespace std::literals::string_literals;

int outliers(configurer* conf) {
    auto signal = conf->get<std::string>("signal");
    auto variables = conf->get<std::vector<std::string>>("variables");
    auto type = conf->get<std::vector<uint32_t>>("type");
    auto target = conf->get<float>("target");

    auto count = static_cast<int64_t>(variables.size());

    for (auto const& t : type)
        if (t > 2) { printf("  type must belong in [0, 2]\n"); exit(1); }

    if (conf->test<std::vector<double>>("lower")
        || conf->test<std::vector<double>>("upper")) { return 0; }

    std::vector<double> lower(count, std::numeric_limits<float>::lowest());
    std::vector<double> upper(count, std::numeric_limits<float>::max());

    TFile* f = new TFile(signal.data(), "read");
    TTree* t = (TTree*)f->Get("e");

    int64_t elements = t->Draw("elePt", "", "goff");

    for (int64_t i = 0; i < count; ++i) {
        auto const& var = variables[i];

        t->SetEstimate(elements);
        t->Draw(var.data(), "", "goff");
        auto data = t->GetV1();

        std::vector<double> values(data, data + elements);
        std::sort(std::begin(values), std::end(values));

        switch (type[i]) {
            case 0:
                upper[i] = values[target * elements];
                break;
            case 1:
                lower[i] = values[(1. - target) * elements];
                break;
            case 2:
                lower[i] = values[((1. - target) / 2.) * elements];
                upper[i] = values[((1. + target) / 2.) * elements];
                break;
        }
    }

    conf->set("lower", std::move(lower));
    conf->set("upper", std::move(upper));

    return 0;
}

int classify(configurer* conf, std::string const& output,
             std::string const& id, float efficiency) {
    auto signal = conf->get<std::string>("signal");
    auto background = conf->get<std::string>("background");
    auto variables = conf->get<std::vector<std::string>>("variables");
    auto lower = conf->get<std::vector<double>>("lower");
    auto upper = conf->get<std::vector<double>>("upper");
    auto type = conf->get<std::vector<uint32_t>>("type");
    auto weight = conf->get<std::string>("weight");

    /* 0: barrel, 1: endcap */
    auto options = conf->get<int32_t>("options");

    /* consistency checks */
    auto size = variables.size();
    if (lower.size() != size || upper.size() != size || type.size() != size) {
        printf("  inconsistent sizes\n"); exit(1); }

    auto count = static_cast<int64_t>(size);
    for (int64_t i = 0; i < count; ++i) {
        if (lower[i] >= upper[i]) {
            printf("  limits on %s: %f >= %f\n",
                variables[i].data(), lower[i], upper[i]);
            exit(1);
        }
    }

    /* setup */
    TMVA::Tools::Instance();

    TFile* fsig = new TFile(signal.data(), "read");
    TTree* tsig = (TTree*)fsig->Get("e");
    TFile* fbkg = new TFile(background.data(), "read");
    TTree* tbkg = (TTree*)fbkg->Get("e");

    TFile* fout = TFile::Open((id + "_"s + output).data(), "recreate");

    TMVA::Factory* factory = new TMVA::Factory(
        "TMVAClassification",
        fout,
        "!V:!Silent:Color:DrawProgressBar:"
        "Transformations=I:"
        "AnalysisType=Classification");

    TMVA::DataLoader* loader = new TMVA::DataLoader("dataset");
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

    TCut barrel = "abs(eleSCEta) < 1.442";
    TCut endcap = "abs(eleSCEta) > 1.556 && abs(eleSCEta) < 2.4";
    TCut region = options ? endcap : barrel;

    TCut match = "gen_index != -1 && abs(mcMomPID[gen_index]) < 25";
    TCut unmatch = "gen_index == -1";
    TCut nonprompt = "gen_index != -1 && abs(mcMomPID[gen_index]) > 24";

    TCut sig_sel = minpt && match && region;
    TCut bkg_sel = minpt && (unmatch || nonprompt) && region;

    loader->PrepareTrainingAndTestTree(
        sig_sel, bkg_sel,
        "nTrain_Signal=10000:nTrain_Background=10000:"
        "SplitMode=Random:NormMode=NumEvents:!V");

    /* settings, variable limits */
    static std::string const bounds[3] = { "FMin"s, "FMax"s, "NotEnforced"s };

    auto settings = "!H:!V:FitMethod=GA:EffSel:Steps=30:Cycles=3:PopSize=400:"s
        + "SC_steps=10:SC_rate=5:SC_factor=0.95:"s;
    for (int64_t i = 0; i < count; ++i) {
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
        factory->GetMethod("dataset", "CutsGA"));
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

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("usage: %s [config] [output]\n", argv[0]);
        return 1;
    }

    auto conf = new configurer(argv[1]);

    outliers(conf);

    classify(conf, argv[2], "veto", 0.95);
    classify(conf, argv[2], "loose", 0.9);
    classify(conf, argv[2], "medium", 0.8);
    classify(conf, argv[2], "tight", 0.7);

    auto variables = conf->get<std::vector<std::string>>("variables");
    auto type = conf->get<std::vector<uint32_t>>("type");
    auto count = static_cast<int64_t>(variables.size());

    for (auto const& id : { "veto"s, "loose"s, "medium"s, "tight"s }) {
        auto lower = conf->get<std::vector<double>>(id + "_lower"s);
        auto upper = conf->get<std::vector<double>>(id + "_upper"s);

        printf("  %s:\n", id.data());
        for (int64_t i = 0; i < count; ++i) {
            printf("%24s: [ ", variables[i].data());
            if (type[i] == 0) { printf("    -inf, "); }
            else { printf("%8.5f, ", lower[i]); }
            if (type[i] == 1) { printf("+inf     ]\n"); }
            else { printf("%8.5f ]\n", upper[i]); }
        }

        printf("\n");
    }

    TMVA::efficiencies("dataset", ("veto_"s + argv[2]).data(), 2, true);
    TMVA::variables("dataset", ("veto_"s + argv[2]).data(),
                    "InputVariables_Id", "comparison", false, true);

    return 0;
}
