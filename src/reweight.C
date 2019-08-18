#include "../git/config/include/configurer.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

using namespace std::literals::string_literals;

int reweight(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto target = conf->get<std::string>("target");
    auto tree = conf->get<std::string>("tree");

    auto variables = conf->get<std::vector<std::string>>("variables");
    auto bins1 = conf->get<std::vector<float>>("bins1");
    auto bins2 = conf->get<std::vector<float>>("bins2");

    auto selection = conf->get<std::string>("selection");

    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get(tree.data());

    TFile* ft = new TFile(target.data(), "read");
    TTree* tt = (TTree*)ft->Get(tree.data());

    TFile* fout = new TFile(output, "recreate");

    const int64_t nvars = static_cast<int64_t>(variables.size());
    /* support for 1d reweighting not implemented yet */
    if (nvars != 2) { return 1; }

    TH2F* hinput = new TH2F("hinput", "",
        bins1.size() - 1, bins1.data(), bins2.size() - 1, bins2.data());
    TH2F* htarget = new TH2F("htarget", "",
        bins1.size() - 1, bins1.data(), bins2.size() - 1, bins2.data());

    auto varexp = variables[1] + ":"s + variables[0];

    t->Draw((varexp + ">>hinput"s).data(), selection.data(), "goff");
    tt->Draw((varexp + ">>htarget"s).data(), selection.data(), "goff");

    TH2F* hweights = (TH2F*)htarget->Clone("hweights");
    hweights->Divide(hinput);

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return reweight(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}

