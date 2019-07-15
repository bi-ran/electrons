#include "../include/etree.h"

#include "../git/config/include/configurer.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

using namespace std::literals::string_literals;

int reweight(configurer* conf, std::string const& output) {
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

    TFile* fout = new TFile(output.data(), "recreate");

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

int append(configurer* conf, std::string const& output) {
    auto input = conf->get<std::string>("input");
    auto tree = conf->get<std::string>("tree");

    TFile* fw = new TFile(output.data(), "read");
    TH2F* hweights = (TH2F*)fw->Get("hweights");

    auto xaxis = hweights->GetXaxis();
    auto yaxis = hweights->GetYaxis();

    TFile* f = new TFile(input.data(), "update");
    TTree* t = (TTree*)f->Get(tree.data());
    auto e = new etree(0, 0, t);

    std::vector<float> weight;
    TBranch* b = t->Branch("weight", &weight);

    int64_t nentries = t->GetEntries();
    for (int64_t i = 0; i < nentries; ++i) {
        weight.clear();

        t->GetEntry(i);
        for (int64_t j = 0; j < e->nEle; ++j)
            weight.push_back(hweights->GetBinContent(
                xaxis->FindBin((*e->elePt)[j]),
                yaxis->FindBin((*e->eleEta)[j])));

        b->Fill();
    }

    t->Write("", TObject::kOverwrite);
    f->Write("", TObject::kOverwrite);
    f->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("usage: %s [config] [output]\n", argv[0]);
        return 1;
    }

    auto conf = new configurer(argv[1]);

    auto apply = conf->get<bool>("apply");

    reweight(conf, argv[2]);
    if (apply) { append(conf, argv[2]); }

    return 0;
}

