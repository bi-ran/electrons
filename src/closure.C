#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/trunk.h"
#include "../git/tricks-and-treats/include/zip.h"

#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

static std::string index_to_string(int64_t i, int64_t j) {
    return std::to_string(i) + "_"s + std::to_string(j);
}

int closure(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto files = conf->get<std::vector<std::string>>("files");
    auto labels = conf->get<std::vector<std::string>>("labels");
    auto legends = conf->get<std::vector<std::string>>("legends");
    auto dcent = conf->get<std::vector<float>>("cent");

    auto icent = new interval(dcent);

    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    std::vector<TFile*> fs(files.size(), nullptr);
    std::vector<history<TH1F>*> hs(files.size(), nullptr);

    zip([&](TFile*& f, std::string const& file, history<TH1F>*& h,
            std::string const& label) {
        f = new TFile(file.data(), "read");
        h = new history<TH1F>(f, label);

        h->apply([](TH1* hist) { hist->Scale(1. / hist->Integral()); });
    }, fs, files, hs, labels);

    auto hratio = new history<TH1F>(*hs[0], "ratio");

    auto hb = new pencil();
    hb->category("type", "bb", "be", "ee");
    hb->category("sample", "data", "mc", "ratio");

    hb->alias("bb", "EB #otimes EB");
    hb->alias("be", "EB #otimes EE");
    hb->alias("ee", "EE #otimes EE");
    hb->alias("data", legends[0]);
    hb->alias("mc", legends[1]);
    hb->alias("ratio", legends[0] + " / " + legends[1]);

    hb->ditto("ratio", "data");

    auto ncents = icent->size();

    /* lambda to customise y-axis range per canvas */
    auto ratio_style = [&](TH1* obj, int64_t index) {
        if (index > ncents) {
            obj->GetYaxis()->SetTitle("ratio");
            obj->SetAxisRange(0.0, 2.0, "Y");
        }
    };

    /* lambda to display mean, sigma from fit */
    auto info_text = [&](int64_t index) {
        int64_t j = (index - 1) % 3;

        char buffer[128];
        sprintf(buffer, "%.0f - %.0f%%", dcent[j] / 2, dcent[j + 1] / 2);

        TLatex* info = new TLatex();
        info->SetTextFont(43);
        info->SetTextSize(12);
        info->DrawLatexNDC(0.675, 0.84, buffer);
    };

    std::vector<std::string> types = { "bb"s, "be"s, "ee"s };

    auto c1 = std::array<paper*, 3>();
    for (int64_t i = 0; i < 3; ++i) {
        c1[i] = new paper("closure_"s + types[i], hb);
        apply_style(c1[i], system + " #sqrt{s} = 5.02 TeV"s, 0., 0.3);
        c1[i]->legend(std::bind(coordinates, 0.135, 0.4, 0.75, 0.04));
        c1[i]->accessory(info_text);
        c1[i]->accessory(std::bind(line_at, _1, 1, 60, 120));
        c1[i]->jewellery(ratio_style);

        c1[i]->add(ncents * 2);
        c1[i]->divide(ncents, 2);

        for (int64_t j = 0; j < ncents; ++j) {
            auto index = hratio->index_for(x{i, j, 0});
            (*hratio)[index]->Divide((*hs[1])[index]);

            for (auto h : { hs[0], hs[1], hratio })
                (*h)[index]->GetFunction(("f_"s + index_to_string(i, j)).data())
                    ->SetBit(TF1::kNotDraw);

            c1[i]->stack(j + 1, (*hs[1])[index], types[i], "mc");
            c1[i]->stack(j + 1, (*hs[0])[index], types[i], "data");
            c1[i]->stack(ncents + j + 1, (*hratio)[index], types[i], "ratio");
        }
    }

    hb->set_binary("sample");
    hb->sketch();

    for (auto c : c1) { c->draw("pdf"); }

    in(output, [&]() {
        hs[0]->save(tag);
        hs[1]->save(tag);
        hratio->save(tag);
    });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return closure(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
