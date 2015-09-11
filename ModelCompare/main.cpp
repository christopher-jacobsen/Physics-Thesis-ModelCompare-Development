//
//  main.cpp
//  ModelCompare
//
//  Created by Christopher Jacobsen on 31/08/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "ModelCompare.h"
#include "RootUtil.h"
#include "common.h"

////////////////////////////////////////////////////////////////////////////////

using namespace RootUtil;
using namespace ModelCompare;

////////////////////////////////////////////////////////////////////////////////

static const ObservableVector Observables1 =
{
    { "PTZ", "P_{T}(Z)",  100,     0,  500, "P_{T}(Z) [GeV/c]", "Events per 5 GeV/c",       [](const HepMC::GenVertex & s, TH1D & h, double w) { FillHistPT( s, h, w, 24);     } },
    { "MWZ", "M(WZ)",     100,     0, 1500, "M(WZ) [GeV/c^2]",  "Events per 15 GeV/c^2",    [](const HepMC::GenVertex & s, TH1D & h, double w) { FillHistM2( s, h, w, 24, 23); } },
    { "ETZ", "#eta(Z)",   100,   -10,   10, "#eta(Z)",          "Events per bin",           [](const HepMC::GenVertex & s, TH1D & h, double w) { FillHistEta(s, h, w, 24);     } },
    { "PHZ", "#phi(Z)",   100, -M_PI, M_PI, "#phi(Z)",          "Events per bin",           [](const HepMC::GenVertex & s, TH1D & h, double w) { FillHistPhi(s, h, w, 24);     } },
};

////////////////////////////////////////////////////////////////////////////////

static const ModelFileVector Models_1E4 =
{
    { "SM_211_1E4.hepmc2g",                 "SM_211",       "SM (2.1.1)",                                   18.3748 },
    { "SM_220_1E4.hepmc2g",                 "SM_220",       "SM (2.2.0)",                                   18.2613 },
    { "SM_AGC_211_1E4.hepmc2g",             "SM_AGC",       "SM-AGC (2.1.1)",                               18.6957 },
    { "SM_UFO_220_1E4.hepmc2g",             "SM_UFO",       "SM-UFO (2.2.0)",                               18.2796 },

    { "EFT_220_cWWW_3E-5_1E4.hepmc2g",      "EFT_cWWW",     "EFT cWWW = 3E-5",                              30.3160 },
    { "EFT_220_cW_5E-5_1E4.hepmc2g",        "EFT_cW",       "EFT cW = 5E-5",                                30.7824 },
    { "EFT_220_cB_9E-4_1E4.hepmc2g",        "EFT_cB",       "EFT cB = 9E-4",                                31.4724 },
    { "EFT_220_all_1E4.hepmc2g",            "EFT_all",      "EFT (all)",                                    42.9091 },

    { "AGC_211_lambda_127E-3_1E4.hepmc2g",  "AGC_lambda",   "AGC #lambda_{#gamma/Z} = 0.127",               30.5172 },
    { "AGC_211_g1_121E-2_1E4.hepmc2g",      "AGC_g1",       "AGC #Deltag1_{Z} = 0.208",                     30.6187 },
    { "AGC_211_kappa_391E-2_1E4.hepmc2g",   "AGC_kappa",    "AGC #Delta#kappa_{#gamma/Z} = 2.91/-0.83",     31.5387 },
    { "AGC_211_all_1E4.hepmc2g",            "AGC_all",      "AGC (all)",                                    43.1840 },
};

static const ModelFileVector Models_1E5 =
{
    { "SM_211_1E5.hepmc2g",                 "SM_211",       "SM (2.1.1)",                                   18.4996 },
    { "SM_220_1E5.hepmc2g",                 "SM_220",       "SM (2.2.0)",                                   18.4850 },
    { "SM_AGC_211_1E5.hepmc2g",             "SM_AGC",       "SM-AGC (2.1.1)",                               18.5791 },
    { "SM_UFO_220_1E5.hepmc2g",             "SM_UFO",       "SM-UFO (2.2.0)",                               18.4768 },

    { "EFT_220_cWWW_3E-5_1E5.hepmc2g",      "EFT_cWWW",     "EFT cWWW = 3E-5",                              30.4574 },
    { "EFT_220_cW_5E-5_1E5.hepmc2g",        "EFT_cW",       "EFT cW = 5E-5",                                30.6749 },
    { "EFT_220_cB_9E-4_1E5.hepmc2g",        "EFT_cB",       "EFT cB = 9E-4",                                31.6042 },
    { "EFT_220_all_1E5.hepmc2g",            "EFT_all",      "EFT (all)",                                    42.8231 },

    { "AGC_211_lambda_127E-3_1E5.hepmc2g",  "AGC_lambda",   "AGC #lambda_{#gamma/Z} = 0.127",               30.5842 },
    { "AGC_211_g1_121E-2_1E5.hepmc2g",      "AGC_g1",       "AGC #Deltag1_{Z} = 0.208",                     30.7988 },
    { "AGC_211_kappa_391E-2_1E5.hepmc2g",   "AGC_kappa",    "AGC #Delta#kappa_{#gamma/Z} = 2.91/-0.83",     31.6995 },
    { "AGC_211_all_1E5.hepmc2g",            "AGC_all",      "AGC (all)",                                    42.7991 },
};

static const ModelFileVector Models_1E6 =
{
    { "SM_211_1E6.hepmc2g",                 "SM_211",       "SM (2.1.1)",                                   18.5248 },
    { "SM_220_1E6.hepmc2g",                 "SM_220",       "SM (2.2.0)",                                   18.5537 },
    { "SM_AGC_211_1E6.hepmc2g",             "SM_AGC",       "SM-AGC (2.1.1)",                               18.5432 },
    { "SM_UFO_220_1E6.hepmc2g",             "SM_UFO",       "SM-UFO (2.2.0)",                               18.5476 },

    { "EFT_220_cWWW_3E-5_1E6.hepmc2g",      "EFT_cWWW",     "EFT cWWW = 3E-5",                              30.4268 },
    { "EFT_220_cW_5E-5_1E6.hepmc2g",        "EFT_cW",       "EFT cW = 5E-5",                                30.7167 },
    { "EFT_220_cB_9E-4_1E6.hepmc2g",        "EFT_cB",       "EFT cB = 9E-4",                                31.6710 },
    { "EFT_220_all_1E6.hepmc2g",            "EFT_all",      "EFT (all)",                                    42.8051 },

    { "AGC_211_lambda_127E-3_1E6.hepmc2g",  "AGC_lambda",   "AGC #lambda_{#gamma/Z} = 0.127",               30.4655 },
    { "AGC_211_g1_121E-2_1E6.hepmc2g",      "AGC_g1",       "AGC #Deltag1_{Z} = 0.208",                     30.7150 },
    { "AGC_211_kappa_391E-2_1E6.hepmc2g",   "AGC_kappa",    "AGC #Delta#kappa_{#gamma/Z} = 2.91/-0.83",     31.7246 },
    { "AGC_211_all_1E6.hepmc2g",            "AGC_all",      "AGC (all)",                                    42.8286 },
};

////////////////////////////////////////////////////////////////////////////////

static const FigureSetupVector Compare1 =
{
    //  base, compare files,    luminosity
    { { "SM_220", "SM_211" },  0.0, { kBlack, kBlue  } },
    { { "SM_211", "SM_AGC" },  0.0, { kBlue,  kGreen } },
    { { "SM_220", "SM_AGC" },  0.0, { kBlack, kGreen } },
    { { "SM_220", "SM_UFO" },  0.0, { kBlack, kRed   } },

    { { "SM_220", "EFT_cWWW"   }, 0.0, { kBlack, kRed } },
    { { "SM_211", "AGC_lambda" }, 0.0, { kBlack, kRed } },
};

static const FigureSetupVector Compare2 =
{
    { { "SM_220", "EFT_cWWW"   }, 1.0, { kBlack, kRed  } },
    { { "SM_211", "AGC_lambda" }, 1.0, { kBlack, kBlue } },
};

static const FigureSetupVector Compare3 =
{
    { { "AGC_lambda", "EFT_cWWW" }, 0.0, { kBlue,  kRed } },
};

static const FigureSetupVector Compare4 =
{
    { { "AGC_g1", "EFT_cW" }, 0.0, { kBlue,  kRed } },
};

static const FigureSetupVector Compare5 =
{
    { { "AGC_kappa", "EFT_cB" }, 0.0, { kBlue,  kRed } },
};

static const FigureSetupVector Compare6 =
{
    { { "AGC_all", "EFT_all" }, 0.0, { kBlue,  kRed } },
};

////////////////////////////////////////////////////////////////////////////////
int main()
{
  //ModelCompare::ModelCompare( "../compare/compare1.root",  Models_1E4, Observables1, Compare1 );
  //ModelCompare::ModelCompare( "../compare/compare2b.root", Models_1E4, Observables1, Compare2 );
  //ModelCompare::ModelCompare( "../compare/compare3.root" , Models_1E6, Observables1, Compare3 );
  //ModelCompare::ModelCompare( "../compare/compare4.root" , Models_1E6, Observables1, Compare4 );
  //ModelCompare::ModelCompare( "../compare/compare5.root" , Models_1E6, Observables1, Compare5 );
    ModelCompare::ModelCompare( "../compare/compare6.root" , Models_1E6, Observables1, Compare6 );

    LogMsgInfo( "Done." );
    return 0;
}
