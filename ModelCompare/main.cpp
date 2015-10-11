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
    // phase-space observables

    { "PTZ", "P_{T}(Z)",  150,     0,  750, "P_{T}(Z) [GeV/c]", "Events per 5 GeV/c",       GETOBS{ GetObs(s,v,c, GetObsPT,   24);     } },
    { "MWZ", "M(WZ)",     150,     0, 3000, "M(WZ) [GeV/c^2]",  "Events per 20 GeV/c^2",    GETOBS{ GetObs(s,v,c, GetObsMass, 24, 23); } },
    { "RAZ", "Y(Z)",      100,    -5,    5, "Y(Z)",             "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsRap,  24);     } },
    { "ETZ", "#eta(Z)",   100,   -10,   10, "#eta(Z)",          "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsEta,  24);     } },
    { "PHZ", "#phi(Z)",   100, -M_PI, M_PI, "#phi(Z)",          "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsPhi,  24);     } },
};

////////////////////////////////////////////////////////////////////////////////

static const ModelFileVector Models_1E4 =
{
    { "SM_211_1E4.hepmc2g",                 "SM_211",       "SM (2.1.1)",                                   18.3748, 0.151908, 10000 },
    { "SM_220_1E4.hepmc2g",                 "SM_220",       "SM (2.2.0)",                                   18.2613, 0.154079, 10000 },
    { "SM_AGC_211_1E4.hepmc2g",             "SM_AGC",       "SM-AGC (2.1.1)",                               18.6957, 0.154950, 10000 },
    { "SM_UFO_220_1E4.hepmc2g",             "SM_UFO",       "SM-UFO (2.2.0)",                               18.2796, 0.154169, 10000 },

    { "EFT_220_cWWW_3E-5_1E4.hepmc2g",      "EFT_cWWW",     "EFT cWWW = 3E-5",                              30.3160, 0.266402, 10000 },
    { "EFT_220_cW_5E-5_1E4.hepmc2g",        "EFT_cW",       "EFT cW = 5E-5",                                30.7824, 0.272221, 10000 },
    { "EFT_220_cB_9E-4_1E4.hepmc2g",        "EFT_cB",       "EFT cB = 9E-4",                                31.4724, 0.260784, 10000 },
    { "EFT_220_all_1E4.hepmc2g",            "EFT_all",      "EFT (all)",                                    42.9091, 0.379962, 10000 },

    { "AGC_211_lambda_127E-3_1E4.hepmc2g",  "AGC_lambda",   "AGC #lambda_{#gamma/Z} = 0.127",               30.5172, 0.267667, 10000 },
    { "AGC_211_g1_121E-2_1E4.hepmc2g",      "AGC_g1",       "AGC #Deltag1_{Z} = 0.208",                     30.6187, 0.270293, 10000 },
    { "AGC_211_kappa_391E-2_1E4.hepmc2g",   "AGC_kappa",    "AGC #Delta#kappa_{#gamma/Z} = 2.91/-0.83",     31.5387, 0.260345, 10000 },
    { "AGC_211_all_1E4.hepmc2g",            "AGC_all",      "AGC (all)",                                    43.1840, 0.386040, 10000 },
};

static const ModelFileVector Models_1E5 =
{
    { "SM_211_1E5.hepmc2g",                 "SM_211",       "SM (2.1.1)",                                   18.4996, 0.0482874, 100000 },
    { "SM_220_1E5.hepmc2g",                 "SM_220",       "SM (2.2.0)",                                   18.4850, 0.0491951, 100000 },
    { "SM_AGC_211_1E5.hepmc2g",             "SM_AGC",       "SM-AGC (2.1.1)",                               18.5791, 0.0487636, 100000 },
    { "SM_UFO_220_1E5.hepmc2g",             "SM_UFO",       "SM-UFO (2.2.0)",                               18.4768, 0.0491667, 100000 },

    { "EFT_220_cWWW_3E-5_1E5.hepmc2g",      "EFT_cWWW",     "EFT cWWW = 3E-5",                              30.4574, 0.0845763, 100000 },
    { "EFT_220_cW_5E-5_1E5.hepmc2g",        "EFT_cW",       "EFT cW = 5E-5",                                30.6749, 0.0858233, 100000 },
    { "EFT_220_cB_9E-4_1E5.hepmc2g",        "EFT_cB",       "EFT cB = 9E-4",                                31.6042, 0.0827277, 100000 },
    { "EFT_220_all_1E5.hepmc2g",            "EFT_all",      "EFT (all)",                                    42.8231, 0.1199580, 100000 },

    { "AGC_211_lambda_127E-3_1E5.hepmc2g",  "AGC_lambda",   "AGC #lambda_{#gamma/Z} = 0.127",               30.5842, 0.0847977, 100000 },
    { "AGC_211_g1_121E-2_1E5.hepmc2g",      "AGC_g1",       "AGC #Deltag1_{Z} = 0.208",                     30.7988, 0.0859008, 100000 },
    { "AGC_211_kappa_391E-2_1E5.hepmc2g",   "AGC_kappa",    "AGC #Delta#kappa_{#gamma/Z} = 2.91/-0.83",     31.6995, 0.0826413, 100000 },
    { "AGC_211_all_1E5.hepmc2g",            "AGC_all",      "AGC (all)",                                    42.7991, 0.1211260, 100000 },
};

static const ModelFileVector Models_1E6 =
{
    { "SM_211_1E6.hepmc2g",                 "SM_211",       "SM (2.1.1)",                                   18.5248, 0.0152862, 1000000 },
    { "SM_220_1E6.hepmc2g",                 "SM_220",       "SM (2.2.0)",                                   18.5537, 0.0156025, 1000000 },
    { "SM_AGC_211_1E6.hepmc2g",             "SM_AGC",       "SM-AGC (2.1.1)",                               18.5432, 0.0153974, 1000000 },
    { "SM_UFO_220_1E6.hepmc2g",             "SM_UFO",       "SM-UFO (2.2.0)",                               18.5476, 0.0155949, 1000000 },

    { "EFT_220_cWWW_3E-5_1E6.hepmc2g",      "EFT_cWWW",     "EFT cWWW = 3E-5",                              30.4268, 0.0267228, 1000000 },
    { "EFT_220_cW_5E-5_1E6.hepmc2g",        "EFT_cW",       "EFT cW = 5E-5",                                30.7167, 0.0271720, 1000000 },
    { "EFT_220_cB_9E-4_1E6.hepmc2g",        "EFT_cB",       "EFT cB = 9E-4",                                31.6710, 0.0262027, 1000000 },
    { "EFT_220_all_1E6.hepmc2g",            "EFT_all",      "EFT (all)",                                    42.8051, 0.0379205, 1000000 },

    { "AGC_211_lambda_127E-3_1E6.hepmc2g",  "AGC_lambda",   "AGC #lambda_{#gamma/Z} = 0.127",               30.4655, 0.0267271, 1000000 },
    { "AGC_211_g1_121E-2_1E6.hepmc2g",      "AGC_g1",       "AGC #Deltag1_{Z} = 0.208",                     30.7150, 0.0271012, 1000000 },
    { "AGC_211_kappa_391E-2_1E6.hepmc2g",   "AGC_kappa",    "AGC #Delta#kappa_{#gamma/Z} = 2.91/-0.83",     31.7246, 0.0261498, 1000000 },
    { "AGC_211_all_1E6.hepmc2g",            "AGC_all",      "AGC (all)",                                    42.8286, 0.0383268, 1000000 },
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
