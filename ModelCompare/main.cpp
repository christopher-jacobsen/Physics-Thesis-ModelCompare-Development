//
//  main.cpp
//  ModelCompare
//
//  Created by Christopher Jacobsen on 31/08/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "common.h"

// Root includes
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>

#include "RootUtil.h"

using namespace RootUtil;

////////////////////////////////////////////////////////////////////////////////

typedef std::function<void(const HepMC::GenVertex & signal, TH1D & hist)> FillHistFunction;

struct Observable
{
    const char *        name;
    const char *        title;
    Int_t               nBins;
    Double_t            min;
    Double_t            max;
    const char *        xAxisTitle;
    const char *        yAxisTitle;
    FillHistFunction    fillFunction;
};

typedef std::vector<Observable> ObservableVector;

////////////////////////////////////////////////////////////////////////////////

struct ModelFile
{
    const char *    fileName;
    const char *    modelName;
    const char *    modelTitle;
    double          crossSection;  // in pb
};

typedef std::vector<ModelFile>  ModelFileVector;

////////////////////////////////////////////////////////////////////////////////

struct FigureSetup
{
  //const char *    name;   // TODO
  //const char *    title;  // TODO

    CStringVector   modelNames;
    double          luminosity  = 0.0;     // in fb^-1
    ColorVector     colors      = DefaultColors;

    FigureSetup() = default;
    FigureSetup( const CStringVector & n )                                  : modelNames(n)                             {}
    FigureSetup( const CStringVector & n, double l )                        : modelNames(n), luminosity(l)              {}
    FigureSetup( const CStringVector & n, double l, const ColorVector & c ) : modelNames(n), luminosity(l), colors(c)   {}

    static const ColorVector DefaultColors;
};

typedef std::vector<FigureSetup> FigureSetupVector;

const ColorVector FigureSetup::DefaultColors =
{
    kBlack,
    kBlue,   kGreen,  kRed,
    kViolet, kOrange, kMagenta

    // do not use:
    // kYellow - close to white
    // kGreen  - too light, bleeds into surrounding colors
    // kCyan   - too light, bleeds into surrounding colors
    // kSpring - too bright to be useful
    // kPink   - too close to kRed
    // kAzure
    // kTeal

    // TODO: try TColor::GetColorDark(kGreen) [need to call TColor::InitializeColors() first, or construct a TColor object]
};

////////////////////////////////////////////////////////////////////////////////

static void ModelCompare( const char * outputFileName, const ModelFileVector & models, const ObservableVector & observables, const FigureSetupVector & figures );

////////////////////////////////////////////////////////////////////////////////

static const ObservableVector Observables =
{
    { "PTZ", "P_{T}(Z)",  100,     0,  500, "P_{T}(Z) [GeV/c]", "Events per 5 GeV/c",       [](const HepMC::GenVertex & s, TH1D & h) { FillHistPT( s, h, 24);     } },
    { "MWZ", "M(WZ)",     100,     0, 1500, "M(WZ) [GeV/c^2]",  "Events per 15 GeV/c^2",    [](const HepMC::GenVertex & s, TH1D & h) { FillHistM2( s, h, 24, 23); } },
    { "ETZ", "#eta(Z)",   100,   -10,   10, "#eta(Z)",          "Events per bin",           [](const HepMC::GenVertex & s, TH1D & h) { FillHistEta(s, h, 24);     } },
    { "PHZ", "#phi(Z)",   100, -M_PI, M_PI, "#phi(Z)",          "Events per bin",           [](const HepMC::GenVertex & s, TH1D & h) { FillHistPhi(s, h, 24);     } },
};

////////////////////////////////////////////////////////////////////////////////

static const ModelFileVector Models_1E4 =
{
    { "SM_211_1E4.hepmc2g",             "SM_211",       "SM (2.1.1)",                   18.3748 },
    { "SM_220_1E4.hepmc2g",             "SM_220",       "SM (2.2.0)",                   18.2613 },
    { "SM_AGC_211_1E4.hepmc2g",         "SM_AGC",       "SM-AGC (2.1.1)",               18.6957 },
    { "SM_UFO_220_1E4.hepmc2g",         "SM_UFO",       "SM-UFO (2.2.0)",               18.2796 },
    { "EFT_cWWW_3E-5_220_1E4.hepmc2g",  "EFT_cWWW",     "EFT cWWW = 3E-5 (2.2.0)",      31.9028 },
};

static const ModelFileVector Models_1E5 =
{
    { "SM_211_1E5.hepmc2g",             "SM_211",       "SM (2.1.1)",                   18.4996 },
    { "SM_220_1E5.hepmc2g",             "SM_220",       "SM (2.2.0)",                   18.4850 },
    { "SM_AGC_211_1E5.hepmc2g",         "SM_AGC",       "SM-AGC (2.1.1)",               18.5791 },
    { "SM_UFO_220_1E5.hepmc2g",         "SM_UFO",       "SM-UFO (2.2.0)",               18.4768 },
};

static const ModelFileVector Models_1E6 =
{
    { "SM_211_1E6.hepmc2g",             "SM_211",       "SM (2.1.1)",                   18.5248 },
    { "SM_220_1E6.hepmc2g",             "SM_220",       "SM (2.2.0)",                   18.5537 },
    { "SM_AGC_211_1E6.hepmc2g",         "SM_AGC",       "SM-AGC (2.1.1)",               18.5432 },
    { "SM_UFO_220_1E6.hepmc2g",         "SM_UFO",       "SM-UFO (2.2.0)",               18.5476 },
};

////////////////////////////////////////////////////////////////////////////////

static const FigureSetupVector Compare1 =
{
    //  base, compare files,    luminosity
    { { "SM_220", "SM_211" },  0.0, { kBlack, kBlue  } },
    { { "SM_211", "SM_AGC" },  0.0, { kBlue,  kGreen } },
    { { "SM_220", "SM_AGC" },  0.0, { kBlack, kGreen } },
    { { "SM_220", "SM_UFO" },  0.0, { kBlack, kRed   } },

  //{ {{ "SM_220", kBlack }, { "EFT_cWWW", kRed }}, 0.0   },
  //{ {{ "SM_220", kBlack }, { "EFT_cWWW", kRed }}, 100.0 },
};

////////////////////////////////////////////////////////////////////////////////
int main()
{
    ModelCompare( "compare/compare1.root", Models_1E4, Observables, Compare1 );

    LogMsgInfo( "Done." );
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
static void WriteCompareFigure( const char * name, const char * title, const ConstTH1DVector & data, ConstTH1DVector & compare, const ColorVector & dataColors )
{
    TCanvas canvas( name, title );

    canvas.Divide(1,2); // divide canvas into an upper and lower pad
    
    // draw upper pad
    {
        canvas.cd(1);

        TH1DVector drawHists = DrawMultipleHist( title, data, dataColors );

        // add a customized legend, different than TPad::BuildLegend
        {
            TLegend * pLegend = new TLegend( 0.5, 0.67, 0.88, 0.88 );  // using default parameters to TPad::BuildLegend

            pLegend->SetMargin( 0.1 );  // reduce width for entry symbol from 25% to 10%

            for (size_t i = 0; i < drawHists.size(); ++i)
            {
                const TH1D * pDrawHist = drawHists[i];
                const TH1D * pDataHist = data[i];

                pLegend->AddEntry( pDrawHist, pDrawHist->GetTitle() );

                if (i != 0)
                {
                    LogMsgInfo( "\n--- %hs : pad 1 ---", FMT_HS(name) );

                    // add Kolmogorov probability
                    {
                        Double_t prob = data[0]->KolmogorovTest( pDataHist );

                        char label[100];
                        sprintf( label, "Kolmogorov = %0.4f", FMT_F(prob) );
                        LogMsgInfo( label );

                        pLegend->AddEntry( (TObject *)nullptr, label, "" );
                    }

                    // add Chi2Test probability
                    {
                        Double_t chi2     = 0;
                        Int_t    ndf      = 0;
                        Int_t    igood    = 0;
                        Double_t prob     = data[0]->Chi2TestX( pDataHist, chi2, ndf, igood );
                        Double_t chi2_ndf = (ndf > 0 ? chi2 / ndf : 0.0);

                        char label[200];
                        sprintf( label, "#chi^{2}/ndf = %.4g/%i = %.3f   p-value = %0.4f", FMT_F(chi2), FMT_I(ndf), FMT_F(chi2_ndf), FMT_F(prob) );
                        LogMsgInfo( label );

                        pLegend->AddEntry( (TObject *)nullptr, label, "" );
                    }
                }
            }

            pLegend->SetBit( kCanDelete );  // inform pad that it can delete this object
            pLegend->Draw();                // add legend to current pad's list of primatives
        }
    }

    // draw lower pad
    {
        canvas.cd(2);

        // draw the histograms
        TH1DVector drawHists = DrawMultipleHist( "", compare );

        // add ticks to top and right
        canvas.GetPad(2)->SetTickx();
        canvas.GetPad(2)->SetTicky();

        // TODO: possibly resize vertical min/max to exclude error bars

        // draw black horizontal line at 1
        {
            TLine * pLine = new TLine( compare[0]->GetXaxis()->GetXmin(), 1.0, compare[0]->GetXaxis()->GetXmax(), 1.0 );
            pLine->SetLineColor(kBlack);
            pLine->SetLineWidth(1);

            pLine->SetBit(kCanDelete);  // inform pad that it can delete this object
            pLine->Draw( "" );          // add line to current pad's list of primatives
        }

        // add a customized legend, different than TPad::BuildLegend
        {
            TLegend * pLegend = new TLegend( 0.12, 0.67, 0.6, 0.88 );  // .1 wider than default and aligned to left

            pLegend->SetMargin( 0.1 );  // reduce width for entry symbol from 25% to 10%

            for (size_t i = 0; i < compare.size(); ++i)
            {
                const TH1D * pDrawHist = drawHists[i];

                pLegend->AddEntry( pDrawHist, pDrawHist->GetTitle() );

                TH1D histLocal( *compare[i] );
                histLocal.SetDirectory( nullptr );  // ensure not owned by any directory

                // count non-empty bins for NDF used below
                Int_t nBins     = histLocal.GetNbinsX();
                Int_t nBinsFull = 0;
                for (Int_t i = 1; i <= nBins; ++i)
                {
                    if (histLocal.GetBinContent(i) != 0.0)
                        ++nBinsFull;
                }

                LogMsgInfo( "\n--- %hs : pad 2 ---", FMT_HS(name) );

                // add a fit to a horizontal line at y=1.0
                {
                    TF1 horz1( "horz1", "1.0" );

                    Double_t chi2     = histLocal.Chisquare( &horz1 );
                    Int_t    ndf      = nBinsFull;
                    Double_t prob     = (ndf > 0 ? TMath::Prob( chi2, ndf ) : 0.0);
                    Double_t chi2_ndf = (ndf > 0 ? chi2 / ndf : 0.0);

                    char label[200];
                    sprintf( label, "Fit to 1: #chi^{2}/ndf = %.4g/%i = %.3f   p-value = %0.4f", FMT_F(chi2), FMT_I(ndf), FMT_F(chi2_ndf), FMT_F(prob) );
                    LogMsgInfo( label );

                    pLegend->AddEntry( (TObject *)nullptr, label, "" );
                }

                // add a fit to a horizontal line at a y=c
                {
                    TF1 horz( "horz", "pol0" );

                    horz.SetParameter( 0, 1.0 );

                    int fitStatus = histLocal.Fit( &horz, "NQM" );
                    if ((int)fitStatus >= 0)
                    {
                        Double_t chi2     = horz.GetChisquare();
                        Int_t    ndf      = horz.GetNDF();
                        Double_t prob     = horz.GetProb();
                        Double_t chi2_ndf = (ndf > 0 ? chi2 / ndf : 0.0);

                        Double_t p0_val = horz.GetParameter(0);
                        Double_t p0_err = horz.GetParError(0);

                        char label[200];
                        sprintf( label, "Fit to c = %.2f#pm%.2g: #chi^{2}/ndf = %.4g/%i = %.3f   p-value = %0.4f",
                                 FMT_F(p0_val), FMT_F(p0_err), FMT_F(chi2), FMT_I(ndf), FMT_F(chi2_ndf), FMT_F(prob) );
                        LogMsgInfo( label );

                        pLegend->AddEntry( (TObject *)nullptr, label, "" );
                    }
                }
            }

            pLegend->SetBit( kCanDelete );  // inform pad that it can delete this object
            pLegend->Draw();                // add legend to current pad's list of primatives
        }
    }

    // write canvas
    canvas.Write();
}

////////////////////////////////////////////////////////////////////////////////
static void LoadHistData( const ModelFileVector & models, const ObservableVector & observables, std::vector<TH1DVector> & hists )
{
    hists.clear();

    for (const ModelFile & model : models)
    {
        TH1DVector data;

        for (const Observable & obs : observables)
        {
            TH1D * pHist = new TH1D( (std::string(model.modelName) + "_" + obs.name).c_str(),
                                     (std::string(model.modelTitle) + " - " + obs.title).c_str(),
                                     obs.nBins, obs.min, obs.max );

            SetupHist( *pHist, obs.xAxisTitle, obs.yAxisTitle );

            data.push_back( pHist );
        }

        auto FillFunc = [&](const HepMC::GenVertex & signal)
        {
            size_t dataIndex = 0;
            for (const Observable & obs : observables)
            {
                obs.fillFunction( signal, *data[dataIndex++] );
            }
        };

        LoadEvents( model.fileName, FillFunc );

        for (const TH1D * pHist : data)
            LogMsgUnderOverflow( *pHist );

        hists.push_back( data );
    }
}

////////////////////////////////////////////////////////////////////////////////
static void CalculateCompareHists( const Observable & obs, const ConstTH1DVector & data, TH1DVector & comp, const ModelFileVector & models, const ColorVector & dataColors )
{
    comp.clear();

    // calculate comparison histograms

    const TH1D * pBase = data.front();

    std::string nameSuffix  = "_vs_" + std::string(models[0].modelName) + "_" + std::string(obs.name);
    std::string titleSuffix = " vs " + std::string(models[0].modelTitle) + " - " + obs.title;

    for ( size_t i = 1; i < data.size(); ++i)
    {
        TH1D * pHist = new TH1D( *data[i] / *pBase );
        comp.push_back(pHist);

        if (i == 0)
        {
            Int_t nbins = pHist->GetNbinsX();
            for (Int_t i = 0; i <= nbins + 1; ++i)
            {
                Double_t error = pHist->GetBinError(i);
                error /= std::sqrt(2);
                pHist->SetBinError( i, error );
                pHist->SetBinContent( i, 1.0 );
            }
        }

        std::string name  = std::string(models[i].modelName)  + nameSuffix;
        std::string title = std::string(models[i].modelTitle) + titleSuffix;

        pHist->SetName(  name .c_str() );
        pHist->SetTitle( title.c_str() );

        Color_t color = dataColors[i];
        pHist->SetLineColor(   color );
        pHist->SetMarkerColor( color );

        pHist->GetYaxis()->SetTitle( "Relative frequency" );
    }
}

////////////////////////////////////////////////////////////////////////////////
static void ModelCompare( const char * outputFileName, const ModelFileVector & models, const ObservableVector & observables, const FigureSetupVector & figures )
{
    // disable automatic histogram addition to current directory
    TH1::AddDirectory(kFALSE);

    LogMsgInfo( "Output file: %hs", FMT_HS(outputFileName) );
    std::unique_ptr<TFile> upOutputFile( new TFile( outputFileName, "RECREATE" ) );
    if (upOutputFile->IsZombie() || !upOutputFile->IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create output file (%hs).", FMT_HS(outputFileName) );
        ThrowError( std::invalid_argument( outputFileName ) );
    }

    // determine which model files are to be loaded

    ModelFileVector loadModels;     // loadModels[model]
    {
        std::set<std::string> loadNames;
        for ( const FigureSetup & figSetup : figures )
            loadNames.insert( figSetup.modelNames.cbegin(), figSetup.modelNames.cend() );

        for ( const std::string & name : loadNames )
        {
            auto MatchModelName = [&name](const ModelFile & elem) -> bool { return strcmp(elem.modelName, name.c_str()) == 0; };
            
            auto itr = std::find_if( models.cbegin(), models.cend(), MatchModelName );
            if (itr == models.cend())
                ThrowError( std::invalid_argument("Model " + name + " not found.") );

            loadModels.push_back( *itr );
        }
    }

    std::vector<TH1DVector> modelData;  // modelData[model][observable]

    // load the model data for each model and observable
    LoadHistData( loadModels, observables, modelData );

    // write observables histograms
    for ( const TH1DVector & data : modelData )
    {
        WriteHists( upOutputFile.get(), data );  // output file takes ownership of histograms
    }

    // for each figure

    for ( const FigureSetup & figSetup : figures )
    {
        // select figure models and data

        ModelFileVector         figModels;  // figModels[model]
        std::vector<TH1DVector> figData;    // figData[model][observable]

        for ( const char * modelName : figSetup.modelNames )
        {
            size_t modelIndex = 0;
            for ( ; modelIndex < loadModels.size(); ++modelIndex )
            {
                if (strcmp( modelName, loadModels[modelIndex].modelName ) == 0)
                    break;
            }
            if (modelIndex == loadModels.size())
                ThrowError( std::logic_error("Internal Error: Required model not loaded.") );

            figModels.push_back( loadModels[ modelIndex ] );
            figData  .push_back( modelData[  modelIndex ] );
        }

        // for each observable

        for (size_t obsIndex = 0; obsIndex < observables.size(); ++obsIndex)
        {
            const Observable & obs = observables[obsIndex];

            ConstTH1DVector obsData;  // obsData[model]
            TH1DVector      obsComp;

            for (size_t modelIndex = 0; modelIndex < figModels.size(); ++modelIndex)
                obsData.push_back( figData[modelIndex][obsIndex] );

            // calculate the comparisons
            CalculateCompareHists( obs, obsData, obsComp, figModels, figSetup.colors );

            // write the comparison hist
            WriteHists( upOutputFile.get(), obsComp );  // output file takes ownership of histograms

            {
                std::string figName  = "fig_" + std::string(obsComp[0]->GetName());
                std::string figTitle = obsComp[0]->GetTitle();

                WriteCompareFigure( figName.c_str(), figTitle.c_str(), obsData, reinterpret_cast<ConstTH1DVector &>(obsComp), figSetup.colors );
            }
        }
    }

    //upOutputFile->Write( 0, TFile::kOverwrite );
    upOutputFile->Close();
}