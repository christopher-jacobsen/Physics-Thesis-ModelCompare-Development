//
//  ModelCompare.cpp
//  ModelCompare
//
//  Created by Christopher Jacobsen on 11/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "ModelCompare.h"

#include "common.h"
#include "RootUtil.h"

// Root includes
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>

////////////////////////////////////////////////////////////////////////////////

using namespace RootUtil;

////////////////////////////////////////////////////////////////////////////////

namespace ModelCompare
{

////////////////////////////////////////////////////////////////////////////////

std::string Observable::BuildHistName(  const char * namePrefix  /*= nullptr*/, const char * nameSuffix  /*= nullptr*/ ) const
{
    std::string sName;

    if (namePrefix && namePrefix[0]) { sName += namePrefix; sName += "_"; }
    sName += this->name;
    if (nameSuffix && nameSuffix[0]) { sName += "_"; sName += nameSuffix; }

    return sName;
}

std::string Observable::BuildHistTitle( const char * titlePrefix /*= nullptr*/, const char * titleSuffix /*= nullptr*/ ) const
{
    std::string sTitle;

    if (titlePrefix && titlePrefix[0]) { sTitle += titlePrefix; sTitle += " - "; }
    sTitle += this->title;
    if (titleSuffix && titleSuffix[0]) { sTitle += " - "; sTitle += titleSuffix; }

    return sTitle;
}

////////////////////////////////////////////////////////////////////////////////

TH1D * DefaultTH1DFactory( const Observable & obs, const char * name, const char * title )
{
    TH1D * pHist = new TH1D( name, title, obs.nBins, obs.xMin, obs.xMax );

    SetupHist( *pHist, obs.xAxisTitle, obs.yAxisTitle );

    return pHist;
}

TH1D * DefaultTProfileFactory( const Observable & obs, const char * name, const char * title )
{
    TProfile * pHist = new TProfile( name, title, obs.nBins, obs.xMin, obs.xMax );

    SetupHist( *pHist, obs.xAxisTitle, obs.yAxisTitle );

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////

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
void ScaleHistToLuminosity( double luminosity, const TH1DVector & hists, const ModelFile & eventFile, bool bApplyCrossSectionError /*= false*/ )
{
    RootUtil::ScaleHistToLuminosity( luminosity, hists, eventFile.crossSectionEvents, eventFile.crossSection,
                                     eventFile.crossSectionError, bApplyCrossSectionError );
}

////////////////////////////////////////////////////////////////////////////////

struct GoodBadHists
{
    TH1DUniquePtr   good;
    TH1DUniquePtr   bad;
};

////////////////////////////////////////////////////////////////////////////////
GoodBadHists HistSplitGoodBadBins( const TH1D * pSource, const TH1D * pCompare = nullptr )
{
    const Double_t GoodStatMinEvents = 10;

    if (!pSource)
        return { nullptr, nullptr };

    if (!pCompare)
        pCompare = pSource;

    const std::string sourceName( pSource->GetName() );

    TH1D * pGood = (TH1D *)pSource->Clone( (sourceName + "_good").c_str() );    // polymorhic clone
    pGood->SetDirectory( nullptr );                                             // ensure not owned by any directory

    TH1D * pBad  = (TH1D *)pSource->Clone( (sourceName + "_bad").c_str() );     // polymorphic clone
    pBad->SetDirectory( nullptr );                                              // ensure not owned by any directory

    size_t nGood(0), nBad(0), nEmpty(0);

    const Int_t nSize = pCompare->GetSize();
    for (Int_t bin = 0; bin < nSize; ++bin)  // include under/overflow bins
    {
        Double_t compEffEntries = GetHistBinEffectiveEntries( *pCompare, bin );

        bool bGood = (compEffEntries >= GoodStatMinEvents * (1.0 - std::numeric_limits<Double_t>::epsilon()));

        const TH1D * pKeep =  bGood ? pGood : pBad;
        TH1D *       pZero = !bGood ? pGood : pBad;

        // if good zero bad and visa versa
        ZeroHistBin( *pZero, bin );

        Double_t     keepEffEntries = GetHistBinEffectiveEntries( *pKeep, bin );
        const char * sType          = nullptr;

        if (keepEffEntries == 0)
        {
            ++nEmpty;
            sType = "Empty";
        }
        else if (bGood)
        {
            ++nGood;
            sType = "Good";
        }
        else
        {
            ++nBad;
            sType = "Bad";
        }

        //LogMsgInfo( "%i Comp: %g (±%g)[%g]        %hs: %g (±%g)[%g]", FMT_I(bin),
        //    FMT_F(pCompare->GetBinContent(bin)), FMT_F(pCompare->GetBinError(bin)), FMT_F(compEffEntries), FMT_HS(sType),
        //    FMT_F(pKeep   ->GetBinContent(bin)), FMT_F(pKeep   ->GetBinError(bin)), FMT_F(keepEffEntries) );
    }

    pGood->ResetStats();
    pBad ->ResetStats();

    LogMsgInfo( "HistSplitGoodBadBins: %hs using %hs -> %u bins: %u good, %u bad, %u empty",
                FMT_HS(pSource->GetName()), FMT_HS(pCompare->GetName()),
                FMT_I(nSize), FMT_U(nGood), FMT_U(nBad), FMT_U(nEmpty) );

    return { TH1DUniquePtr(pGood), TH1DUniquePtr(pBad) };
}

////////////////////////////////////////////////////////////////////////////////
std::list<GoodBadHists> HistSplitGoodBadBins( const ConstTH1DVector & hists )
{
    std::list<GoodBadHists> result;

    for (const TH1D * pOrigHist : hists)
    {
        GoodBadHists goodBad = HistSplitGoodBadBins( pOrigHist );
        result.push_back( std::move(goodBad) );
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
std::string GetLabel_FitToHorzLineAtOne( const TH1D & hist )
{
    TF1 horz1( "horz1", "1.0" );

    Chi2Result res;
    res.chi2     = hist.Chisquare( &horz1 );
    res.ndf      = (Int_t)HistNonEmptyBinCount(hist);
    res.prob     = (res.ndf > 0 ? TMath::Prob( res.chi2, res.ndf ) : 0.0);
    res.chi2_ndf = (res.ndf > 0 ? res.chi2 / res.ndf : 0.0);

    std::string label = "Fit to 1: " + res.Label();
    return label;
}

////////////////////////////////////////////////////////////////////////////////
std::string GetLabel_FitToHorzLineAtConstant( const TH1D & hist )
{
    TF1 horz( "horz", "pol0" );
    horz.SetParameter( 0, 1.0 );

    TH1DUniquePtr pFitHist{ (TH1D *)hist.Clone() };     // clone hist as Fit is not const

    int fitStatus = pFitHist->Fit( &horz, "NQM" );
    if ((int)fitStatus < 0)
        return "";

    Chi2Result res;
    res.chi2     = horz.GetChisquare();
    res.ndf      = horz.GetNDF();
    res.prob     = horz.GetProb();
    res.chi2_ndf = (res.ndf > 0 ? res.chi2 / res.ndf : 0.0);

    Double_t p0_val = horz.GetParameter(0);
    Double_t p0_err = horz.GetParError(0);

    char prefix[50];
    sprintf( prefix, "Fit to c = %.2f#pm%.2g: ", FMT_F(p0_val), FMT_F(p0_err) );

    std::string label = prefix + res.Label();
    return label;
}

////////////////////////////////////////////////////////////////////////////////
void WriteCompareFigure( const char * name, const char * title, const ConstTH1DVector & data, const ConstTH1DVector & compare, const ColorVector & dataColors )
{
    TCanvas canvas( name, title );

    canvas.Divide(1,2); // divide canvas into an upper and lower pad
    
    // draw upper pad
    {
        LogMsgInfo( "\n--- %hs : pad 1 ---", FMT_HS(name) );

        canvas.cd(1);

        // draw the histograms
        TH1DVector drawHists = DrawMultipleHist( title, data, dataColors );  // drawHists are owned by the current pad

        // determine good/bad histograms
        std::list<GoodBadHists> goodBadData = HistSplitGoodBadBins( ToConstTH1DVector(drawHists) );

        // draw bad hists
        for (const auto & gb : goodBadData)
        {
            gb.bad->SetMarkerStyle( kOpenCircle );

            if (gb.bad->GetEffectiveEntries() != 0)
                gb.bad->DrawCopy( "SAME" );     // draw copy so object persists after goodBadData goes out of scope
        }

        // add a customized legend, different than TPad::BuildLegend
        {
            TLegend * pLegend = new TLegend( 0.4, 0.53, 0.88, 0.88 );  // using default parameters to TPad::BuildLegend

            pLegend->SetMargin( 0.1 );  // reduce width for entry symbol from 25% to 10%

            auto gbitr = goodBadData.cbegin();

            const TH1D * pBaseAll  = drawHists[0];
            const TH1D * pBaseGood = gbitr->good.get();

            for (size_t i = 0; i < drawHists.size(); ++i, ++gbitr)
            {
                const TH1D * pDrawHist = drawHists[i];

                pLegend->AddEntry( pDrawHist, pDrawHist->GetTitle() );

                if (i != 0)
                {
                    const TH1D * pCompAll  = pDrawHist;
                    const TH1D * pCompGood = gbitr->good.get();

                    // add Kolmogorov probability
                    {
                        Double_t probAll  = KolmogorovTest_NonEmptyBins( *pBaseAll,  *pCompAll  );
                        Double_t probGood = KolmogorovTest_NonEmptyBins( *pBaseGood, *pCompGood );

                        char label[100];
                        sprintf( label, "Kolmogorov = all:%0.4f good:%0.4f", FMT_F(probAll), FMT_F(probGood) );
                        LogMsgInfo( label );

                        pLegend->AddEntry( (TObject *)nullptr, label, "" );
                    }

                    // add Chi2Test probability
                    {
                        Chi2Result chi2All;
                        chi2All.Chi2Test( *pBaseAll, *pCompAll );

                        std::string labelAll = std::string("All ") + chi2All.Label();
                        LogMsgInfo( labelAll );
                        pLegend->AddEntry( (TObject *)nullptr, labelAll.c_str(), "" );

                        Chi2Result chi2Good;
                        chi2Good.Chi2Test( *pBaseGood, *pCompGood );

                        std::string labelGood = std::string("Good ") + chi2Good.Label();
                        LogMsgInfo( labelGood );
                        pLegend->AddEntry( (TObject *)nullptr, labelGood.c_str(), "" );
                    }
                }
            }

            pLegend->SetBit( kCanDelete );  // inform pad that it can delete this object
            pLegend->Draw();                // add legend to current pad's list of primatives
        }
    }

    // draw lower pad
    {
        LogMsgInfo( "\n--- %hs : pad 2 ---", FMT_HS(name) );

        canvas.cd(2);

        // draw the histograms
        TH1DVector drawHists = DrawMultipleHist( "", compare );

        // determine good/bad histograms
        std::list<GoodBadHists> goodBadCompare;
        {
            size_t i = 1;
            for (const TH1D * pHist : drawHists)
            {
                GoodBadHists goodBad1 = HistSplitGoodBadBins( pHist,               data[0]   );
                GoodBadHists goodBad2 = HistSplitGoodBadBins( goodBad1.good.get(), data[i++] );

                goodBad2.bad->Add( goodBad1.bad.get() );  // add the two bad hists together

                goodBadCompare.push_back( std::move(goodBad2) );
            }
        }

        // draw bad hists
        for (const auto & gb : goodBadCompare)
        {
            gb.bad->SetMarkerStyle( kOpenCircle );

            if (gb.bad->GetEffectiveEntries() != 0)
                gb.bad->DrawCopy( "SAME" );     // draw copy so object persists after goodBadData goes out of scope
        }

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
            TLegend * pLegend = new TLegend( 0.12, 0.53, 0.7, 0.88 );  // .1 wider than default and aligned to left

            pLegend->SetMargin( 0.1 );  // reduce width for entry symbol from 25% to 10%

            std::string label;

            auto gbitr = goodBadCompare.cbegin();

            for (size_t i = 0; i < compare.size(); ++i, ++gbitr)
            {
                const TH1D * pDrawHist = drawHists[i];

                pLegend->AddEntry( pDrawHist, pDrawHist->GetTitle() );

                const TH1D * pCompAll  = pDrawHist;
                const TH1D * pCompGood = gbitr->good.get();

                // add a fit to a horizontal line at y=1.0
                {
                    label = "All " + GetLabel_FitToHorzLineAtOne( *pCompAll );
                    LogMsgInfo( label );
                    pLegend->AddEntry( (TObject *)nullptr, label.c_str(), "" );

                    label = "Good " + GetLabel_FitToHorzLineAtOne( *pCompGood );
                    LogMsgInfo( label );
                    pLegend->AddEntry( (TObject *)nullptr, label.c_str(), "" );
                }

                // add a fit to a horizontal line at a y=c
                {
                    label = "All " + GetLabel_FitToHorzLineAtConstant( *pCompAll );
                    LogMsgInfo( label );
                    pLegend->AddEntry( (TObject *)nullptr, label.c_str(), "" );

                    label = "Good " + GetLabel_FitToHorzLineAtConstant( *pCompGood );
                    LogMsgInfo( label );
                    pLegend->AddEntry( (TObject *)nullptr, label.c_str(), "" );
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
void LoadHistData( const ModelFileVector & models, const ObservableVector & observables, std::vector<TH1DVector> & hists )
{
    hists.clear();

    for (const ModelFile & model : models)
    {
        TH1DVector data;

        for (const Observable & obs : observables)
        {
            TH1D * pHist = obs.MakeHist( model.modelName, model.modelTitle );
            data.push_back( pHist );
        }

        auto FillFunc = [&](const HepMC::GenVertex & signal)
        {
            size_t dataIndex = 0;
            for (const Observable & obs : observables)
            {
                obs.fillFunction( *data[dataIndex++], 1.0, signal );
            }
        };

        LoadEvents( model.fileName, FillFunc );

        hists.push_back( data );
    }
}

////////////////////////////////////////////////////////////////////////////////
void CalculateCompareHists( const Observable & obs, const ConstTH1DVector & data, TH1DVector & comp, const ModelFileVector & models, const ColorVector & dataColors )
{
    comp.clear();

    // calculate comparison histograms

    std::unique_ptr<const TH1D> upBase( ConvertTProfileToTH1D( data.front(), false ) );

    std::string nameSuffix  = "_vs_" + std::string(models[0].modelName) + "_" + std::string(obs.name);
    std::string titleSuffix = " vs " + std::string(models[0].modelTitle) + " - " + obs.title;

    for ( size_t i = 1; i < data.size(); ++i)
    {
        TH1D * pHist = ConvertTProfileToTH1D( data[i], false );

        pHist->Divide( upBase.get() );

        comp.push_back(pHist);

        /* add if comparing to oneself - to correct the error for being increased by sqrt(2)
        if (i == 0)
        {
            Int_t nbins = pHist->GetNbinsX();
            for (Int_t bin = 0; bin <= nbins + 1; ++bin)
            {
                Double_t error = pHist->GetBinError(bin);
                error /= std::sqrt(2);
                pHist->SetBinError( bin, error );
                pHist->SetBinContent( bin, 1.0 );
            }
        }
        */

        std::string name  = std::string(models[i].modelName)  + nameSuffix;
        std::string title = std::string(models[i].modelTitle) + titleSuffix;

        pHist->SetName(  name .c_str() );
        pHist->SetTitle( title.c_str() );

        Color_t color = dataColors[i];
        pHist->SetLineColor(   color );
        pHist->SetMarkerColor( color );

        pHist->GetYaxis()->SetTitle( "Ratio" );
    }
}

////////////////////////////////////////////////////////////////////////////////
void ModelCompare( const char * outputFileName, const ModelFileVector & models, const ObservableVector & observables, const FigureSetupVector & figures )
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

    std::vector< TH1DUniquePtr > tempHists;

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
        LogMsgHistUnderOverflow( ToConstTH1DVector(data) );
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

        // adjust for luminosity
        if (figSetup.luminosity > 0)
        {
            double luminosity = figSetup.luminosity;  // fb^-1

            for (size_t modelIndex = 0; modelIndex < figModels.size(); ++modelIndex)
            {
                double crossSection = figModels[modelIndex].crossSection * 1000; // fb
                double nEntries     = figData[modelIndex][0]->GetEntries();
                size_t nEvents      = (size_t)nEntries;

                if (((double)nEvents != nEntries) || !nEvents)
                    ThrowError( "Non-integral number of entries: " + std::to_string(nEntries) );

                double scale = luminosity * crossSection / nEvents;

                TH1DVector & obsData = figData[modelIndex];
                for ( TH1D * & pHist : obsData )
                {
                    if (pHist->GetEntries() != nEntries)
                        ThrowError( "Inconsistent number of entries: " + std::to_string(pHist->GetEntries()) + " expected: " + std::to_string(nEntries) );

                    TH1D * pClone = dynamic_cast<TH1D *>( pHist->Clone() );

                    pClone->Scale( scale );
                    pClone->SetDirectory( nullptr );    // ensure not owned by any directory

                    pHist = pClone;  // replace histogram in obsData

                    tempHists.push_back( TH1DUniquePtr(pClone) );
                }
            }
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

                WriteCompareFigure( figName.c_str(), figTitle.c_str(), obsData, ToConstTH1DVector(obsComp), figSetup.colors );
            }
        }
    }

    //upOutputFile->Write( 0, TFile::kOverwrite );
    upOutputFile->Close();
}

////////////////////////////////////////////////////////////////////////////////

} // namespace ModelCompare
