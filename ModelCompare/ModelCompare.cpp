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
void WriteCompareFigure( const char * name, const char * title, const ConstTH1DVector & data, const ConstTH1DVector & compare, const ColorVector & dataColors )
{
    // ensure local hists are deleted
    std::list< std::unique_ptr<TH1D> > cleanupHists;

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
                        Double_t prob     = data[0]->Chi2TestX( pDataHist, chi2, ndf, igood, "WW" );
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
                const TH1D * pCompHist = compare[i];

                pLegend->AddEntry( pDrawHist, pDrawHist->GetTitle() );

                // count non-empty bins for NDF used below
                Int_t nBins     = pCompHist->GetNbinsX();
                Int_t nBinsFull = 0;
                for (Int_t i = 1; i <= nBins; ++i)
                {
                    if (pCompHist->GetBinContent(i) != 0.0)
                        ++nBinsFull;
                }

                LogMsgInfo( "\n--- %hs : pad 2 ---", FMT_HS(name) );

                // add a fit to a horizontal line at y=1.0
                {
                    TF1 horz1( "horz1", "1.0" );

                    Double_t chi2     = pCompHist->Chisquare( &horz1 );
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

                    TH1D * pFitHist = (TH1D *)pCompHist->Clone();               // clone hist as Fit is not const
                    pFitHist->SetDirectory( nullptr );                          // ensure not owned by any directory
                    cleanupHists.push_back( std::unique_ptr<TH1D>(pFitHist) );  // ensure cleaned up

                    int fitStatus = pFitHist->Fit( &horz, "NQM" );
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

    const TH1D * pBase = data.front();

    std::string nameSuffix  = "_vs_" + std::string(models[0].modelName) + "_" + std::string(obs.name);
    std::string titleSuffix = " vs " + std::string(models[0].modelTitle) + " - " + obs.title;

    for ( size_t i = 1; i < data.size(); ++i)
    {
        TH1D * pHist = (TH1D *)data[i]->Clone();  // polymorphic clone
        pHist->SetDirectory( nullptr );           // ensure not owned by any directory

        pHist->Divide( pBase );

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

        /* debug - looking at small bins
        {
            Int_t nbins = pHist->GetNbinsX();
            for (Int_t bin = 1; bin <= nbins; ++bin)
            {
                Double_t binContent = pHist->GetBinContent(bin);
                if ((binContent > 0) && (binContent < 0.1))
                {
                    LogMsgInfo( "%hs bin %i: %g (±%g) / %g (±%g) = %g (±%g)",
                                FMT_HS(pHist->GetName()), FMT_I(bin),
                                FMT_F(data[i]->GetBinContent(bin)), FMT_F(data[i]->GetBinError(bin)),
                                FMT_F(pBase  ->GetBinContent(bin)), FMT_F(pBase  ->GetBinError(bin)),
                                FMT_F(pHist  ->GetBinContent(bin)), FMT_F(pHist  ->GetBinError(bin)) );
                }
            }
        }
		*/
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

    std::vector< std::unique_ptr<TH1D> > tempHists;

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
        LogMsgUnderOverflow( ToConstTH1DVector(data) );
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

                    tempHists.push_back( std::unique_ptr<TH1D>(pClone) );
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
