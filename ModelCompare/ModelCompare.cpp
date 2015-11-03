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
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
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
void Observable::FillHist( TH1D & hist, double weight, const HepMC::GenVertex & signal ) const
{
    if (hist.InheritsFrom(TProfile::Class()))
    {
        double values[2] = { };
        getFunction( signal, values, 2 );

        static_cast<TProfile &>(hist).Fill( values[0], values[1], weight );
    }
    else
    {
        double value(0);
        getFunction( signal, &value, 1 );

        hist.Fill( value, weight );
    }
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
    TProfile * pHist = new TProfile( name, title, obs.nBins, obs.xMin, obs.xMax, "" );  // sdom

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
GoodBadHists HistSplitGoodBadBins( const TH1D * pSource, const TH1D * pCompare /*= nullptr*/ )
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
std::list<GoodBadHists> HistSplitGoodBadBins( const ConstTH1DVector & hists, const ConstTH1DVector & compare )
{
    std::list<GoodBadHists> result;

    auto compItr = compare.cbegin();
    for (const TH1D * pSource : hists)
    {
        const TH1D * pComp = (compItr != compare.end()) ? *compItr++ : nullptr;

        GoodBadHists goodBad = HistSplitGoodBadBins( pSource, pComp );
        result.push_back( std::move(goodBad) );
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
std::string GetChi2ResultString( const Chi2Result & res )
{
    // output probablity as %.4f, but trim trailing zeros

    Double_t prob(res.prob);

    prob = std::round(prob * 1E4) / 1E4;    // only show first 4 decimals

    return StringFormat(
                "#chi^{2}/ndf = %.4g / %i = %.4g  p-value = %.4g",
                FMT_F(res.chi2), FMT_I(res.ndf), FMT_F(res.chi2_ndf),
                FMT_F(prob) );
}

////////////////////////////////////////////////////////////////////////////////
std::string GetChi2ResultString( const Chi2Result & res1, const Chi2Result & res2 )
{
    // output probablity as %.4f, but trim trailing zeros

    Double_t prob1(res1.prob), prob2(res2.prob);

    prob1 = std::round(prob1 * 1E4) / 1E4;  // only show first 4 decimals
    prob2 = std::round(prob2 * 1E4) / 1E4;  // only show first 4 decimals

    return StringFormat(
                "#chi^{2}/ndf = %.4g[%.4g] / %i[%i] = %.4g[%.4g]  p-value = %.4g[%.4g]",
                FMT_F(res1.chi2),       FMT_F(res2.chi2),
                FMT_I(res1.ndf),        FMT_I(res2.ndf),
                FMT_F(res1.chi2_ndf),   FMT_F(res2.chi2_ndf),
                FMT_F(prob1),           FMT_F(prob2) );
}

////////////////////////////////////////////////////////////////////////////////
Chi2Result FitToHorzLineAtOne( const TH1D & hist )
{
    TF1 horz1( "horz1", "1.0" );

    Chi2Result res;
    res.chi2     = hist.Chisquare( &horz1 );    // skips bins with zero error
    res.ndf      = (Int_t)HistErrorBinCount(hist);
    res.prob     = (res.ndf > 0 ? TMath::Prob( res.chi2, res.ndf ) : 0.0);
    res.chi2_ndf = (res.ndf > 0 ? res.chi2 / res.ndf : 0.0);

    return res;
}

////////////////////////////////////////////////////////////////////////////////
std::string GetLabel_FitToHorzLineAtOne( const TH1D & hist )
{
    Chi2Result res = FitToHorzLineAtOne( hist );

    std::string label = "Fit to 1: " + GetChi2ResultString( res );
    return label;
}

////////////////////////////////////////////////////////////////////////////////
std::string GetLabel_FitToHorzLineAtOne( const TH1D & hist1, const TH1D & hist2 )
{
    Chi2Result res1 = FitToHorzLineAtOne( hist1 );
    Chi2Result res2 = FitToHorzLineAtOne( hist2 );

    std::string label = "Fit to 1: " + GetChi2ResultString( res1, res2 );
    return label;
}

////////////////////////////////////////////////////////////////////////////////
Chi2Result FitToHorzLineAtConstant( const TH1D & hist, Double_t & cValue, Double_t & cError )
{
    cValue = cError = 0;

    TF1 horz( "horz", "pol0" );
    horz.SetParameter( 0, 1.0 );

    TH1DUniquePtr pFitHist{ (TH1D *)hist.Clone() };     // clone hist as Fit is not const

    int fitStatus = pFitHist->Fit( &horz, "NQM" );      // skips bins with zero error
    if ((fitStatus < 0) || (fitStatus % 1000 != 0))     // ignore improve (M) errors
        return Chi2Result();

    Chi2Result res;
    res.chi2     = horz.GetChisquare();
    res.ndf      = horz.GetNDF();
    res.prob     = horz.GetProb();
    res.chi2_ndf = (res.ndf > 0 ? res.chi2 / res.ndf : 0.0);

    cValue = horz.GetParameter(0);
    cError = horz.GetParError(0);

    return res;
}

////////////////////////////////////////////////////////////////////////////////
std::string GetLabel_FitToHorzLineAtConstant( const TH1D & hist )
{
    Double_t cValue(0), cError(0);

    Chi2Result res = FitToHorzLineAtConstant( hist, cValue, cError );

    std::string label = StringFormat( "Fit to c = %.2g#pm%.2g: ",
                                      FMT_F(cValue), FMT_F(cError) );

    label += GetChi2ResultString( res );

    return label;
}

////////////////////////////////////////////////////////////////////////////////
std::string GetLabel_FitToHorzLineAtConstant( const TH1D & hist1, const TH1D & hist2 )
{

    Double_t cValue1(0), cError1(0);
    Double_t cValue2(0), cError2(0);

    Chi2Result res1 = FitToHorzLineAtConstant( hist1, cValue1, cError1 );
    Chi2Result res2 = FitToHorzLineAtConstant( hist2, cValue2, cError2 );

    std::string label = StringFormat( "Fit to c = %.2g#pm%.2g[%.2g#pm%.2g]: ",
                                      FMT_F(cValue1), FMT_F(cError1),
                                      FMT_F(cValue2), FMT_F(cError2) );

    label += GetChi2ResultString( res1, res2 );

    return label;
}

////////////////////////////////////////////////////////////////////////////////
void HistScaleTextTicks( TH1D & hist, Float_t vert, Float_t horz = 1 )
{
    TAxis * xAxis = hist.GetXaxis();
    TAxis * yAxis = hist.GetYaxis();

    if (vert <= 0) vert = 1;
    if (horz <= 0) horz = 1;

    if (vert != 1)
    {
        // titles, labels, and x-ticks are vertically-sized

        xAxis->SetTitleSize(  xAxis->GetTitleSize()  * vert );
        xAxis->SetLabelSize(  xAxis->GetLabelSize()  * vert );
        xAxis->SetTickLength( xAxis->GetTickLength() * vert );

        yAxis->SetTitleSize(  yAxis->GetTitleSize()  * vert );
        yAxis->SetLabelSize(  yAxis->GetLabelSize()  * vert );

        // drawn title offset scales with both title size and offset
        // correct y-axis title-offset so title does not move even though it is larger
        yAxis->SetTitleOffset( yAxis->GetTitleOffset() / vert );
    }

    if (horz != 1)
    {
        // y-ticks are horizontally-sized

        yAxis->SetTickLength( yAxis->GetTickLength() * horz );
    }
}

////////////////////////////////////////////////////////////////////////////////
void HistScaleTextTicks( const TH1DVector & hists, Float_t vert, Float_t horz = 1 )
{
    for (TH1D * pHist : hists)
    {
        if (pHist)
            HistScaleTextTicks( *pHist, vert, horz );
    }
}

////////////////////////////////////////////////////////////////////////////////
void WriteCompareFigure( const char * name, const char * title, const ConstTH1DVector & data, const ConstTH1DVector & compare, const ColorVector & dataColors,
                         const ConstTH1DVector & rawData )
{
    auto SetupCompareHists = []( const TH1DVector & hists ) -> void
    {
        for (TH1D * pHist : hists)
        {
            if (pHist)
            {
                //pHist->SetLineWidth(2);
                pHist->GetXaxis()->CenterTitle();
                pHist->GetYaxis()->CenterTitle();
            }
        }
    };

    /////

    const Double_t LowerPadFraction = 1.0/3.0;
    const Double_t UpperPadFraction = 1.0 - LowerPadFraction;

    TCanvas canvas( name, title );

    // divide the canvas into two pads
    {
        canvas.SetMargin(0, 0, 0, 0);   // clear margins before division, so that entire canvas is "owned" by subpads
        canvas.Divide(1,2,0,0);         // divide canvas into an upper and lower pad, with no space between pads

        TVirtualPad * pPad = nullptr;

        // setup upper pad
        {
            pPad = canvas.GetPad(1);

            pPad->SetPad( 0, LowerPadFraction, 1, 1 );  // xlow, ylow, xup, yup

            pPad->UseCurrentStyle();    // restore margins to default after division
            pPad->SetBottomMargin(0);   // remove bottom margin

            pPad->SetTopMargin( Float_t(pPad->GetTopMargin() / UpperPadFraction) );         // increase top margin
        }

        // setup lower pad
        {
            pPad = canvas.GetPad(2);

            pPad->SetPad( 0, 0, 1, LowerPadFraction );  // xlow, ylow, xup, yup

            pPad->UseCurrentStyle();    // restore margins to default after division
            pPad->SetTopMargin(0);      // remove top margin

            pPad->SetBottomMargin( Float_t(pPad->GetBottomMargin() / LowerPadFraction) );   // increase bottom margin
        }

        //Double_t xlow, ylow, xup, yup;
        //pPad->GetPadPar(xlow, ylow, xup, yup);
        //LogMsgInfo( "Before: x:%f->%f y:%f->%f L:%f R:%f B:%f T:%f", FMT_F(xlow), FMT_F(xup), FMT_F(ylow), FMT_F(yup),
        //            FMT_F(pPad->GetLeftMargin()), FMT_F(pPad->GetRightMargin()), FMT_F(pPad->GetBottomMargin()), FMT_F(pPad->GetTopMargin()) );
    }

    // draw upper pad
    {
        LogMsgInfo( "\n--- %hs : pad 1 ---", FMT_HS(name) );

        canvas.cd(1);

        // draw the histograms
        TH1DVector drawHists = DrawMultipleHist( title, data, dataColors );  // drawHists are owned by the current pad

        SetupCompareHists( drawHists );

        HistScaleTextTicks( drawHists, 1/UpperPadFraction );

        // determine good/bad histograms
        std::list<GoodBadHists> goodBadData = HistSplitGoodBadBins( ToConstTH1DVector(drawHists), rawData );

        // draw bad hists
        for (const auto & gb : goodBadData)
        {
            gb.bad->SetMarkerStyle( kOpenCircle );

            if (gb.bad->GetEffectiveEntries() != 0)
                gb.bad->DrawCopy( "SAME" );     // draw copy so object persists after goodBadData goes out of scope
        }

        // add a customized legend, different than TPad::BuildLegend
        {
            TLegend * pLegend = new TLegend( 0.33, 0.67, 0.88, 0.88 );  // using default position used by TPad::BuildLegend
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

                        std::string label = StringFormat( "Kolmogorov = %.3g[%.3g]", FMT_F(probGood), FMT_F(probAll) );
                        LogMsgInfo( label );
                        pLegend->AddEntry( (TObject *)nullptr, label.c_str(), "" );
                    }

                    // add Chi2Test probability
                    {
                        Chi2Result chi2All;
                        chi2All.Chi2Test( *pBaseAll, *pCompAll );       // supports both TH1D and TProfile

                        Chi2Result chi2Good;
                        chi2Good.Chi2Test( *pBaseGood, *pCompGood );    // supports both TH1D and TProfile

                        std::string label = GetChi2ResultString( chi2Good, chi2All );
                        LogMsgInfo( label );
                        pLegend->AddEntry( (TObject *)nullptr, label.c_str(), "" );
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

        SetupCompareHists( drawHists );

        HistScaleTextTicks( drawHists, 1/LowerPadFraction );

        // determine good/bad histograms
        std::list<GoodBadHists> goodBadCompare;
        {
            size_t i = 1;
            for (const TH1D * pHist : drawHists)
            {
                GoodBadHists goodBad1 = HistSplitGoodBadBins( pHist,               rawData[0]   );
                GoodBadHists goodBad2 = HistSplitGoodBadBins( goodBad1.good.get(), rawData[i++] );

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
            //TLegend * pLegend = new TLegend( 0.12, 0.53, 0.7, 0.88 );  // .2 wider than default and aligned to left
            //pLegend->SetMargin( 0.1 );  // reduce width for entry symbol from 25% to 10%

            TLegend legendStyle( 0.12, 0.74, 0.6, 0.88 );

            TPaveText * pTextBox = new TPaveText;
            *static_cast<TAttText *>( pTextBox ) = legendStyle;
            *static_cast<TPave *>(    pTextBox ) = legendStyle;
            pTextBox->SetMargin( 0.01 );

            auto gbitr = goodBadCompare.cbegin();

            for (size_t i = 0; i < compare.size(); ++i, ++gbitr)
            {
                const TH1D * pDrawHist = drawHists[i];

                const TH1D * pCompAll  = pDrawHist;
                const TH1D * pCompGood = gbitr->good.get();

                // add a fit to a horizontal line at y=1.0
                {
                    std::string label = GetLabel_FitToHorzLineAtOne( *pCompGood, *pCompAll );
                    LogMsgInfo( label );
                    pTextBox->AddText( label.c_str() );
                }

                // add a fit to a horizontal line at a y=c
                {
                    std::string label = GetLabel_FitToHorzLineAtConstant( *pCompGood, *pCompAll );
                    LogMsgInfo( label );
                    pTextBox->AddText( label.c_str() );
                }
            }

            pTextBox->SetBit( kCanDelete );  // inform pad that it can delete this object
            pTextBox->Draw();                // add legend to current pad's list of primatives
        }
    }

    // write canvas
    canvas.Write();
}

////////////////////////////////////////////////////////////////////////////////
bool LoadCacheHist( const char * cacheFileName, TH1D * & pHist )
{
    if (!cacheFileName || !cacheFileName[0] || !pHist)
        return false;

    std::unique_ptr<TH1D> pHistCache( LoadHist( cacheFileName, pHist->GetName() ) );
    if (!pHistCache)
        return false;

    if (pHistCache->Class() != pHist->Class())
        return false;

    try
    {
        struct MyTH1 : public TH1
        {
            using TH1::CheckConsistency;
        };

        if (MyTH1::CheckConsistency( pHistCache.get(), pHist ))
        {
            delete pHist;
            pHist = pHistCache.release();
            return true;
        }
    }
    catch (...) { }

    return false;
}

////////////////////////////////////////////////////////////////////////////////
void LoadHistData( const ModelFileVector & models, const ObservableVector & observables, std::vector<TH1DVector> & hists, const char * cacheFileName /*= nullptr*/ )
{
    hists.clear();

    for (const ModelFile & model : models)
    {
        bool bLoadEvents = false;

        TH1DVector data;
        TH1DVector load;

        for (const Observable & obs : observables)
        {
            TH1D * pHist = obs.MakeHist( model.modelName, model.modelTitle );

            if (LoadCacheHist( cacheFileName, pHist ))
            {
                LogMsgInfo( "Loaded %hs from cache", FMT_HS(pHist->GetName()) );
                load.push_back( nullptr );  // skip this histogram
            }
            else
            {
                load.push_back( pHist );
                bLoadEvents = true;
            }

            data.push_back( pHist );
        }

        auto FillFunc = [&](const HepMC::GenVertex & signal)
        {
            size_t obsIndex = 0;
            for (const Observable & obs : observables)
            {
                TH1D * pHist = load[obsIndex++];
                if (pHist)
                    obs.FillHist( *pHist, 1.0, signal );
            }
        };

        if (bLoadEvents)
        {
            LoadEvents( model.fileName, FillFunc, model.maxLoadEvents );

            if (cacheFileName && cacheFileName[0])
                SaveHists( cacheFileName, ToConstTH1DVector(data) );
        }
        
        hists.push_back( data );
    }
}

////////////////////////////////////////////////////////////////////////////////
void CalculateCompareHists( const Observable & obs, const ConstTH1DVector & data, TH1DVector & comp, const ModelFileVector & models, const ColorVector & dataColors )
{
    comp.clear();

    // calculate comparison histograms

    std::unique_ptr<const TH1D> upBase( ConvertTProfileToTH1D( data.front(), false ) );

    std::string nameSuffix  = "_vs_" + std::string(models[0].modelName)  + "_"   + std::string(obs.name);
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
void ModelCompare( const char * outputFileName,
                   const ModelFileVector & models, const ObservableVector & observables,
                   const FigureSetupVector & figures,
                   const char * cacheFileName /*= nullptr*/ )
{
    // disable automatic histogram addition to current directory
    TH1::AddDirectory(kFALSE);
    // enable automatic sumw2 for every histogram
    TH1::SetDefaultSumw2(kTRUE);

    // modify the global style
    gStyle->SetPaperSize( TStyle::kA4 );
    gStyle->SetTitleOffset( 1.3, "xyz" ); // increase title offsets a little more
    gStyle->SetPadTopMargin(   0.03 );
    gStyle->SetPadRightMargin( 0.03 );
    gStyle->SetPadLeftMargin(  0.09 );
    gStyle->SetOptTitle( kFALSE );

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
    LoadHistData( loadModels, observables, modelData, cacheFileName );

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

                WriteCompareFigure( figName.c_str(), figTitle.c_str(), obsData, ToConstTH1DVector(obsComp), figSetup.colors, obsData );
            }
        }
    }

    //upOutputFile->Write( 0, TFile::kOverwrite );
    upOutputFile->Close();
}

////////////////////////////////////////////////////////////////////////////////

} // namespace ModelCompare
