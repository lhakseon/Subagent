#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TFile.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TPaveStats.h>
#include <TPaletteAxis.h>
#include <TGraphErrors.h>

#include "./SW_func1_EM.h"
#include "./SW_func2_EM.h"

void draw_median_fit(const std::string& deta_dphi, int nth_sw, const std::string& eta_region);

void Draw_SWs(){


  for(int i = 0; i < 6; i++){
     draw_median_fit( "phi", i+1, "1");
  }

  //for(int i = 0; i < 4; i++){
  //   draw_median_fit( "eta", i+1, "1");
  //}
}


void draw_median_fit( const std::string& deta_dphi, int nth_sw, const std::string& eta_region){
    
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    c1->SetRightMargin(0.1);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.15);
    c1->SetGridy();
    c1->SetGridx();
    c1->SetLogz();
    c1->cd();

    TString median = "EM_d";

    TString nth_sw_; nth_sw_.Form("%d",nth_sw);
    median = median + deta_dphi + "_SW_" + nth_sw_ + "_eta_region_" + eta_region + "_median";

    TString dist = "EMD";
    dist = dist + deta_dphi + "_" + nth_sw_ + "_eta_region_" + eta_region;
   
    c1->Update();

    TF1 *median_up = new TF1("median_up","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) +[5])", 10., 100.);
    TF1 *median_down = new TF1("median_down","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) +[5])", 10., 100.);

    for( int i=0; i < 5; i++){
        median_up->  SetParameter(i,SW_func2_dphi(i));
        median_down->SetParameter(i,SW_func2_dphi(i));
    }

    median_up->SetParameter(5, 0.02);
    median_up->GetYaxis()->SetDecimals(3);
    median_up->GetXaxis()->SetTitle("L1 EG E_{T} (crystal) [GeV]");
    median_up->GetXaxis()->SetTitleSize(0.05);
    median_up->GetYaxis()->SetTitle("#Delta#phi [rad.]");
    median_up->GetYaxis()->SetTitleSize(0.05);
    median_up->GetYaxis()->SetTitleOffset(1.35);
    median_up->GetYaxis()->SetRangeUser(-0.05, 0.09);
    median_up->Draw("l");

    median_down->SetParameter(5, -0.02);
    median_down->Draw("lsame");

    
    TF1 *sw_up = new TF1("sw_up","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4])) + ( ( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4])) -( [5]*pow(x,0) + [6]*pow(x,[7])*exp(-pow(x,[8]) +[9])) ) / 2.", 10., 100.);
   
    int phi_or_eta = 0;
    TString str = deta_dphi;
    if( str.CompareTo("phi") == 0 ) phi_or_eta = 1; 
    if( str.CompareTo("eta") == 0 ) phi_or_eta = 2; 
    sw_up->SetParameter(0, SW_func1(phi_or_eta, nth_sw-1, 0, 0));
    sw_up->SetParameter(1, SW_func1(phi_or_eta, nth_sw-1, 0, 1));
    sw_up->SetParameter(2, SW_func1(phi_or_eta, nth_sw-1, 0, 2));
    sw_up->SetParameter(3, SW_func1(phi_or_eta, nth_sw-1, 0, 3));
    sw_up->SetParameter(4, SW_func1(phi_or_eta, nth_sw-1, 0, 4));

    sw_up->SetParameter(5, SW_func1(phi_or_eta, nth_sw-1, 1, 0));
    sw_up->SetParameter(6, SW_func1(phi_or_eta, nth_sw-1, 1, 1));
    sw_up->SetParameter(7, SW_func1(phi_or_eta, nth_sw-1, 1, 2));
    sw_up->SetParameter(8, SW_func1(phi_or_eta, nth_sw-1, 1, 3));
    sw_up->SetParameter(9, SW_func1(phi_or_eta, nth_sw-1, 1, 4));;
    sw_up->SetLineColor(kBlack);
    sw_up->Draw("lsame");

    sw_up->GetYaxis()->SetDecimals(3);
    sw_up->GetXaxis()->SetTitle("L1 EG E_{T} (crystal) [GeV]");
    sw_up->GetXaxis()->SetTitleSize(0.05);
    sw_up->GetYaxis()->SetTitle("#Delta#phi [rad.]");
    sw_up->GetYaxis()->SetTitleSize(0.05);
    sw_up->GetYaxis()->SetTitleOffset(1.35);
    sw_up->GetYaxis()->SetRangeUser(-0.005, 0.01);
    
    TF1 *sw_down = new TF1("sw_down","( [5]*pow(x,0) + [6]*pow(x,[7])*exp(-pow(x,[8])+[9])) - ( ( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4])) -( [5]*pow(x,0) + [6]*pow(x,[7])*exp(-pow(x,[8]) +[9])) ) / 2.", 10., 100.);
   
    sw_down->SetParameter(0, SW_func1(phi_or_eta, nth_sw-1, 0, 0));
    sw_down->SetParameter(1, SW_func1(phi_or_eta, nth_sw-1, 0, 1));
    sw_down->SetParameter(2, SW_func1(phi_or_eta, nth_sw-1, 0, 2));
    sw_down->SetParameter(3, SW_func1(phi_or_eta, nth_sw-1, 0, 3));
    sw_down->SetParameter(4, SW_func1(phi_or_eta, nth_sw-1, 0, 4));

    sw_down->SetParameter(5, SW_func1(phi_or_eta, nth_sw-1, 1, 0));
    sw_down->SetParameter(6, SW_func1(phi_or_eta, nth_sw-1, 1, 1));
    sw_down->SetParameter(7, SW_func1(phi_or_eta, nth_sw-1, 1, 2));
    sw_down->SetParameter(8, SW_func1(phi_or_eta, nth_sw-1, 1, 3));
    sw_down->SetParameter(9, SW_func1(phi_or_eta, nth_sw-1, 1, 4));;
    sw_down->SetLineColor(kBlack);
    sw_down->Draw("lsame");
 
    c1->SaveAs(dist + ".pdf");
    delete  c1; 
}
