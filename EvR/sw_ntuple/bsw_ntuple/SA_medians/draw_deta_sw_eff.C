void draw_deta_sw_eff(int region, int et){

 TString et10 = "Et10_11_R";
 TString et50 = "Et50_51_R";

 TString sw[4] = {"L123_deta_", "L124_deta_", "L134_deta_", "L234_deta_"};

 TCanvas *c1 = new TCanvas("c1","c1",900,500);
 gStyle->SetOptStat(0);
 gStyle->SetLineWidth(1); // axis width, default is 1
 c1->SetTopMargin(0.05);
 c1->SetBottomMargin(0.12);
 c1->SetRightMargin(0.05);
 c1->SetLeftMargin(0.2);
 c1->SetGrid();
 c1->SetTicky(1);
 c1->SetTickx(1);
 c1->cd();

 TH1D* h_denom; 
 TH1D* h_nom;

 TLegend* leg = new TLegend(0.45, 0.20, 0.9, 0.7,"","brNDC");
 leg->SetTextSize(0.04);
 leg->SetLineWidth(1);
 leg->SetBorderSize(0);
 leg->SetFillStyle(0);

 for(int i = 0; i < 4; i++){

     c1->cd();
     TString region_;
     region_.Form("%d", region);

     TString denom_name; 
     TString nom_name;

     if(et == 10){
       denom_name = sw[i] + et10 + region_+"_denom";
       nom_name =   sw[i] + et10 + region_+"_nom";
     }
     if(et == 50){
       denom_name = sw[i] + et50 + region_+"_denom";
       nom_name =   sw[i] + et50 + region_+"_nom";
     }

     h_denom = (TH1D*)gDirectory->Get(denom_name);
     h_nom = (TH1D*)gDirectory->Get(nom_name);

     h_nom->Divide(h_denom);
     h_nom->SetTitle("");

     h_nom->GetYaxis()->SetTitleOffset(1.2);
     h_nom->GetYaxis()->SetTitleSize(0.05);
     h_nom->SetMarkerStyle(20);
     h_nom->SetMarkerSize(0.6);
     h_nom->GetYaxis()->SetTitle("Efficiency");
     h_nom->GetXaxis()->SetTitle("#Delta #eta cut");
     h_nom->GetYaxis()->SetRangeUser(0.0, 1.05);
     h_nom->GetXaxis()->SetTitleSize(0.05);
     //h_nom->GetXaxis()->SetLabelSize(0.05);
     //h_nom->GetYaxis()->SetLabelSize(0.05);
     h_nom->SetMarkerColor(i+1);
     if(i==9) h_nom->SetMarkerColor(i+2);
     h_nom->SetMarkerStyle(20);
     if(i==9) h_nom->SetLineColor(i+2);
     if(i==0) h_nom->Draw(); 
     if(i!=0) h_nom->Draw("same"); 
  
     if(et == 10) leg->AddEntry(h_nom, sw[i] + et10 + region_, "pe");
     if(et == 50) leg->AddEntry(h_nom, sw[i] + et50 + region_, "pe");
     leg->Draw();
 }

 TString region_;
 region_.Form("%d", region);
 if(et==10) c1->SaveAs( et10 + region_ + "_deta_eff.pdf");
 if(et==50) c1->SaveAs( et50 + region_ +"_deta_eff.pdf");
 delete c1;

}
