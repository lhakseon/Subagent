void draw_roi()
{

 TFile* infile1 = NULL;
 TString filename1 = "/Volumes/Samsung_T3/Pixel_plots/CMSSW_9_2_0/sw_ntuple/roi_plots/roi_median.root";
 infile1 = new TFile(filename1);


 TCanvas *c1 = new TCanvas("c1","c1",950,600);
 gStyle->SetOptStat(0);
 c1->SetRightMargin(0.1);
 c1->SetLeftMargin(0.15);
 c1->SetBottomMargin(0.15);
 c1->SetGridy();
 c1->SetGridx();
 c1->SetTicky(1);
 c1->SetTickx(1);

 TString pixel = "Pixel_2";
 TGraph *roi; 

 roi = (TGraph*)infile1->Get(pixel+"eta_region_1_median");
 
 roi->SetTitle("");
 roi->SetMarkerSize(0.5);
 roi->GetXaxis()->SetTitle("L1 EG E_{T} (crystal) [GeV]");
 roi->GetXaxis()->SetTitleSize(0.05);
 roi->GetXaxis()->SetRangeUser(10., 100.);
 roi->GetYaxis()->SetTitle("#Delta#phi [rad.]");
 roi->GetYaxis()->SetTitleSize(0.05);
 roi->GetYaxis()->SetRangeUser(-0.1, 0.1);
 roi->Draw("apeZ");

 for(int i = 1;i < 3; i++){
    TString eta_region;
    eta_region.Form("%d", i + 1);
    roi = (TGraph*)infile1->Get(pixel+"eta_region_"+eta_region+"_median");
    roi->SetMarkerColor(i+1);
    roi->Draw("samepleZ");
 }

 c1->SaveAs(pixel+"_roi.pdf");

}
