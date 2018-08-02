void draw_roi_pixelDependance()
{

 TFile* infile1 = NULL;
 TString filename1 = "/Volumes/Samsung_T3/Pixel_plots/CMSSW_9_2_0/sw_ntuple/EM_SA_medians/roi_median.root";
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

 TString nth_sw = "EG_Pix_SW_1";
 TGraphErrors *roi; 

 int eta_ = 5;
 TString eta_region;
 eta_region.Form("%d", eta_);

 roi = (TGraphErrors*)infile1->Get(nth_sw+"_eta_region_"+eta_region+"_median");
 
 roi->SetTitle("");
 roi->SetMarkerSize(0.5);
 roi->GetXaxis()->SetTitle("L1 EG E_{T} (crystal) [GeV]");
 roi->GetXaxis()->SetTitleSize(0.05);
 roi->GetXaxis()->SetRangeUser(10., 100.);
 roi->GetYaxis()->SetTitle("#Delta#phi [rad.]");
 roi->GetYaxis()->SetTitleSize(0.05);
 roi->GetYaxis()->SetRangeUser(-0.15, 0.15);
 roi->Draw("apeZ");

 for(int i = 2;i < 7; i++){
    int nth_sw = i;

    TString nth_sw_;
    nth_sw_.Form("%d", nth_sw);

    roi = (TGraphErrors*)infile1->Get("EG_Pix_SW_"+ nth_sw_ + "_eta_region_"+eta_region+"_median");
    roi->SetMarkerColor(i+1);
    roi->Draw("samepleZ");
    cout << "nth_sw " << nth_sw_ << " printed..." << endl;
 }

 c1->SaveAs("Eta_"+eta_region+"_swDependence.pdf");

}
