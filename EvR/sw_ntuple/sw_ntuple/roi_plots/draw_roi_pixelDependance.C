void draw_roi_pixelDependance()
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

 TString pixel = "Pixel_5";
 TGraphErrors *roi; 

 int eta_ = 5;
 TString eta_region;
 eta_region.Form("%d", eta_);

 roi = (TGraphErrors*)infile1->Get(pixel+"eta_region_"+eta_region+"_median");
 
 roi->SetTitle("");
 roi->SetMarkerSize(0.5);
 roi->GetXaxis()->SetTitle("L1 EG E_{T} (crystal) [GeV]");
 roi->GetXaxis()->SetTitleSize(0.05);
 roi->GetXaxis()->SetRangeUser(10., 100.);
 roi->GetYaxis()->SetTitle("#Delta#phi [rad.]");
 roi->GetYaxis()->SetTitleSize(0.05);
 roi->GetYaxis()->SetRangeUser(-0.1, 0.05);
 roi->Draw("apeZ");

 for(int i = 2;i < 5; i++){
    int nth_pixel = i;
    if(eta_ == 2 && nth_pixel == 4) nth_pixel=nth_pixel+1; // D1 

    if(eta_ == 3 && nth_pixel == 3) nth_pixel=nth_pixel+2; // D1
    if(eta_ == 3 && nth_pixel == 4) nth_pixel=nth_pixel+2; // D2

    if(eta_ == 4 && nth_pixel == 2) nth_pixel=nth_pixel+3; // D1
    if(eta_ == 4 && nth_pixel == 3) nth_pixel=nth_pixel+3; // D2
    if(eta_ == 4 && nth_pixel == 4) nth_pixel=nth_pixel+3; // D2

    if(eta_ == 5 && nth_pixel == 2) nth_pixel=nth_pixel+4; // D1
    if(eta_ == 5 && nth_pixel == 3) nth_pixel=nth_pixel+4; // D2
    if(eta_ == 5 && nth_pixel == 4) nth_pixel=nth_pixel+4; // D2

    TString pixel_;
    pixel_.Form("%d", nth_pixel);

    roi = (TGraphErrors*)infile1->Get("Pixel_"+ pixel_ + "eta_region_"+eta_region+"_median");
    roi->SetMarkerColor(i+1);
    roi->Draw("samepleZ");
    cout << "pixel " << pixel_ << " printed..." << endl;
 }

 c1->SaveAs("eta_"+eta_region+"_roi_pixelDependence.pdf");

}
