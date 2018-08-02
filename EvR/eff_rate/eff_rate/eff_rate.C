void eff_rate(){


 TFile* eff = new TFile("/cms/ldap_home/lhakseon/L1PT/t802/eff_hist_width.root");
 TFile* bkg_eff = new TFile("/cms/ldap_home/lhakseon/L1PT/t802/bkg_eff_hist_width.root");

 TGraphErrors* h_eff = (TGraphErrors*)eff->Get("pix");
 TGraphErrors* h_bkg_eff = (TGraphErrors*)bkg_eff->Get("pix");

 TGraphErrors* h_eg_eff = (TGraphErrors*)eff->Get("eg");
 TGraphErrors* h_eg_bkg_eff = (TGraphErrors*)bkg_eff->Get("eg");

 double *x = h_eff->GetX();

 double *eff_ = h_eff->GetY();
 double *eff_err_ = h_eff->GetEY();

 double *bkg_eff_ = h_bkg_eff->GetY();
 double *bkg_eff_err_ = h_bkg_eff->GetEY();

 double *eg_eff_ = h_eg_eff->GetY();
 double *eg_eff_err_ = h_eg_eff->GetEY();

 double *eg_bkg_eff_ = h_eg_bkg_eff->GetY();
 double *eg_bkg_eff_err_ = h_eg_bkg_eff->GetEY();

 double x_bkg_eff[27], x_bkg_eff_err[27];
 double y_eff[27], y_eff_err[27];

 double x_eg_bkg_eff[1], x_eg_bkg_eff_err[1];
 double y_eg_eff[1], y_eg_eff_err[1];

 x_eg_bkg_eff[0] = eg_bkg_eff_[0];
 x_eg_bkg_eff_err[0] = eg_bkg_eff_err_[0];
 
 y_eg_eff[0] = eg_eff_[0];
 y_eg_eff_err[0] = eg_eff_err_[0];

 for(int i = 0; i < 27; i++){

   cout << "eff: " << x[i] << " eff: " << eff_[i] << "+/-" << eff_err_[i] << " bkg eff: " << bkg_eff_[i] << "+/-" << bkg_eff_err_[i] << endl;

   x_bkg_eff[i] = bkg_eff_[i];
   x_bkg_eff_err[i] = bkg_eff_err_[i];
   
   y_eff[i] = eff_[i];
   y_eff_err[i] = eff_err_[i]; 

 }

   //TGraphErrors* eff_rate = new TGraphErrors(27, x_bkg_eff, y_eff, x_bkg_eff_err, y_eff_err);
   TGraphErrors* eff_rate = new TGraphErrors(27, y_eff, x_bkg_eff, y_eff_err, x_bkg_eff_err);
   TGraphErrors* eg_eff_rate = new TGraphErrors(1, y_eg_eff, x_eg_bkg_eff, y_eg_eff_err, x_eg_bkg_eff_err);

   TCanvas *c1 = new TCanvas("c1","c1",950,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(1); // axis width, default is 1
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.12);
   c1->SetRightMargin(0.03);
   c1->SetLeftMargin(0.2);
   c1->SetGrid();
   c1->SetTicky(1);
   c1->SetTickx(1);
   c1->SetLogy();
   c1->cd();

   eff_rate->SetTitle("");
   eff_rate->GetXaxis()->SetTitleOffset(1.);
   eff_rate->GetXaxis()->SetTitleSize(0.055);
   eff_rate->GetXaxis()->SetNdivisions(510);
   eff_rate->GetXaxis()->SetLabelSize(0.05);
   eff_rate->GetXaxis()->SetTitle("Efficiency");

   eff_rate->GetYaxis()->SetTitle("Rate (kHz)");
   eff_rate->GetYaxis()->SetNdivisions(506);
   eff_rate->GetYaxis()->SetLabelSize(0.05);
   //eff_rate->GetXaxis()->SetRangeUser(0.0, 0.1);
   //eff_rate->GetYaxis()->SetRangeUser(0., 1.2);

   eff_rate->GetXaxis()->SetRangeUser(.9, 1.0);
   eff_rate->GetYaxis()->SetRangeUser(9., 500.);
   eff_rate->GetYaxis()->SetTitleOffset(1.2);
   eff_rate->GetYaxis()->SetTitleSize(0.055);

   eff_rate->SetMarkerColor(4);
   eff_rate->SetLineColor(4);
   eff_rate->SetLineWidth(2);
   eff_rate->SetMarkerStyle(29);
   eff_rate->SetMarkerSize(1.8);
   eff_rate->Draw("ape");

   eg_eff_rate->SetMarkerColor(2);
   eg_eff_rate->SetLineColor(2);
   eg_eff_rate->Draw("same pe");

   c1->Print("Eff_Calorimeter.png");

}
