void draw_median_fit(const std::string& median){

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  c1->SetRightMargin(0.1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetGridy();
  c1->SetGridx();
  c1->cd();

  cout << "string " << median << endl; 
  TGraphErrors* median_;
  median_ = (TGraphErrors*)gDirectory->Get((median).c_str() );

  median_->GetYaxis()->SetRangeUser(-0.1, 0.05);
  median_->SetTitle("");
  median_->Draw("apeZ");

}

void doFits(const std::string& median){

 TGraph* median_;
 median_ = (TGraph*)gDirectory->Get((median).c_str());

 TString nth_sw = (median).c_str();

 TF1 *median_fitFunc = new TF1("funcUp","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10., 100.);
 median_fitFunc->SetLineColor(kRed);
 median_fitFunc->SetLineStyle(2);

 median_->Fit(median_fitFunc,"0W");
 median_fitFunc->Draw("lsame");

 TCanvas* Current = gPad->GetCanvas();
 Current->SaveAs(nth_sw+"_roi.pdf");

 delete Current;

 // save the fit parameters
 ofstream fit_result;
 char fit_parameter[50];
 sprintf(fit_parameter, "./roi.txt");

 fit_result.open(fit_parameter);

 fit_result << endl;
 fit_result << "if( dp_de_dr == 0 && i == 6  && up_down == 0){" <<endl;
 for( int i=0; i < 5; i++){
     fit_result << "p[" << i << "] = " << median_fitFunc->GetParameter(i) << ";" << endl;
 }
 //fit_result << "return p[0]*pow(x,0) + p[1]*pow(x,p[2])*exp(-pow(x,p[3])+p[4]);" << endl;
 fit_result << "}" << endl;
 fit_result << endl;

}
