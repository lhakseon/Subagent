
double SW_func1( int dp_de_dr, int i, int up_down, int nth_par){
  double p[5];

  if( dp_de_dr == 1 && i == 0 && up_down == 0){
  p[0] = 0.000975278;
  p[1] = 0.0133;
  p[2] = -0.388549;
  p[3] = 0.250104;
  p[4] = 1.19152;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 0 && up_down == 1){
  p[0] = -0.00193819;
  p[1] = 0.0133738;
  p[2] = -0.271429;
  p[3] = 0.0942481;
  p[4] = 0.566483;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 1 && up_down == 0){
  p[0] = -0.000521219;
  p[1] = 0.0178876;
  p[2] = -0.106844;
  p[3] = 0.194591;
  p[4] = 0.749352;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 1 && up_down == 1){
  p[0] = -0.00287549;
  p[1] = 0.017569;
  p[2] = -0.249616;
  p[3] = 0.0878684;
  p[4] = 0.770004;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 2 && up_down == 0){
  p[0] = 0.00114637;
  p[1] = 0.0208553;
  p[2] = -0.199312;
  p[3] = 0.267869;
  p[4] = 1.27016;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 2 && up_down == 1){
  p[0] = -6.83996e-05;
  p[1] = 0.028109;
  p[2] = -0.121331;
  p[3] = 0.348718;
  p[4] = 0.971395;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 3 && up_down == 0){
  p[0] = 0.0005629;
  p[1] = 0.0199326;
  p[2] = -0.14783;
  p[3] = 0.257888;
  p[4] = 0.91573;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 3 && up_down == 1){
  p[0] = -0.0030603;
  p[1] = 0.0225859;
  p[2] = -0.367993;
  p[3] = -0.0434304;
  p[4] = 0.524482;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 4 && up_down == 0){
  p[0] = 0.00142613;
  p[1] = 0.0214952;
  p[2] = -0.197876;
  p[3] = 0.272952;
  p[4] = 1.27689;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 4 && up_down == 1){
  p[0] = -6.44681e-05;
  p[1] = 0.0246838;
  p[2] = -0.122363;
  p[3] = 0.317283;
  p[4] = 0.898769;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 5 && up_down == 0){
  p[0] = 0.00240379;
  p[1] = 0.0184934;
  p[2] = -0.303992;
  p[3] = 0.225401;
  p[4] = 1.47589;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 5 && up_down == 1){
  p[0] = -0.000177651;
  p[1] = 0.0282997;
  p[2] = -0.116189;
  p[3] = 0.347584;
  p[4] = 0.921057;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 6 && up_down == 0){
  p[0] = 0.00139361;
  p[1] = 0.0153155;
  p[2] = -0.42477;
  p[3] = 0.236272;
  p[4] = 1.49141;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 6 && up_down == 1){
  p[0] = -0.00229594;
  p[1] = 0.0153785;
  p[2] = -0.237471;
  p[3] = 0.0914285;
  p[4] = 0.580549;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 7 && up_down == 0){
  p[0] = 0.00145951;
  p[1] = 0.0191037;
  p[2] = -0.214232;
  p[3] = 0.233397;
  p[4] = 1.15065;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 7 && up_down == 1){
  p[0] = -0.000230043;
  p[1] = 0.0263307;
  p[2] = -0.111647;
  p[3] = 0.333713;
  p[4] = 0.729226;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 8 && up_down == 0){
  p[0] = 0.00279709;
  p[1] = 0.0171982;
  p[2] = -0.276193;
  p[3] = 0.211584;
  p[4] = 1.0763;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 8 && up_down == 1){
  p[0] = -0.000925045;
  p[1] = 0.0284484;
  p[2] = -0.0918623;
  p[3] = 0.316092;
  p[4] = 0.555009;
  return p[nth_par];
  }
  
  if( dp_de_dr == 1 && i == 9 && up_down == 0){
  p[0] = 0.000780196;
  p[1] = 0.0134127;
  p[2] = -0.107403;
  p[3] = 0.0681098;
  p[4] = 0.658363;
  return p[nth_par];
  }
  if( dp_de_dr == 1 && i == 9 && up_down == 1){
  p[0] = -0.00262261;
  p[1] = 0.0123806;
  p[2] = 0.12444;
  p[3] = 0.291816;
  p[4] = 0.688309;
  return p[nth_par];
  }
  
  if( dp_de_dr == 2 && i == 0 && up_down == 0){
  p[0] = -0.000260009;
  p[1] = 0.00907667;
  p[2] = -0.0173921;
  p[3] = 0.0218334;
  p[4] = -0.073583;
  return p[nth_par];
  }
  if( dp_de_dr == 2 && i == 0 && up_down == 1){
  p[0] = -0.00265178;
  p[1] = 0.00279172;
  p[2] = -0.0650998;
  p[3] = 0.17266;
  p[4] = -0.271637;
  return p[nth_par];
  }
  
  if( dp_de_dr == 2 && i == 1 && up_down == 0){
  p[0] = -0.000410777;
  p[1] = 0.00876401;
  p[2] = -0.0204496;
  p[3] = 0.0287879;
  p[4] = -0.0873883;
  return p[nth_par];
  }
  if( dp_de_dr == 2 && i == 1 && up_down == 1){
  p[0] = -0.00232068;
  p[1] = 0.00338415;
  p[2] = -0.0626242;
  p[3] = 0.158658;
  p[4] = -0.261477;
  return p[nth_par];
  }
  
  if( dp_de_dr == 2 && i == 2 && up_down == 0){
  p[0] = -0.000360437;
  p[1] = 0.00894169;
  p[2] = -0.0180675;
  p[3] = 0.0246982;
  p[4] = -0.07712;
  return p[nth_par];
  }
  if( dp_de_dr == 2 && i == 2 && up_down == 1){
  p[0] = -0.00246151;
  p[1] = 0.00311131;
  p[2] = -0.0637055;
  p[3] = 0.164837;
  p[4] = -0.266012;
  return p[nth_par];
  }
  
  if( dp_de_dr == 2 && i == 3 && up_down == 0){
  p[0] = -0.000155345;
  p[1] = 0.00946475;
  p[2] = -0.0108639;
  p[3] = 0.0123213;
  p[4] = -0.0469508;
  return p[nth_par];
  }
  if( dp_de_dr == 2 && i == 3 && up_down == 1){
  p[0] = -0.00316137;
  p[1] = 0.000798782;
  p[2] = 0.619524;
  p[3] = 0.149595;
  p[4] = -1.72321;
  return p[nth_par];
  }

 return -1;
}

void draw_median_fit(const std::string& deta_dphi, const std::string& nth_sw, const std::string& eta_region){
    
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
    
    TString median = "Pix_Pix_d";
    median = median + deta_dphi + "_SW_" + nth_sw + "_eta_region_" + eta_region + "_median";
    cout << "string " << median << endl;
    TGraphErrors* median_;
    median_ = (TGraphErrors*)gDirectory->Get((median));
    
    
    TString dist = "SAD";
    dist = dist + deta_dphi + "_" + nth_sw + "_eta_region_" + eta_region;
    cout << dist << endl;
    
    TH2F* SA_dist;
    SA_dist = (TH2F*)gDirectory->Get((dist));
    
    SA_dist->Draw("COLZ");
    c1->Update();
    TPaletteAxis *palette = (TPaletteAxis*)SA_dist->GetListOfFunctions()->FindObject("palette");
    float x2_palette = palette->GetX2NDC();
    float x1_palette = palette->GetX1NDC();
    palette->SetX2NDC( x2_palette - (x2_palette-x1_palette) * 0.15 );
    palette->SetX1NDC( x1_palette + (x2_palette-x1_palette) * 0.15 );
    
    SA_dist->GetYaxis()->SetDecimals(3);
    SA_dist->GetXaxis()->SetTitle("L1 EG E_{T} (crystal) [GeV]");
    SA_dist->GetXaxis()->SetTitleSize(0.05);
    //SA_dist->GetYaxis()->SetTitle("#Delta#phi (" + directory[nth_sw] + ") [rad.]");
    SA_dist->GetYaxis()->SetTitleSize(0.05);
    SA_dist->GetYaxis()->SetTitleOffset(1.35);
    SA_dist->GetYaxis()->SetRangeUser(-0.005, 0.01);
    
    median_->GetYaxis()->SetRangeUser(-0.005, 0.01);
    median_->SetTitle("");
    //median_->Draw("same peZ");
    
    TF1 *median_fitFunc = new TF1("funcUp","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);
    
    median_fitFunc->SetParLimits(0, -0.005, 0.005);
    median_fitFunc->SetParLimits(1, -0.001, 5);
    median_fitFunc->SetParLimits(2, -1., 0.001);
    median_fitFunc->SetParLimits(3, -0.5, 0.5);
    
    median_fitFunc->SetLineColor(kBlack);
    median_fitFunc->SetLineStyle(1);
    
    median_->Fit(median_fitFunc,"0W");
    median_fitFunc->Draw("lsame");
    
    TF1 *median_up = new TF1("median_up","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) + [4] )", 10., 100.);
    TF1 *median_down = new TF1("median_down","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) + [4] )", 10., 100.);
    for( int i=0; i < 4; i++){
        median_up->SetParameter(i,median_fitFunc->GetParameter(i));
        median_down->SetParameter(i,median_fitFunc->GetParameter(i));
    }
    median_up->SetParameter(4, 0.001);
    median_up->Draw("lsame");
    
    median_down->SetParameter(4, -0.001);
    median_down->Draw("lsame");
    
    TF1 *sw_down = new TF1("sw_down","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]))", 10., 100.);
    //p[0] = 0.000975278;
    //p[1] = 0.0133;
    //p[2] = -0.388549;
    //p[3] = 0.250104;
    //p[4] = 1.19152;
    
    //p[0] = -0.000521219;
    //p[1] = 0.0178876;
    //p[2] = -0.106844;
    //p[3] = 0.194591;
    //p[4] = 0.749352;
    //sw_down->SetParameter(0, 0.000975278);
    //sw_down->SetParameter(1, 0.0133);
    //sw_down->SetParameter(2, -0.388549);
    //sw_down->SetParameter(3, 0.250104);
    //sw_down->SetParameter(4, 1.19152);
    
    sw_down->SetParameter(0, -0.000521219);
    sw_down->SetParameter(1, 0.0178876);
    sw_down->SetParameter(2, -0.106844);
    sw_down->SetParameter(3, 0.194591);
    sw_down->SetParameter(4, 0.749352);
    
    sw_down->SetLineColor(kBlack);
    sw_down->Draw("lsame");
    
    TF1 *sw_up = new TF1("sw_down","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3]) +[4]))", 10., 100.);
    //p[0] = -0.00193819;
    //p[1] = 0.0133738;
    //p[2] = -0.271429;
    //p[3] = 0.0942481;
    //p[4] = 0.566483;
    //sw_up->SetParameter(0, -0.00193819);
    //sw_up->SetParameter(1, 0.0133738);
    //sw_up->SetParameter(2, -0.271429);
    //sw_up->SetParameter(3, 0.0942481);
    //sw_up->SetParameter(4, 0.566483);;
    
    //p[0] = -0.00287549;
    //p[1] = 0.017569;
    //p[2] = -0.249616;
    //p[3] = 0.0878684;
    //p[4] = 0.770004;
    
    sw_up->SetParameter(0,  -0.00287549);
    sw_up->SetParameter(1, 0.017569);
    sw_up->SetParameter(2, -0.249616);
    sw_up->SetParameter(3, 0.0878684);
    sw_up->SetParameter(4, 0.770004);;
    
    sw_up->SetLineColor(kBlack);
    sw_up->Draw("lsame");
    
    SW_func1(1, 0, 0, 0);

    cout << "test: " << SW_func1(1, 0, 0, 0) << endl;    

    delete  c1; 
}


// make plot of median points and their fit result
void doFits(const std::string& median){

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

 median_->GetYaxis()->SetRangeUser(-0.015, 0.015);
 median_->SetTitle("");
 median_->Draw("apeZ");

 TString nth_sw = (median).c_str();

 TF1 *median_fitFunc = new TF1("funcUp","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);
 //TF1 *median_fitFunc = new TF1("funcUp","( [0]*pow(x,0) + [1]*exp(-pow(x,[2])) )", 10., 100.);

 median_fitFunc->SetParLimits(0, -0.005, 0.005);
 median_fitFunc->SetParLimits(1, -0.001, 5);
 median_fitFunc->SetParLimits(2, -1., 0.001);
 median_fitFunc->SetParLimits(3, -0.5, 0.5);

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
 for( int i=0; i < 4; i++){
     fit_result << "p[" << i << "] = " << median_fitFunc->GetParameter(i) << ";" << endl;
     cout <<  "p[" << i << "] = " << median_fitFunc->GetParameter(i) << ";" << endl;
 }
 //fit_result << "return p[0]*pow(x,0) + p[1]*pow(x,p[2])*exp(-pow(x,p[3])+p[4]);" << endl;
 fit_result << "}" << endl;
 fit_result << endl;

}
