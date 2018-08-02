#define sw_ntuple_eff_cxx
#include "sw_ntuple_eff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <TLorentzVector.h>

// measure efficiency of each signal windows at 10 and 50 GeV

void sw_ntuple_eff::Loop()
{

   TFile *file = new TFile("eff_check.root","recreate");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float low_et = 0.;
   float high_et = 0.;

   int binSize = 90;
   vector<float> x[10][5], median[10][5], x_err[10][5], median_err[10][5];
   vector<float> deta_x[4][5], deta_median[4][5], deta_x_err[4][5], deta_median_err[4][5];

    
   for( int nth = 0; nth < 2; nth++){
      if(nth == 1)  nth = 40;
      low_et = 10. + nth;
      high_et = low_et + 1.;

      cout << "et: " << low_et << endl;


      vector<float> dPhi_pixV_pixV[10][5];
      vector<float> dEta_pixV_pixV[4][5];
      float abs_pterr;     
      float closest_egEt;  
      float closest_egEta; 
      float closest_egPhi; 
      float closest_eg_dr; 

      double dphi_width_[20] = {0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02};
      double deta_width_[20] = {0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02}; 


      for(int ith_width = 0; ith_width < 20; ith_width++){ 
      Long64_t nbytes = 0, nb = 0;
      cout << "dphi_width: " << dphi_width_[ith_width] << endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //for (Long64_t jentry=0; jentry<500000;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         // if (Cut(ientry) < 0) continue;

         //find closest egamma object to the gen electron
         int egN_ = egCrysClusterEt->size();
         int cl3d_N_ = cl3d_pt->size();
         if(egN_==0 && cl3d_N_==0) continue; // skip if there is no l1 eg objects at all 

         // barrel ECAL
         float closest_dr = 9999.;
         int closest_eg = 0;
         int egCount = 0;

         for(int i=0; i < egN_;i++){

            float dPhi = deltaPhi(propgenElPartPhi->at(0), egCrysClusterPhi->at(i));

            float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-egCrysClusterEta->at(i),2));
            if(egCrysClusterEt->at(i) > high_et || egCrysClusterEt->at(i) < low_et) continue;
            egCount++;
            if(current_dr < closest_dr){
              closest_dr = current_dr;
              closest_eg = i;
            }
         }// end of loop to find the closest egamma to gen electron 

         // HGCAL 3D cluster
         float closest_cl3d_dr = 9999.;
         int closest_cl3d = 0;
         int cl3d_Count = 0;

         for(int i=0; i < cl3d_N_;i++){

            float dPhi = deltaPhi(propgenElPartPhi->at(0), cl3d_phi->at(i));

            float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-cl3d_eta->at(i),2));
            if(cl3d_pt->at(i) > high_et || cl3d_pt->at(i) < low_et) continue;
            cl3d_Count++;
            if(current_dr < closest_cl3d_dr){
              closest_cl3d_dr = current_dr;
              closest_cl3d = i;
            }
         }// end of loop to find the closest egamma to gen electron 

         // skip if there is no barrel eg/ hgcal eg object in the current Et range
         if(egCount == 0 && cl3d_Count == 0) continue;
    
         TVector3 emvector;

         if( closest_dr < closest_cl3d_dr ){
            abs_pterr = fabs(genPartPt->at(0)-egCrysClusterEt->at(closest_eg))/genPartPt->at(0);
            closest_egEt = egCrysClusterEt->at(closest_eg);
            closest_egEta = egCrysClusterEta->at(closest_eg);
            closest_egPhi = egCrysClusterPhi->at(closest_eg);
            closest_eg_dr = closest_dr;
            emvector.SetXYZ(egCrysClusterGx->at(closest_eg),egCrysClusterGy->at(closest_eg), egCrysClusterGz->at(closest_eg));
         }
         else{
            abs_pterr = fabs(genPartPt->at(0)-cl3d_pt->at(closest_cl3d))/genPartPt->at(0);
            closest_egEt =  cl3d_pt->at(closest_cl3d);
            closest_egEta = cl3d_eta->at(closest_cl3d);
            closest_egPhi = cl3d_phi->at(closest_cl3d);
            closest_eg_dr = closest_cl3d_dr;
            emvector.SetXYZ(cl3d_x->at(closest_cl3d),cl3d_y->at(closest_cl3d), cl3d_z->at(closest_cl3d));
         }

         if(closest_eg_dr > 0.1) continue;

         // store pixel hits into vectors
         std::vector<TVector3> first_layer_hits;
         std::vector<TVector3> second_layer_hits;
         std::vector<TVector3> third_layer_hits;
         std::vector<TVector3> fourth_layer_hits;

         std::vector<TVector3> first_disk_hits;
         std::vector<TVector3> second_disk_hits;
         std::vector<TVector3> third_disk_hits;
         std::vector<TVector3> fourth_disk_hits;
         std::vector<int> hitted_layers;

         int layers[9] = {};  // initialize as 0 for each event, save number of hits for each pixel layer 
         layers[0] = 1; // beam spot

         std::vector<TVector3> PiXTRK_first_hits;
         std::vector<TVector3> PiXTRK_second_hits;
         std::vector<TVector3> PiXTRK_third_hits;
         std::vector<TVector3> PiXTRK_fourth_hits;

         int PiXTRK_layers[5] = {};  // initialize as 0 for each event, save number of hits for each pixel layer 
         PiXTRK_layers[0] = 1; // beam spot
 
         int eta_region = 0;
         if( fabs(closest_egEta) < 0.8 ) eta_region =1;
         if( fabs(closest_egEta) < 1.4 && fabs(closest_egEta) > 0.8 ) eta_region =2;
         if( fabs(closest_egEta) < 1.8 && fabs(closest_egEta) > 1.4 ) eta_region =3;
         if( fabs(closest_egEta) < 2.4 && fabs(closest_egEta) > 1.8 ) eta_region =4;
         if( fabs(closest_egEta) < 3.0 && fabs(closest_egEta) > 2.4 ) eta_region =5;

         if( fabs(closest_egEta) > 3.0 ) continue;

         SetSingalBoundary(eta_region, closest_egEt, dphi_width_[ith_width], deta_width_[ith_width]);

         //set roi dphi cut
         float upper_roi = ROI_func(eta_region, closest_egEt) + 0.055; 
         float lower_roi = ROI_func(eta_region, closest_egEt) - 0.055; 

         int bpix_size = bRecHitGx->size();
         for(int a=0; a<bpix_size; a++){
           TVector3 current_hit;
           current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
           double Dphi = deltaPhi(current_hit.Phi(), closest_egPhi);

           if( Dphi > upper_roi || Dphi < lower_roi ) continue;

            if( bRecHitLayer->at(a) == 1 ){
               layers[1]++;
               first_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            }
            if( bRecHitLayer->at(a) == 2 ){
               layers[2]++;
               second_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            }
            if( bRecHitLayer->at(a) == 3 ){
               layers[3]++;
               third_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            }
            if( bRecHitLayer->at(a) == 4 ){
               layers[4]++;
               fourth_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            }
         }

         int fpix_size = fRecHitGx->size();
         for(int a=0; a<fpix_size; a++){

           TVector3 current_hit;
           current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
           double Dphi = deltaPhi(current_hit.Phi(), closest_egPhi);

           if( Dphi > upper_roi || Dphi < lower_roi ) continue;

            if( fRecHitDisk->at(a) == 1 ){
               layers[5]++;
               first_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            if( fRecHitDisk->at(a) == 2 ){
               layers[6]++;
               second_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            if( fRecHitDisk->at(a) == 3 ){
               layers[7]++;
               third_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            if( fRecHitDisk->at(a) == 4 ){
               layers[8]++;
               fourth_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
         }

       if( fabs(closest_egEta) < 0.8 ){

         PiXTRK_first_hits = first_layer_hits;
         PiXTRK_layers[1] = PiXTRK_first_hits.size();

         PiXTRK_second_hits = second_layer_hits;
         PiXTRK_layers[2] = PiXTRK_second_hits.size();

         PiXTRK_third_hits = third_layer_hits;
         PiXTRK_layers[3] = PiXTRK_third_hits.size();

         PiXTRK_fourth_hits = fourth_layer_hits;
         PiXTRK_layers[4] = PiXTRK_fourth_hits.size();
       } 

       if( fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4 ){

         PiXTRK_first_hits = first_layer_hits;
         PiXTRK_layers[1] = PiXTRK_first_hits.size();

         PiXTRK_second_hits = second_layer_hits;
         PiXTRK_layers[2] = PiXTRK_second_hits.size();

         PiXTRK_third_hits = third_layer_hits;
         PiXTRK_layers[3] = PiXTRK_third_hits.size();

         PiXTRK_fourth_hits = first_disk_hits;
         PiXTRK_layers[4] = PiXTRK_fourth_hits.size();
       }

       if( fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8 ){

         PiXTRK_first_hits = first_layer_hits;
         PiXTRK_layers[1] = PiXTRK_first_hits.size();

         PiXTRK_second_hits = second_layer_hits;
         PiXTRK_layers[2] = PiXTRK_second_hits.size();

         PiXTRK_third_hits = first_disk_hits;
         PiXTRK_layers[3] = PiXTRK_third_hits.size();

         PiXTRK_fourth_hits = second_disk_hits;
         PiXTRK_layers[4] = PiXTRK_fourth_hits.size();
       }

       if( fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4 ){

         PiXTRK_first_hits = first_layer_hits;
         PiXTRK_layers[1] = PiXTRK_first_hits.size();

         PiXTRK_second_hits = first_disk_hits;
         PiXTRK_layers[2] = PiXTRK_second_hits.size();

         PiXTRK_third_hits = second_disk_hits;
         PiXTRK_layers[3] = PiXTRK_third_hits.size();

         PiXTRK_fourth_hits = third_disk_hits;
         PiXTRK_layers[4] = PiXTRK_fourth_hits.size();
       }

       if( fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0 ){

         PiXTRK_first_hits = first_disk_hits;
         PiXTRK_layers[1] = PiXTRK_first_hits.size();

         PiXTRK_second_hits = second_disk_hits;
         PiXTRK_layers[2] = PiXTRK_second_hits.size();

         PiXTRK_third_hits = third_disk_hits;
         PiXTRK_layers[3] = PiXTRK_third_hits.size();

         PiXTRK_fourth_hits = fourth_disk_hits;
         PiXTRK_layers[4] = PiXTRK_fourth_hits.size();
       }

        int n_pixels = 0; // initialize as 0 for each event
        hitted_layers.push_back(0); // 0 mean beam spot i.e., (0,0,0)
        for( int i=1; i < 5; i++){
           if( PiXTRK_layers[i] != 0 ){
             hitted_layers.push_back(i); //check if hits on each barrel or disk exists
             n_pixels++;
           }
        }

       bool L012_dphi_Et10_11_denom[5] = {false}, L012_dphi_Et10_11_nom[5] = {false};
       bool L012_dphi_Et50_51_denom[5] = {false}, L012_dphi_Et50_51_nom[5] = {false};

       bool L013_dphi_Et10_11_denom[5] = {false}, L013_dphi_Et10_11_nom[5] = {false};
       bool L013_dphi_Et50_51_denom[5] = {false}, L013_dphi_Et50_51_nom[5] = {false};

       bool L014_dphi_Et10_11_denom[5] = {false}, L014_dphi_Et10_11_nom[5] = {false};
       bool L014_dphi_Et50_51_denom[5] = {false}, L014_dphi_Et50_51_nom[5] = {false};

       bool L023_dphi_Et10_11_denom[5] = {false}, L023_dphi_Et10_11_nom[5] = {false};
       bool L023_dphi_Et50_51_denom[5] = {false}, L023_dphi_Et50_51_nom[5] = {false};

       bool L024_dphi_Et10_11_denom[5] = {false}, L024_dphi_Et10_11_nom[5] = {false};
       bool L024_dphi_Et50_51_denom[5] = {false}, L024_dphi_Et50_51_nom[5] = {false};

       bool L034_dphi_Et10_11_denom[5] = {false}, L034_dphi_Et10_11_nom[5] = {false};
       bool L034_dphi_Et50_51_denom[5] = {false}, L034_dphi_Et50_51_nom[5] = {false};

       bool L123_dphi_Et10_11_denom[5] = {false}, L123_dphi_Et10_11_nom[5] = {false};
       bool L123_dphi_Et50_51_denom[5] = {false}, L123_dphi_Et50_51_nom[5] = {false};

       bool L124_dphi_Et10_11_denom[5] = {false}, L124_dphi_Et10_11_nom[5] = {false};
       bool L124_dphi_Et50_51_denom[5] = {false}, L124_dphi_Et50_51_nom[5] = {false};

       bool L134_dphi_Et10_11_denom[5] = {false}, L134_dphi_Et10_11_nom[5] = {false};
       bool L134_dphi_Et50_51_denom[5] = {false}, L134_dphi_Et50_51_nom[5] = {false};

       bool L234_dphi_Et10_11_denom[5] = {false}, L234_dphi_Et10_11_nom[5] = {false};
       bool L234_dphi_Et50_51_denom[5] = {false}, L234_dphi_Et50_51_nom[5] = {false};


       bool L123_deta_Et10_11_denom[5] = {false}, L123_deta_Et10_11_nom[5] = {false};
       bool L123_deta_Et50_51_denom[5] = {false}, L123_deta_Et50_51_nom[5] = {false};

       bool L124_deta_Et10_11_denom[5] = {false}, L124_deta_Et10_11_nom[5] = {false};
       bool L124_deta_Et50_51_denom[5] = {false}, L124_deta_Et50_51_nom[5] = {false};

       bool L134_deta_Et10_11_denom[5] = {false}, L134_deta_Et10_11_nom[5] = {false};
       bool L134_deta_Et50_51_denom[5] = {false}, L134_deta_Et50_51_nom[5] = {false};

       bool L234_deta_Et10_11_denom[5] = {false}, L234_deta_Et10_11_nom[5] = {false};
       bool L234_deta_Et50_51_denom[5] = {false}, L234_deta_Et50_51_nom[5] = {false};

       if( n_pixels >= 3 ){
            for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){ // hitted_layers.begin()+1: to avoid using beam spot
               for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
                   for ( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){

                      for( int k=0; k < PiXTRK_layers[*first_hit]; k++){
                        for( int i=0; i < PiXTRK_layers[*second_hit]; i++){
                           for( int j=0; j < PiXTRK_layers[*third_hit]; j++){
                               double dPhi = 0.;
                               double dEta = 0.;

                               if( *first_hit == 0 && *second_hit == 1 && *third_hit == 2 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_first_hits[i], PiXTRK_second_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L012_dphi_Et10_11_denom[0] = true; //FillHist("L012_dphi_Et10_11_R1_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);   
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et10_11_nom[0] = true; //FillHist("L012_dphi_Et10_11_R1_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
                                     }
                                     if(low_et == 50.){
                                       L012_dphi_Et50_51_denom[0] = true; // FillHist("L012_dphi_Et50_11_R1_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);   
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et50_51_nom[0] = true; // FillHist("L012_dphi_Et50_11_R1_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
                                      }

                                  }
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L012_dphi_Et10_11_denom[1] = true; 
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et10_11_nom[1] = true; 
                                     }
                                     if(low_et == 50.){
                                       L012_dphi_Et50_51_denom[1] = true;
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et50_51_nom[1] = true;
                                      }
                                  }
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L012_dphi_Et10_11_denom[2] = true; 
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et10_11_nom[2] = true;
                                     }
                                     if(low_et == 50.){
                                       L012_dphi_Et50_51_denom[2] = true;
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et50_51_nom[2] = true;
                                      }
                                  }
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L012_dphi_Et10_11_denom[3] = true; 
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et10_11_nom[3] = true;
                                     }
                                     if(low_et == 50.){
                                       L012_dphi_Et50_51_denom[3] = true;
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et50_51_nom[3] = true;
                                      }
                                  }
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L012_dphi_Et10_11_denom[4] = true; 
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et10_11_nom[4] = true;
                                     }
                                     if(low_et == 50.){
                                       L012_dphi_Et50_51_denom[4] = true;
                                       if(L012_DPhi_cut1 > dPhi && L012_DPhi_cut2 < dPhi) L012_dphi_Et50_51_nom[4] = true;
                                      }
                                  }
                               }// L012
                               if( *first_hit == 0 && *second_hit == 1 && *third_hit == 3 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_first_hits[i], PiXTRK_third_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L013_dphi_Et10_11_denom[0] = true; 
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et10_11_nom[0] = true;
                                     }
                                     if(low_et == 50.){
                                       L013_dphi_Et50_51_denom[0] = true;
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et50_51_nom[0] = true;
                                      }

                                  }
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L013_dphi_Et10_11_denom[1] = true; 
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et10_11_nom[1] = true;
                                     }
                                     if(low_et == 50.){
                                       L013_dphi_Et50_51_denom[1] = true;
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et50_51_nom[1] = true;
                                      }
                                  }
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L013_dphi_Et10_11_denom[2] = true;
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et10_11_nom[2] = true;
                                     }
                                     if(low_et == 50.){
                                       L013_dphi_Et50_51_denom[2] = true;
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et50_51_nom[2] = true;
                                      }
                                  }
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L013_dphi_Et10_11_denom[3] = true;
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et10_11_nom[3] = true;
                                     }
                                     if(low_et == 50.){
                                       L013_dphi_Et50_51_denom[3] = true;
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et50_51_nom[3] = true;
                                      }
                                  }
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L013_dphi_Et10_11_denom[4] = true;
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et10_11_nom[4] = true;
                                     }
                                     if(low_et == 50.){
                                       L013_dphi_Et50_51_denom[4] = true;
                                       if(L013_DPhi_cut1 > dPhi && L013_DPhi_cut2 < dPhi) L013_dphi_Et50_51_nom[4] = true;
                                      }
                                  }
                               }// L013
                               if( *first_hit == 0 && *second_hit == 1 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_first_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L014_dphi_Et10_11_denom[0] = true; 
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et10_11_nom[0] = true;
                                     }   
                                     if(low_et == 50.){
                                       L014_dphi_Et50_51_denom[0] = true;
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et50_51_nom[0] = true;
                                      }   

                                  }   
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L014_dphi_Et10_11_denom[1] = true; 
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et10_11_nom[1] = true;
                                     }   
                                     if(low_et == 50.){
                                       L014_dphi_Et50_51_denom[1] = true;
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et50_51_nom[1] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L014_dphi_Et10_11_denom[2] = true;
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et10_11_nom[2] = true;
                                     }   
                                     if(low_et == 50.){
                                       L014_dphi_Et50_51_denom[2] = true;
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et50_51_nom[2] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L014_dphi_Et10_11_denom[3] = true;
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et10_11_nom[3] = true;
                                     }   
                                     if(low_et == 50.){
                                       L014_dphi_Et50_51_denom[3] = true;
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et50_51_nom[3] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L014_dphi_Et10_11_denom[4] = true;
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et10_11_nom[4] = true;
                                     }   
                                     if(low_et == 50.){
                                       L014_dphi_Et50_51_denom[4] = true;
                                       if(L014_DPhi_cut1 > dPhi && L014_DPhi_cut2 < dPhi) L014_dphi_Et50_51_nom[4] = true;
                                      }   
                                  }   
                               }
                               if( *first_hit == 0 && *second_hit == 2 && *third_hit == 3 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_second_hits[i], PiXTRK_third_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L023_dphi_Et10_11_denom[0] = true; 
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et10_11_nom[0] = true;
                                     }   
                                     if(low_et == 50.){
                                       L023_dphi_Et50_51_denom[0] = true;
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et50_51_nom[0] = true;
                                      }   

                                  }   
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L023_dphi_Et10_11_denom[1] = true; 
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et10_11_nom[1] = true;
                                     }   
                                     if(low_et == 50.){
                                       L023_dphi_Et50_51_denom[1] = true;
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et50_51_nom[1] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L023_dphi_Et10_11_denom[2] = true;
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et10_11_nom[2] = true;
                                     }   
                                     if(low_et == 50.){
                                       L023_dphi_Et50_51_denom[2] = true;
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et50_51_nom[2] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L023_dphi_Et10_11_denom[3] = true;
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et10_11_nom[3] = true;
                                     }   
                                     if(low_et == 50.){
                                       L023_dphi_Et50_51_denom[3] = true;
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et50_51_nom[3] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L023_dphi_Et10_11_denom[4] = true;
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et10_11_nom[4] = true;
                                     }   
                                     if(low_et == 50.){
                                       L023_dphi_Et50_51_denom[4] = true;
                                       if(L023_DPhi_cut1 > dPhi && L023_DPhi_cut2 < dPhi) L023_dphi_Et50_51_nom[4] = true;
                                      }   
                                  }   
                               }
                               if( *first_hit == 0 && *second_hit == 2 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L024_dphi_Et10_11_denom[0] = true; 
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et10_11_nom[0] = true;
                                     }   
                                     if(low_et == 50.){
                                       L024_dphi_Et50_51_denom[0] = true;
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et50_51_nom[0] = true;
                                      }   

                                  }   
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L024_dphi_Et10_11_denom[1] = true; 
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et10_11_nom[1] = true;
                                     }   
                                     if(low_et == 50.){
                                       L024_dphi_Et50_51_denom[1] = true;
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et50_51_nom[1] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L024_dphi_Et10_11_denom[2] = true;
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et10_11_nom[2] = true;
                                     }   
                                     if(low_et == 50.){
                                       L024_dphi_Et50_51_denom[2] = true;
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et50_51_nom[2] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L024_dphi_Et10_11_denom[3] = true;
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et10_11_nom[3] = true;
                                     }   
                                     if(low_et == 50.){
                                       L024_dphi_Et50_51_denom[3] = true;
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et50_51_nom[3] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L024_dphi_Et10_11_denom[4] = true;
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et10_11_nom[4] = true;
                                     }   
                                     if(low_et == 50.){
                                       L024_dphi_Et50_51_denom[4] = true;
                                       if(L024_DPhi_cut1 > dPhi && L024_DPhi_cut2 < dPhi) L024_dphi_Et50_51_nom[4] = true;
                                      }   
                                  }   
                               }
                               if( *first_hit == 0 && *second_hit == 3 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L034_dphi_Et10_11_denom[0] = true; 
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et10_11_nom[0] = true;
                                     }   
                                     if(low_et == 50.){
                                       L034_dphi_Et50_51_denom[0] = true;
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et50_51_nom[0] = true;
                                      }   

                                  }   
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L034_dphi_Et10_11_denom[1] = true; 
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et10_11_nom[1] = true;
                                     }   
                                     if(low_et == 50.){
                                       L034_dphi_Et50_51_denom[1] = true;
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et50_51_nom[1] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L034_dphi_Et10_11_denom[2] = true;
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et10_11_nom[2] = true;
                                     }   
                                     if(low_et == 50.){
                                       L034_dphi_Et50_51_denom[2] = true;
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et50_51_nom[2] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L034_dphi_Et10_11_denom[3] = true;
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et10_11_nom[3] = true;
                                     }   
                                     if(low_et == 50.){
                                       L034_dphi_Et50_51_denom[3] = true;
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et50_51_nom[3] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L034_dphi_Et10_11_denom[4] = true;
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et10_11_nom[4] = true;
                                     }   
                                     if(low_et == 50.){
                                       L034_dphi_Et50_51_denom[4] = true;
                                       if(L034_DPhi_cut1 > dPhi && L034_DPhi_cut2 < dPhi) L034_dphi_Et50_51_nom[4] = true;
                                      }   
                                  }   
                               }

                               if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){
                                  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_third_hits[j]);
                                  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_third_hits[j]);

                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L123_dphi_Et10_11_denom[0] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et10_11_nom[0] = true;
                                     }
                                     if(low_et == 50.){
                                       L123_dphi_Et50_51_denom[0] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et50_51_nom[0] = true;
                                      }

                                  }
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L123_dphi_Et10_11_denom[1] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et10_11_nom[1] = true;
                                     }
                                     if(low_et == 50.){
                                       L123_dphi_Et50_51_denom[1] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et50_51_nom[1] = true;
                                      }
                                  }
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L123_dphi_Et10_11_denom[2] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et10_11_nom[2] = true;
                                     }   
                                     if(low_et == 50.){
                                       L123_dphi_Et50_51_denom[2] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et50_51_nom[2] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L123_dphi_Et10_11_denom[3] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et10_11_nom[3] = true;
                                     }   
                                     if(low_et == 50.){
                                       L123_dphi_Et50_51_denom[3] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et50_51_nom[3] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L123_dphi_Et10_11_denom[4] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et10_11_nom[4] = true;
                                     }   
                                     if(low_et == 50.){
                                       L123_dphi_Et50_51_denom[4] = true;
                                       if(L123_DPhi_cut1 > dPhi && L123_DPhi_cut2 < dPhi) L123_dphi_Et50_51_nom[4] = true;
                                      }   
                                  }  

                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L123_deta_Et10_11_denom[0] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et10_11_nom[0] = true;
                                     }
                                     if(low_et == 50.){
                                       L123_deta_Et50_51_denom[0] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et50_51_nom[0] = true;
                                      }

                                  }
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L123_deta_Et10_11_denom[1] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et10_11_nom[1] = true;
                                     }
                                     if(low_et == 50.){
                                       L123_deta_Et50_51_denom[1] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et50_51_nom[1] = true;
                                      }
                                  }
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L123_deta_Et10_11_denom[2] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et10_11_nom[2] = true;
                                     }  
                                     if(low_et == 50.){
                                       L123_deta_Et50_51_denom[2] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et50_51_nom[2] = true;
                                      }  
                                  }  
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L123_deta_Et10_11_denom[3] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et10_11_nom[3] = true;
                                     }  
                                     if(low_et == 50.){
                                       L123_deta_Et50_51_denom[3] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et50_51_nom[3] = true;
                                      }  
                                  }  
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L123_deta_Et10_11_denom[4] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et10_11_nom[4] = true;
                                     }  
                                     if(low_et == 50.){
                                       L123_deta_Et50_51_denom[4] = true;
                                       if(L123_DEta_cut1 > dEta && L123_DEta_cut2 < dEta) L123_deta_Et50_51_nom[4] = true;
                                      }  
                                  } 
                               }
                               if( *first_hit == 1 && *second_hit == 2 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);
                                  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);

                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L124_dphi_Et10_11_denom[0] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et10_11_nom[0] = true;
                                     }   
                                     if(low_et == 50.){
                                       L124_dphi_Et50_51_denom[0] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et50_51_nom[0] = true;
                                      }   

                                  }   
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L124_dphi_Et10_11_denom[1] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et10_11_nom[1] = true;
                                     }   
                                     if(low_et == 50.){
                                       L124_dphi_Et50_51_denom[1] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et50_51_nom[1] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L124_dphi_Et10_11_denom[2] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et10_11_nom[2] = true;
                                     }   
                                     if(low_et == 50.){
                                       L124_dphi_Et50_51_denom[2] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et50_51_nom[2] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L124_dphi_Et10_11_denom[3] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et10_11_nom[3] = true;
                                     }   
                                     if(low_et == 50.){
                                       L124_dphi_Et50_51_denom[3] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et50_51_nom[3] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L124_dphi_Et10_11_denom[4] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et10_11_nom[4] = true;
                                     }   
                                     if(low_et == 50.){
                                       L124_dphi_Et50_51_denom[4] = true;
                                       if(L124_DPhi_cut1 > dPhi && L124_DPhi_cut2 < dPhi) L124_dphi_Et50_51_nom[4] = true;
                                      }   
                                  }   

                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L124_deta_Et10_11_denom[0] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et10_11_nom[0] = true;
                                     }    
                                     if(low_et == 50.){
                                       L124_deta_Et50_51_denom[0] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et50_51_nom[0] = true;
                                      }    

                                  }    
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L124_deta_Et10_11_denom[1] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et10_11_nom[1] = true;
                                     }    
                                     if(low_et == 50.){
                                       L124_deta_Et50_51_denom[1] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et50_51_nom[1] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L124_deta_Et10_11_denom[2] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et10_11_nom[2] = true;
                                     }    
                                     if(low_et == 50.){
                                       L124_deta_Et50_51_denom[2] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et50_51_nom[2] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L124_deta_Et10_11_denom[3] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et10_11_nom[3] = true;
                                     }    
                                     if(low_et == 50.){
                                       L124_deta_Et50_51_denom[3] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et50_51_nom[3] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L124_deta_Et10_11_denom[4] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et10_11_nom[4] = true;
                                     }    
                                     if(low_et == 50.){
                                       L124_deta_Et50_51_denom[4] = true;
                                       if(L124_DEta_cut1 > dEta && L124_DEta_cut2 < dEta) L124_deta_Et50_51_nom[4] = true;
                                      }    
                                  }    
                               }
                               if( *first_hit == 1 && *second_hit == 3 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L134_dphi_Et10_11_denom[0] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et10_11_nom[0] = true;
                                     }   
                                     if(low_et == 50.){
                                       L134_dphi_Et50_51_denom[0] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et50_51_nom[0] = true;
                                      }   

                                  }   
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L134_dphi_Et10_11_denom[1] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et10_11_nom[1] = true;
                                     }   
                                     if(low_et == 50.){
                                       L134_dphi_Et50_51_denom[1] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et50_51_nom[1] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L134_dphi_Et10_11_denom[2] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et10_11_nom[2] = true;
                                     }   
                                     if(low_et == 50.){
                                       L134_dphi_Et50_51_denom[2] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et50_51_nom[2] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L134_dphi_Et10_11_denom[3] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et10_11_nom[3] = true;
                                     }   
                                     if(low_et == 50.){
                                       L134_dphi_Et50_51_denom[3] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et50_51_nom[3] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L134_dphi_Et10_11_denom[4] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et10_11_nom[4] = true;
                                     }   
                                     if(low_et == 50.){
                                       L134_dphi_Et50_51_denom[4] = true;
                                       if(L134_DPhi_cut1 > dPhi && L134_DPhi_cut2 < dPhi) L134_dphi_Et50_51_nom[4] = true;
                                      }   
                                  }   

                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L134_deta_Et10_11_denom[0] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et10_11_nom[0] = true;
                                     }    
                                     if(low_et == 50.){
                                       L134_deta_Et50_51_denom[0] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et50_51_nom[0] = true;
                                      }    

                                  }    
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L134_deta_Et10_11_denom[1] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et10_11_nom[1] = true;
                                     }    
                                     if(low_et == 50.){
                                       L134_deta_Et50_51_denom[1] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et50_51_nom[1] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L134_deta_Et10_11_denom[2] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et10_11_nom[2] = true;
                                     }    
                                     if(low_et == 50.){
                                       L134_deta_Et50_51_denom[2] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et50_51_nom[2] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L134_deta_Et10_11_denom[3] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et10_11_nom[3] = true;
                                     }    
                                     if(low_et == 50.){
                                       L134_deta_Et50_51_denom[3] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et50_51_nom[3] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L134_deta_Et10_11_denom[4] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et10_11_nom[4] = true;
                                     }    
                                     if(low_et == 50.){
                                       L134_deta_Et50_51_denom[4] = true;
                                       if(L134_DEta_cut1 > dEta && L134_DEta_cut2 < dEta) L134_deta_Et50_51_nom[4] = true;
                                      }    
                                  }  
                               }
                               if( *first_hit == 2 && *second_hit == 3 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi( PiXTRK_second_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  dEta = StandaloneDEta( PiXTRK_second_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L234_dphi_Et10_11_denom[0] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et10_11_nom[0] = true;
                                     }   
                                     if(low_et == 50.){
                                       L234_dphi_Et50_51_denom[0] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et50_51_nom[0] = true;
                                      }   

                                  }   
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L234_dphi_Et10_11_denom[1] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et10_11_nom[1] = true;
                                     }   
                                     if(low_et == 50.){
                                       L234_dphi_Et50_51_denom[1] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et50_51_nom[1] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L234_dphi_Et10_11_denom[2] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et10_11_nom[2] = true;
                                     }   
                                     if(low_et == 50.){
                                       L234_dphi_Et50_51_denom[2] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et50_51_nom[2] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L234_dphi_Et10_11_denom[3] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et10_11_nom[3] = true;
                                     }   
                                     if(low_et == 50.){
                                       L234_dphi_Et50_51_denom[3] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et50_51_nom[3] = true;
                                      }   
                                  }   
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L234_dphi_Et10_11_denom[4] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et10_11_nom[4] = true;
                                     }   
                                     if(low_et == 50.){
                                       L234_dphi_Et50_51_denom[4] = true;
                                       if(L234_DPhi_cut1 > dPhi && L234_DPhi_cut2 < dPhi) L234_dphi_Et50_51_nom[4] = true;
                                      }   
                                  }   

                                  if(fabs(closest_egEta) < 0.8) {
                                     if(low_et == 10.){
                                        L234_deta_Et10_11_denom[0] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et10_11_nom[0] = true;
                                     }    
                                     if(low_et == 50.){
                                       L234_deta_Et50_51_denom[0] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et50_51_nom[0] = true;
                                      }    

                                  }    
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {
                                     if(low_et == 10.){
                                        L234_deta_Et10_11_denom[1] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et10_11_nom[1] = true;
                                     }    
                                     if(low_et == 50.){
                                       L234_deta_Et50_51_denom[1] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et50_51_nom[1] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {
                                     if(low_et == 10.){
                                        L234_deta_Et10_11_denom[2] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et10_11_nom[2] = true;
                                     }    
                                     if(low_et == 50.){
                                       L234_deta_Et50_51_denom[2] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et50_51_nom[2] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {
                                     if(low_et == 10.){
                                        L234_deta_Et10_11_denom[3] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et10_11_nom[3] = true;
                                     }    
                                     if(low_et == 50.){
                                       L234_deta_Et50_51_denom[3] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et50_51_nom[3] = true;
                                      }    
                                  }    
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {

                                     if(low_et == 10.){
                                        L234_deta_Et10_11_denom[4] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et10_11_nom[4] = true;
                                     }    
                                     if(low_et == 50.){
                                       L234_deta_Et50_51_denom[4] = true;
                                       if(L234_DEta_cut1 > dEta && L234_DEta_cut2 < dEta) L234_deta_Et50_51_nom[4] = true;
                                      }    
                                  }   
                               }





                           }// jth hit in third pixel layer/disk
                        }// ith hit in second pixel layer/disk 
                      }// kth hit in first pixel layer/disk
                 }
               }
             }
       }

       for(int i = 0; i < 5; i++){
          TString region_; 
          region_.Form("%d", i+1); 
          if(L012_dphi_Et10_11_denom[i]) FillHist("L012_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L012_dphi_Et10_11_nom[i]) FillHist("L012_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L012_dphi_Et50_51_denom[i]) FillHist("L012_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L012_dphi_Et50_51_nom[i]) FillHist("L012_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }

       for(int i = 0; i < 5; i++){
          TString region_; 
          region_.Form("%d", i+1); 
          if(L013_dphi_Et10_11_denom[i]) FillHist("L013_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L013_dphi_Et10_11_nom[i]) FillHist("L013_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L013_dphi_Et50_51_denom[i]) FillHist("L013_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L013_dphi_Et50_51_nom[i]) FillHist("L013_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }  

       for(int i = 0; i < 5; i++){
          TString region_; 
          region_.Form("%d", i+1); 
          if(L014_dphi_Et10_11_denom[i]) FillHist("L014_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L014_dphi_Et10_11_nom[i]) FillHist("L014_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L014_dphi_Et50_51_denom[i]) FillHist("L014_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L014_dphi_Et50_51_nom[i]) FillHist("L014_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }  

       for(int i = 0; i < 5; i++){
          TString region_; 
          region_.Form("%d", i+1); 
          if(L023_dphi_Et10_11_denom[i]) FillHist("L023_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L023_dphi_Et10_11_nom[i]) FillHist("L023_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L023_dphi_Et50_51_denom[i]) FillHist("L023_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L023_dphi_Et50_51_nom[i]) FillHist("L023_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }  

       for(int i = 0; i < 5; i++){
          TString region_; 
          region_.Form("%d", i+1); 
          if(L024_dphi_Et10_11_denom[i]) FillHist("L024_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L024_dphi_Et10_11_nom[i]) FillHist("L024_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L024_dphi_Et50_51_denom[i]) FillHist("L024_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L024_dphi_Et50_51_nom[i]) FillHist("L024_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }  

       for(int i = 0; i < 5; i++){
          TString region_; 
          region_.Form("%d", i+1); 
          if(L034_dphi_Et10_11_denom[i]) FillHist("L034_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L034_dphi_Et10_11_nom[i]) FillHist("L034_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L034_dphi_Et50_51_denom[i]) FillHist("L034_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L034_dphi_Et50_51_nom[i]) FillHist("L034_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }  

       for(int i = 0; i < 5; i++){
          TString region_; 
          region_.Form("%d", i+1); 
          if(L123_dphi_Et10_11_denom[i]) FillHist("L123_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L123_dphi_Et10_11_nom[i]) FillHist("L123_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L123_dphi_Et50_51_denom[i]) FillHist("L123_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L123_dphi_Et50_51_nom[i]) FillHist("L123_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }  

       for(int i = 0; i < 5; i++){
          TString region_;
          region_.Form("%d", i+1);
          if(L123_deta_Et10_11_denom[i]) FillHist("L123_deta_Et10_11_R"+region_+"_denom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L123_deta_Et10_11_nom[i]) FillHist("L123_deta_Et10_11_R"+region_+"_nom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L123_deta_Et50_51_denom[i]) FillHist("L123_deta_Et50_51_R"+region_+"_denom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L123_deta_Et50_51_nom[i]) FillHist("L123_deta_Et50_51_R"+region_+"_nom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }

       for(int i = 0; i < 5; i++){
          TString region_;
          region_.Form("%d", i+1);
          if(L124_dphi_Et10_11_denom[i]) FillHist("L124_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L124_dphi_Et10_11_nom[i]) FillHist("L124_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L124_dphi_Et50_51_denom[i]) FillHist("L124_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L124_dphi_Et50_51_nom[i]) FillHist("L124_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }

       for(int i = 0; i < 5; i++){
          TString region_;
          region_.Form("%d", i+1);
          if(L124_deta_Et10_11_denom[i]) FillHist("L124_deta_Et10_11_R"+region_+"_denom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L124_deta_Et10_11_nom[i]) FillHist("L124_deta_Et10_11_R"+region_+"_nom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L124_deta_Et50_51_denom[i]) FillHist("L124_deta_Et50_51_R"+region_+"_denom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L124_deta_Et50_51_nom[i]) FillHist("L124_deta_Et50_51_R"+region_+"_nom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }


       for(int i = 0; i < 5; i++){
          TString region_;
          region_.Form("%d", i+1);
          if(L134_dphi_Et10_11_denom[i]) FillHist("L134_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L134_dphi_Et10_11_nom[i]) FillHist("L134_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L134_dphi_Et50_51_denom[i]) FillHist("L134_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L134_dphi_Et50_51_nom[i]) FillHist("L134_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }

       for(int i = 0; i < 5; i++){
          TString region_;
          region_.Form("%d", i+1);
          if(L134_deta_Et10_11_denom[i]) FillHist("L134_deta_Et10_11_R"+region_+"_denom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L134_deta_Et10_11_nom[i]) FillHist("L134_deta_Et10_11_R"+region_+"_nom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L134_deta_Et50_51_denom[i]) FillHist("L134_deta_Et50_51_R"+region_+"_denom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L134_deta_Et50_51_nom[i]) FillHist("L134_deta_Et50_51_R"+region_+"_nom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }


       for(int i = 0; i < 5; i++){
          TString region_;
          region_.Form("%d", i+1);
          if(L234_dphi_Et10_11_denom[i]) FillHist("L234_dphi_Et10_11_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L234_dphi_Et10_11_nom[i]) FillHist("L234_dphi_Et10_11_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L234_dphi_Et50_51_denom[i]) FillHist("L234_dphi_Et50_51_R"+region_+"_denom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L234_dphi_Et50_51_nom[i]) FillHist("L234_dphi_Et50_51_R"+region_+"_nom", dphi_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }
       for(int i = 0; i < 5; i++){
          TString region_;
          region_.Form("%d", i+1);
          if(L234_deta_Et10_11_denom[i]) FillHist("L234_deta_Et10_11_R"+region_+"_denom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L234_deta_Et10_11_nom[i]) FillHist("L234_deta_Et10_11_R"+region_+"_nom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);

          if(L234_deta_Et50_51_denom[i]) FillHist("L234_deta_Et50_51_R"+region_+"_denom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
          if(L234_deta_Et50_51_nom[i]) FillHist("L234_deta_Et50_51_R"+region_+"_nom", deta_width_[ith_width]+0.0001, 1., 0.00, 0.02, 20);
       }

      }// event loop
      }

   
  }// et loop

  file->Write();

}
