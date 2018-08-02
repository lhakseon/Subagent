#define sw_ntuple_cxx
#include "sw_ntuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <TLorentzVector.h>

void sw_ntuple::Loop()
{

   TFile *file = new TFile("roi_median_.root","recreate");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float low_et = 0.;
   float high_et = 0.;

   int binSize = 90;
   vector<float> x[10][5], median[10][5], x_err[10][5], median_err[10][5];
   vector<float> deta_x[4][5], deta_median[4][5], deta_x_err[4][5], deta_median_err[4][5];

   // for 2D distribution plots
   TH2F* SADphi_dist[10][5];
   TH2F* SADeta_dist[4][5];

   for(int j = 0; j < 10; j++){
      for(int i = 0; i < 5; i++){

         TString nth_eg_sw_;
         nth_eg_sw_.Form("%d", j + 1);
         TString eta_region;
         eta_region.Form("%d", i + 1);

         SADphi_dist[j][i] = new TH2F("","", 90,10,100,100,-0.02,0.02);
         SADphi_dist[j][i]->SetName( "SADphi_"+nth_eg_sw_+"_eta_region_"+eta_region);
      }
   }

   for(int j = 0; j < 4; j++){
      for(int i = 0; i < 5; i++){

         TString nth_eg_sw_;
         nth_eg_sw_.Form("%d", j + 1);
         TString eta_region;
         eta_region.Form("%d", i + 1);

         SADeta_dist[j][i] = new TH2F("","", 90,10,100,100,-0.02,0.02);
         SADeta_dist[j][i]->SetName( "SADeta_"+nth_eg_sw_+"_eta_region_"+eta_region);
      }
   }
 
   for( int nth = 0; nth < 90; nth++){
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
 
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //for (Long64_t jentry=0; jentry<10000;jentry++) {
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
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[0][0].push_back(dPhi); SADphi_dist[0][0]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[0][1].push_back(dPhi); SADphi_dist[0][1]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[0][2].push_back(dPhi); SADphi_dist[0][2]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[0][3].push_back(dPhi); SADphi_dist[0][3]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[0][4].push_back(dPhi); SADphi_dist[0][4]->Fill(closest_egEt,dPhi);}
                               }
                               if( *first_hit == 0 && *second_hit == 1 && *third_hit == 3 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_first_hits[i], PiXTRK_third_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[1][0].push_back(dPhi); SADphi_dist[1][0]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[1][1].push_back(dPhi); SADphi_dist[1][1]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[1][2].push_back(dPhi); SADphi_dist[1][2]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[1][3].push_back(dPhi); SADphi_dist[1][3]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[1][4].push_back(dPhi); SADphi_dist[1][4]->Fill(closest_egEt,dPhi);}
                               }
                               if( *first_hit == 0 && *second_hit == 1 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_first_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[2][0].push_back(dPhi); SADphi_dist[2][0]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[2][1].push_back(dPhi); SADphi_dist[2][1]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[2][2].push_back(dPhi); SADphi_dist[2][2]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[2][3].push_back(dPhi); SADphi_dist[2][3]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[2][4].push_back(dPhi); SADphi_dist[2][4]->Fill(closest_egEt,dPhi);}
                               }
                               if( *first_hit == 0 && *second_hit == 2 && *third_hit == 3 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_second_hits[i], PiXTRK_third_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[3][0].push_back(dPhi); SADphi_dist[3][0]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[3][1].push_back(dPhi); SADphi_dist[3][1]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[3][2].push_back(dPhi); SADphi_dist[3][2]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[3][3].push_back(dPhi); SADphi_dist[3][3]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[3][4].push_back(dPhi); SADphi_dist[3][4]->Fill(closest_egEt,dPhi);}
                               }
                               if( *first_hit == 0 && *second_hit == 2 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[4][0].push_back(dPhi); SADphi_dist[4][0]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[4][1].push_back(dPhi); SADphi_dist[4][1]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[4][2].push_back(dPhi); SADphi_dist[4][2]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[4][3].push_back(dPhi); SADphi_dist[4][3]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[4][4].push_back(dPhi); SADphi_dist[4][4]->Fill(closest_egEt,dPhi);}
                               }
                               if( *first_hit == 0 && *second_hit == 3 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi_BS( PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[5][0].push_back(dPhi); SADphi_dist[5][0]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[5][1].push_back(dPhi); SADphi_dist[5][1]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[5][2].push_back(dPhi); SADphi_dist[5][2]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[5][3].push_back(dPhi); SADphi_dist[5][3]->Fill(closest_egEt,dPhi);}
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[5][4].push_back(dPhi); SADphi_dist[5][4]->Fill(closest_egEt,dPhi);}
                               }

                               if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){
                                  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_third_hits[j]);
                                  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_third_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {
                                    dPhi_pixV_pixV[6][0].push_back(dPhi); dEta_pixV_pixV[0][0].push_back(dEta); 
                                    SADphi_dist[6][0]->Fill(closest_egEt,dPhi); SADeta_dist[0][0]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[6][1].push_back(dPhi); dEta_pixV_pixV[0][1].push_back(dEta);
                                    SADphi_dist[6][1]->Fill(closest_egEt,dPhi); SADeta_dist[0][1]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[6][2].push_back(dPhi); dEta_pixV_pixV[0][2].push_back(dEta);
                                    SADphi_dist[6][2]->Fill(closest_egEt,dPhi); SADeta_dist[0][2]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[6][3].push_back(dPhi); dEta_pixV_pixV[0][3].push_back(dEta);
                                    SADphi_dist[6][3]->Fill(closest_egEt,dPhi); SADeta_dist[0][3]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[6][4].push_back(dPhi); dEta_pixV_pixV[0][4].push_back(dEta);
                                    SADphi_dist[6][4]->Fill(closest_egEt,dPhi); SADeta_dist[0][4]->Fill(closest_egEt,dEta);
                                  }
                               }
                               if( *first_hit == 1 && *second_hit == 2 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);
                                  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_second_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[7][0].push_back(dPhi); dEta_pixV_pixV[1][0].push_back(dEta);
                                    SADphi_dist[7][0]->Fill(closest_egEt,dPhi); SADeta_dist[1][0]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[7][1].push_back(dPhi); dEta_pixV_pixV[1][1].push_back(dEta);
                                    SADphi_dist[7][1]->Fill(closest_egEt,dPhi); SADeta_dist[1][1]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[7][2].push_back(dPhi); dEta_pixV_pixV[1][2].push_back(dEta);
                                    SADphi_dist[7][2]->Fill(closest_egEt,dPhi); SADeta_dist[1][2]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[7][3].push_back(dPhi); dEta_pixV_pixV[1][3].push_back(dEta);
                                    SADphi_dist[7][3]->Fill(closest_egEt,dPhi); SADeta_dist[1][3]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[7][4].push_back(dPhi); dEta_pixV_pixV[1][4].push_back(dEta);
                                    SADphi_dist[7][3]->Fill(closest_egEt,dPhi); SADeta_dist[1][4]->Fill(closest_egEt,dEta);
                                  }
                               }
                               if( *first_hit == 1 && *second_hit == 3 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi( PiXTRK_first_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  dEta = StandaloneDEta( PiXTRK_first_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[8][0].push_back(dPhi); dEta_pixV_pixV[2][0].push_back(dEta);
                                    SADphi_dist[8][0]->Fill(closest_egEt,dPhi); SADeta_dist[2][0]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[8][1].push_back(dPhi); dEta_pixV_pixV[2][1].push_back(dEta);
                                    SADphi_dist[8][1]->Fill(closest_egEt,dPhi); SADeta_dist[2][1]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[8][2].push_back(dPhi); dEta_pixV_pixV[2][2].push_back(dEta);
                                    SADphi_dist[8][2]->Fill(closest_egEt,dPhi); SADeta_dist[2][2]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[8][3].push_back(dPhi); dEta_pixV_pixV[2][3].push_back(dEta);
                                    SADphi_dist[8][3]->Fill(closest_egEt,dPhi); SADeta_dist[2][3]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[8][4].push_back(dPhi); dEta_pixV_pixV[2][4].push_back(dEta);
                                    SADphi_dist[8][4]->Fill(closest_egEt,dPhi); SADeta_dist[2][4]->Fill(closest_egEt,dEta);
                                  }
                               }
                               if( *first_hit == 2 && *second_hit == 3 && *third_hit == 4 ){
                                  dPhi = StandaloneDPhi( PiXTRK_second_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  dEta = StandaloneDEta( PiXTRK_second_hits[k], PiXTRK_third_hits[i], PiXTRK_fourth_hits[j]);
                                  if(fabs(closest_egEta) < 0.8) {dPhi_pixV_pixV[9][0].push_back(dPhi); dEta_pixV_pixV[3][0].push_back(dEta);
                                    SADphi_dist[9][0]->Fill(closest_egEt,dPhi); SADeta_dist[3][0]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 0.8 && fabs(closest_egEta) < 1.4) {dPhi_pixV_pixV[9][1].push_back(dPhi); dEta_pixV_pixV[3][1].push_back(dEta);
                                    SADphi_dist[9][1]->Fill(closest_egEt,dPhi); SADeta_dist[3][1]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 1.4 && fabs(closest_egEta) < 1.8) {dPhi_pixV_pixV[9][2].push_back(dPhi); dEta_pixV_pixV[3][2].push_back(dEta);
                                    SADphi_dist[9][2]->Fill(closest_egEt,dPhi); SADeta_dist[3][2]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 1.8 && fabs(closest_egEta) < 2.4) {dPhi_pixV_pixV[9][3].push_back(dPhi); dEta_pixV_pixV[3][3].push_back(dEta);
                                    SADphi_dist[9][3]->Fill(closest_egEt,dPhi); SADeta_dist[3][3]->Fill(closest_egEt,dEta);
                                  }
                                  if(fabs(closest_egEta) > 2.4 && fabs(closest_egEta) < 3.0) {dPhi_pixV_pixV[9][4].push_back(dPhi); dEta_pixV_pixV[3][4].push_back(dEta);
                                    SADphi_dist[9][4]->Fill(closest_egEt,dPhi); SADeta_dist[3][4]->Fill(closest_egEt,dEta);
                                  }
                               }

                           }// jth hit in third pixel layer/disk
                        }// ith hit in second pixel layer/disk 
                      }// kth hit in first pixel layer/disk
                 }
               }
             }
       }

      }// event loop


      for(int nth_eg_sw = 0; nth_eg_sw < 10; nth_eg_sw++){
        for(int i = 0; i < 5; i++){
           std::sort (dPhi_pixV_pixV[nth_eg_sw][i].begin(), dPhi_pixV_pixV[nth_eg_sw][i].end());

           if( dPhi_pixV_pixV[nth_eg_sw][i].size() != 0 ){
             x[nth_eg_sw][i].push_back(low_et + 0.5);
             median[nth_eg_sw][i].push_back(getMedian(dPhi_pixV_pixV[nth_eg_sw][i]));
             x_err[nth_eg_sw][i].push_back(0.);
             median_err[nth_eg_sw][i].push_back(getMedianErr(dPhi_pixV_pixV[nth_eg_sw][i]));
           }
        }
      }

      for(int nth_eg_sw = 0; nth_eg_sw < 4; nth_eg_sw++){
        for(int i = 0; i < 5; i++){
           std::sort (dEta_pixV_pixV[nth_eg_sw][i].begin(), dEta_pixV_pixV[nth_eg_sw][i].end());

           if( dEta_pixV_pixV[nth_eg_sw][i].size() != 0 ){
             deta_x[nth_eg_sw][i].push_back(low_et + 0.5);
             deta_median[nth_eg_sw][i].push_back(getMedian(dEta_pixV_pixV[nth_eg_sw][i]));
             deta_x_err[nth_eg_sw][i].push_back(0.);
             deta_median_err[nth_eg_sw][i].push_back(getMedianErr(dEta_pixV_pixV[nth_eg_sw][i]));
           }
        }
      }
   
  }// et loop

  TGraphErrors* dphi_median_gr[10][5];
  for(int j = 0; j < 10; j++){
   for(int i = 0; i < 5; i++){
     int dphi_point_size = x[j][i].size();
     dphi_median_gr[j][i] = new TGraphErrors(dphi_point_size, &x[j][i][0], &median[j][i][0], &x_err[j][i][0], &median_err[j][i][0]);

    dphi_median_gr[j][i]->SetMarkerStyle(24);
    dphi_median_gr[j][i]->SetMarkerSize(0.3);
    TString nth_eg_sw_;
    nth_eg_sw_.Form("%d", j + 1);
    TString eta_region;
    eta_region.Form("%d", i + 1);
    dphi_median_gr[j][i]->SetName( "Pix_Pix_dphi_SW_"+nth_eg_sw_+"_eta_region_"+eta_region+"_median");

    dphi_median_gr[j][i]->Write();
   }
  }

  TGraphErrors* deta_median_gr[4][5];
  for(int j = 0; j < 4; j++){
   for(int i = 0; i < 5; i++){
     int deta_point_size = x[j][i].size();
     deta_median_gr[j][i] = new TGraphErrors(deta_point_size, &deta_x[j][i][0], &deta_median[j][i][0], &deta_x_err[j][i][0], &deta_median_err[j][i][0]);

    deta_median_gr[j][i]->SetMarkerStyle(24);
    deta_median_gr[j][i]->SetMarkerSize(0.3);
    TString nth_eg_sw_;
    nth_eg_sw_.Form("%d", j + 1); 
    TString eta_region;
    eta_region.Form("%d", i + 1); 
    deta_median_gr[j][i]->SetName( "Pix_Pix_deta_SW_"+nth_eg_sw_+"_eta_region_"+eta_region+"_median");

    deta_median_gr[j][i]->Write();
   }   
  }

   // save all 2D plots
   for(int j = 0; j < 10; j++){
      for(int i = 0; i < 5; i++){
         
         SADphi_dist[j][i]->Write();
      }
   }

   for(int j = 0; j < 4; j++){
      for(int i = 0; i < 5; i++){

         SADeta_dist[j][i]->Write();
      }
   }

                    
  file->Write();

}
