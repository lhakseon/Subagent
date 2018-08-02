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

   TFile *file = new TFile("roi_median.root","recreate");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float low_et = 0.;
   float high_et = 0.;

   int binSize = 90;
   //float x[5][binSize], median[5][binSize], x_err[5][binSize], median_err[5][binSize];
   vector<float> x[8][5], median[8][5], x_err[8][5], median_err[8][5];
 
   for( int nth = 0; nth < 90; nth++){
      low_et = 10. + nth;
      high_et = low_et + 1.;
      cout << "et: " << low_et << endl;

      vector<float> dPhi_l1eg_pix[8][5];
      float abs_pterr;     
      float closest_egEt;  
      float closest_egEta; 
      float closest_egPhi; 
      float closest_eg_dr; 
 
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //for (Long64_t jentry=0; jentry<50000;jentry++) {
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

         int bpix_size = bRecHitGx->size();
         for(int a=0; a<bpix_size; a++){

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

        int n_pixels = 0; // initialize as 0 for each event
        hitted_layers.push_back(0); // 0 mean beam spot i.e., (0,0,0)
        for( int i=1; i < 9; i++){
           if( layers[i] != 0 ){
             hitted_layers.push_back(i); //check if hits on each barrel or disk exists
             n_pixels++;
           }
        }

       if( n_pixels >= 1 ){ // save if there is at least one pixel hit
            for( std::vector<int>::iterator first_hit = hitted_layers.begin()+1; first_hit != hitted_layers.end(); first_hit++){ // hitted_layers.begin()+1: to avoid using beam spot
               for( int k=0; k < layers[*first_hit]; k++){
                    double dPhi = 0.;
                    double R = 0.;

                    TVector3 pixel_vector;

                    if( *first_hit == 1 ) pixel_vector = first_layer_hits[k];
                    if( *first_hit == 2 ) pixel_vector = second_layer_hits[k];
                    if( *first_hit == 3 ) pixel_vector = third_layer_hits[k];
                    if( *first_hit == 4 ) pixel_vector = fourth_layer_hits[k];
                    if( *first_hit == 5 ) pixel_vector = first_disk_hits[k];
                    if( *first_hit == 6 ) pixel_vector = second_disk_hits[k];
                    if( *first_hit == 7 ) pixel_vector = third_disk_hits[k];
                    if( *first_hit == 8 ) pixel_vector = fourth_disk_hits[k];

                    dPhi = deltaPhi(pixel_vector.Phi(), closest_egPhi);

                    if( *first_hit == 1  ) {
                       if( fabs(closest_egEta) < 0.8 ) dPhi_l1eg_pix[0][0].push_back(dPhi); 
                       if( fabs(closest_egEta) >= 0.8 && fabs(closest_egEta) < 1.4 ) dPhi_l1eg_pix[0][1].push_back(dPhi); 
                       if( fabs(closest_egEta) >= 1.4 && fabs(closest_egEta) < 1.8 ) dPhi_l1eg_pix[0][2].push_back(dPhi); 
                       if( fabs(closest_egEta) >= 1.8 && fabs(closest_egEta) < 2.4 ) dPhi_l1eg_pix[0][3].push_back(dPhi); 
                       if( fabs(closest_egEta) >= 2.4 && fabs(closest_egEta) < 3.0 ) dPhi_l1eg_pix[0][4].push_back(dPhi); 
                    }
                    if( *first_hit == 2  ) { 
                       if( fabs(closest_egEta) < 0.8 ) dPhi_l1eg_pix[1][0].push_back(dPhi);
                       if( fabs(closest_egEta) >= 0.8 && fabs(closest_egEta) < 1.4 ) dPhi_l1eg_pix[1][1].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.4 && fabs(closest_egEta) < 1.8 ) dPhi_l1eg_pix[1][2].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.8 && fabs(closest_egEta) < 2.4 ) dPhi_l1eg_pix[1][3].push_back(dPhi);
                       if( fabs(closest_egEta) >= 2.4 && fabs(closest_egEta) < 3.0 ) dPhi_l1eg_pix[1][4].push_back(dPhi);
                    }
                    if( *first_hit == 3  ) { 
                        if( fabs(closest_egEta) < 0.8 ) dPhi_l1eg_pix[2][0].push_back(dPhi);
                       if( fabs(closest_egEta) >= 0.8 && fabs(closest_egEta) < 1.4 ) dPhi_l1eg_pix[2][1].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.4 && fabs(closest_egEta) < 1.8 ) dPhi_l1eg_pix[2][2].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.8 && fabs(closest_egEta) < 2.4 ) dPhi_l1eg_pix[2][3].push_back(dPhi);
                       if( fabs(closest_egEta) >= 2.4 && fabs(closest_egEta) < 3.0 ) dPhi_l1eg_pix[2][4].push_back(dPhi);
                    }
                    if( *first_hit == 4  ) { 
                       if( fabs(closest_egEta) < 0.8 ) dPhi_l1eg_pix[3][0].push_back(dPhi);
                       if( fabs(closest_egEta) >= 0.8 && fabs(closest_egEta) < 1.4 ) dPhi_l1eg_pix[3][1].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.4 && fabs(closest_egEta) < 1.8 ) dPhi_l1eg_pix[3][2].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.8 && fabs(closest_egEta) < 2.4 ) dPhi_l1eg_pix[3][3].push_back(dPhi);
                       if( fabs(closest_egEta) >= 2.4 && fabs(closest_egEta) < 3.0 ) dPhi_l1eg_pix[3][4].push_back(dPhi);
                    }
                    if( *first_hit == 5  ) { 
                       if( fabs(closest_egEta) < 0.8 ) dPhi_l1eg_pix[4][0].push_back(dPhi);
                       if( fabs(closest_egEta) >= 0.8 && fabs(closest_egEta) < 1.4 ) dPhi_l1eg_pix[4][1].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.4 && fabs(closest_egEta) < 1.8 ) dPhi_l1eg_pix[4][2].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.8 && fabs(closest_egEta) < 2.4 ) dPhi_l1eg_pix[4][3].push_back(dPhi);
                       if( fabs(closest_egEta) >= 2.4 && fabs(closest_egEta) < 3.0 ) dPhi_l1eg_pix[4][4].push_back(dPhi);
                    }
                    if( *first_hit == 6  ) { 
                       if( fabs(closest_egEta) < 0.8 ) dPhi_l1eg_pix[5][0].push_back(dPhi);
                       if( fabs(closest_egEta) >= 0.8 && fabs(closest_egEta) < 1.4 ) dPhi_l1eg_pix[5][1].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.4 && fabs(closest_egEta) < 1.8 ) dPhi_l1eg_pix[5][2].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.8 && fabs(closest_egEta) < 2.4 ) dPhi_l1eg_pix[5][3].push_back(dPhi);
                       if( fabs(closest_egEta) >= 2.4 && fabs(closest_egEta) < 3.0 ) dPhi_l1eg_pix[5][4].push_back(dPhi);
                    }
                    if( *first_hit == 7  ) { 
                       if( fabs(closest_egEta) < 0.8 ) dPhi_l1eg_pix[6][0].push_back(dPhi);
                       if( fabs(closest_egEta) >= 0.8 && fabs(closest_egEta) < 1.4 ) dPhi_l1eg_pix[6][1].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.4 && fabs(closest_egEta) < 1.8 ) dPhi_l1eg_pix[6][2].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.8 && fabs(closest_egEta) < 2.4 ) dPhi_l1eg_pix[6][3].push_back(dPhi);
                       if( fabs(closest_egEta) >= 2.4 && fabs(closest_egEta) < 3.0 ) dPhi_l1eg_pix[6][4].push_back(dPhi);
                    }
                    if( *first_hit == 8  ) { 
                       if( fabs(closest_egEta) < 0.8 ) dPhi_l1eg_pix[7][0].push_back(dPhi);
                       if( fabs(closest_egEta) >= 0.8 && fabs(closest_egEta) < 1.4 ) dPhi_l1eg_pix[7][1].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.4 && fabs(closest_egEta) < 1.8 ) dPhi_l1eg_pix[7][2].push_back(dPhi);
                       if( fabs(closest_egEta) >= 1.8 && fabs(closest_egEta) < 2.4 ) dPhi_l1eg_pix[7][3].push_back(dPhi);
                       if( fabs(closest_egEta) >= 2.4 && fabs(closest_egEta) < 3.0 ) dPhi_l1eg_pix[7][4].push_back(dPhi);
                    }
               }
             }
       }
      }// event loop
   
    for(int pixel = 0; pixel < 8; pixel++){ 
      for(int i = 0; i < 5; i++){ 
         std::sort (dPhi_l1eg_pix[pixel][i].begin(), dPhi_l1eg_pix[pixel][i].end());

         if( dPhi_l1eg_pix[pixel][i].size() != 0 ){
           x[pixel][i].push_back(low_et); median[pixel][i].push_back(getMedian(dPhi_l1eg_pix[pixel][i]));
           x_err[pixel][i].push_back(0.); median_err[pixel][i].push_back(getMedianErr(dPhi_l1eg_pix[pixel][i]));
         }
      }
   }
  }// et loop


  TGraphErrors* median_gr[8][5];
  for(int j = 0; j < 8; j++){
   for(int i = 0; i < 5; i++){
     int point_size = x[j][i].size();
     median_gr[j][i] = new TGraphErrors(point_size, &x[j][i][0], &median[j][i][0], &x_err[j][i][0], &median_err[j][i][0]);

    median_gr[j][i]->SetMarkerStyle(24);
    median_gr[j][i]->SetMarkerSize(0.3);
    TString pixel_;
    pixel_.Form("%d", j + 1);
    TString eta_region;
    eta_region.Form("%d", i + 1);
    median_gr[j][i]->SetTitle(pixel_+"th_Pixel_"+eta_region+"_median");
    median_gr[j][i]->SetName( "Pixel_"+pixel_+"eta_region_"+eta_region+"_median");

    median_gr[j][i]->Write();
   }
  }

  file->Write();
}
