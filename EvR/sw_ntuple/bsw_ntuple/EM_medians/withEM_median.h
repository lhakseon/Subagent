#ifndef withEM_median_h
#define withEM_median_h


double SW_func2_dphi_v2(int region, double eget){
 double p[5] = {};

 if(region == 1){
   p[0] = 0.00460341;
   p[1] = 0.342664;
   p[2] = -0.0331638;
   p[3] = 0.352677;
   p[4] = 0.52688;
 }

 if(region == 2){
   p[0] = 0.00433023;
   p[1] = 0.345418;
   p[2] = -0.0812881;
   p[3] = 0.340587;
   p[4] = 0.537429;
 }

 if(region == 3){
   p[0] = 0.000144747;
   p[1] = 0.28076;
   p[2] = -0.0873595;
   p[3] = 0.284821;
   p[4] = 0.368885;
 }

 if(region == 4){
   p[0] = -0.00142974;
   p[1] = 0.288398;
   p[2] = -0.818192;
   p[3] = -0.0417389;
   p[4] = 0.709082;
 }

 if(region == 5){
   p[0] = -0.00408571;
   p[1] = 0.145543;
   p[2] = -0.172213;
   p[3] = 0.193732;
   p[4] = 0.0531367;
 }

 return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
}

double SW_func2_deta_v2(double eget){

 double p[1] = {};

 p[0] = -1.67124e-05;

 return p[0]*pow(eget,0);

}



#endif
