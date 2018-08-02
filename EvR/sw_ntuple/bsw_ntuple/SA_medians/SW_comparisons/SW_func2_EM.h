#ifndef withEM_SingleCrys_h
#define withEM_SingleCrys_h

double SW_func2_dphi(int nth_par){
  double p[5];

   p[0] = 0.00460341;
   p[1] = 0.342664;
   p[2] = -0.0331638;
   p[3] = 0.352677;
   p[4] = 0.52688;
   return p[nth_par];

}

double SW_func2_deta(){

 double p[1] = {};

 p[0] = -1.67124e-05;

 return p[0];

}

#endif
