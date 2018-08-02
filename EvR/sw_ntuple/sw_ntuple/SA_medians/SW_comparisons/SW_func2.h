#ifndef withoutEM_SingleCrys_2_h
#define withoutEM_SingleCrys_2_h

double SW_func2_dphi( int nth, int nth_par){
  double p[4];

 if( nth == 0){
p[0] = -0.000134002;
p[1] = 0.0451102;
p[2] = -0.479288;
p[3] = 0.179786;
return p[nth_par];
 }
 
 if( nth == 1){
p[0] = -0.000140397;
p[1] = 0.0713432;
p[2] = -0.415087;
p[3] = 0.207985;
return p[nth_par];
 }
 
 if( nth == 2){
p[0] = 0.000214606;
p[1] = 0.0727317;
p[2] = -0.145126;
p[3] = 0.303335;
return p[nth_par];
 }
 
 if( nth == 3){
p[0] = 9.79928e-05;
p[1] = 0.0656966;
p[2] = -0.343087;
p[3] = 0.247032;
return p[nth_par];
 }
 
 if( nth == 4){
p[0] = 0.000577854;
p[1] = 0.0665647;
p[2] = -0.0685691;
p[3] = 0.335524;
return p[nth_par];
 }

 if( nth == 5){
p[0] = 0.000521113;
p[1] = 0.0684741;
p[2] = -0.0659904;
p[3] = 0.338502;
return p[nth_par];
 }
 
 if( nth == 6){
p[0] = 0.000347119;
p[1] = 0.0422779;
p[2] = -0.0955098;
p[3] = 0.342277;
return p[nth_par];
 }
 
 if( nth == 7){
p[0] = -7.06036e-05;
p[1] = 0.0452167;
p[2] = -0.0410123;
p[3] = 0.297338;
return p[nth_par];
 }
 
 if( nth == 8){
p[0] = -0.00041432;
p[1] = 0.0435384;
p[2] = -0.0539955;
p[3] = 0.275643;
return p[nth_par];
 }
 
 if( nth == 9){
p[0] = 0.000309466;
p[1] = 0.0384042;
p[2] = -0.0855676;
p[3] = 0.33162;
return p[nth_par];
 }

}

double SW_func2_deta(){

 double p[0] = {};

 p[0] = 1.81345e-07;
 
 return p[0];
 
}

#endif
