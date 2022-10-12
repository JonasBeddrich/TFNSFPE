using namespace std;
using namespace mfem;

const double alpha = 0.5; 

// this avoids a lot of implicit type casting 
const double zero(const Vector &x){
  return 0; 
}


// #define Experiment1 
// #define Experiment2 
#define Experiment3 

#if defined(Experiment1)
  const std::string scenario = "Exp1"; 
  const int N = 21; 
  double t_final = 10;

  double dt = 0.05;
  int n_x = 10; 

  const double xi = 1; 
  const double chi = 1; 
  const double eps = 1; 

  void u(const Vector &x, Vector &u){
    u(0) = x(1); 
    u(1) = 0; 
  }

  double du1dx1(const Vector &x){
    return 0; 
  }

  // Here this is not correct. But it is needed to reproduce the results ... a little bit odd. 
  // Not sure who made a mistake here. 
  double du1dx2(const Vector &x){
    return 1;  
  }

  double du2dx1(const Vector &x){
    return 0;  
  }    

  double du2dx2(const Vector &x){
    return 0; 
  }
#endif

#if defined(Experiment2) 
  const std::string scenario = "Exp2"; 
  const int N = 21; 
  double t_final = 10;

  double dt = 0.1;
  int n_x = 16; 

  const double xi = 1; 
  const double chi = 1; 
  const double eps = 0.01; 

  void u(const Vector &x, Vector &u){
    u(0) = x(1) * (1 - x(1)); 
    u(1) = 0; 
  }

  double du1dx1(const Vector &x){
    return 0; 
  }

  double du1dx2(const Vector &x){
    return 1 - 2 * x(1);  
  }

  double du2dx1(const Vector &x){
    return 0;  
  }    

  double du2dx2(const Vector &x){
    return 0; 
  }
#endif

#if defined(Experiment3)
  const std::string scenario = "Exp3";
  const int N = 10; 
  double t_final = 15; // reasonable close to steady state

  double dt = 0.05;
  int n_x = 4; 

  const double kappa = 0.5; 

  const double xi = 1; 
  const double chi = 1; 

  // this is not specified in the paper ... 
  const double eps = 0.01; 

  void u(const Vector &x, Vector &u){
    u(0) = kappa * x(0); 
    u(1) = - kappa * x(1); 
  }

  double du1dx1(const Vector &x){
    return kappa; 
  }

  double du1dx2(const Vector &x){
    return 0;  
  }

  double du2dx1(const Vector &x){
    return 0;  
  }    

  double du2dx2(const Vector &x){
    return - kappa; 
  } 
#endif 

// Initial condition
void phi0_function(const Vector &x, Vector &y){
  int dim = x.Size(); 
  int vector_size = y.Size(); 

  for (int i = 0; i < vector_size; i++){
    y(i) = 0; 
  }

  y(0) = 0.07957747154594769;
  y(2) = -0.042202327319864355;
  y(4) = 0.027411215668371434;
  y(6) = -0.018767176437757598;
  y(8) = 0.0131663145651048;
  y(20) = -0.042202327319864355;
  y(22) = 0.02238116387229778;
  y(24) = -0.014536992359754286;
  y(26) = 0.009952798292147094;
  y(28) = -0.0069824927341656075;
  y(40) = 0.027411215668371434;
  y(42) = -0.014536992359754286;
  y(44) = 0.009442053508625628;
  y(46) = -0.006464532119806319;
  y(48) = 0.004535262067145761;
  y(60) = -0.018767176437757598;
  y(62) = 0.009952798292147094;
  y(64) = -0.006464532119806319;
  y(66) = 0.004425962582168262;
  y(68) = -0.0031050816729665523;
  y(80) = 0.0131663145651048;
  y(82) = -0.0069824927341656075;
  y(84) = 0.004535262067145761;
  y(86) = -0.0031050816729665523;
  y(88) = 0.0021784034584109418;
}
