// #define alpha_07
// #define alpha_08
// #define alpha_09
#define alpha_1

// #define Experiment4
// #define Experiment5
#define Experiment6

using namespace std;
using namespace mfem;

const int dim = 2;
const int n_modes = 20; // number of modes for the rational approximation
const int N = 10; // degree of Hermite decomposition 
const int vector_size = N*N;
const double alpha = 0.5; 

void u(const Vector &x, Vector &u){
    u(0) = x(1) * (1 - x(1)); 
    u(1) = 0; 
}

double du1dx1(const Vector &x){
    return 0.; 
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

// this is not nice but it avoids a lot of implicit type casting 
const double zero(const Vector &x){
    return 0; 
}

#if defined(alpha_1)
// double dt = 0.05; 
/* 
The original values produces a visible error within the first time steps. 
Probably this is due to the initial conditions (not given in the paper). 
*/
double dt = 0.01;
int plot_frequency = 1; 
#endif 

#if defined(alpha_09)
double dt = 0.005;
int plot_frequency = 2; 
#endif 

#if defined(alpha_08)
double dt = 0.001;
int plot_frequency = 10; 
#endif 

#if defined(alpha_07)
double dt = 0.0001;
int plot_frequency = 100; 
#endif 

#if defined(Experiment4)
const std::string scenario = "Exp4"; 
double t_final = 2;
const char *mesh_file = "../src/test.mesh";

// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 0.; 

void u_BC(const Vector &x, double t, Vector&u){
    if (x(1) > 0.5){ // top boundary 
        u(0) = 16 * x(0) * x(0) * (1-x(0)) * (1-x(0)); 
    } else { // bottom boundary 
        u(0) = 0; 
    }
    u(1) = 0; 
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
}
#endif 

#if defined(Experiment5)

#endif 

#if defined(Experiment6)
const std::string scenario = "Exp6"; 
double t_final = 4;
const char *mesh_file = "../src/square-disc.mesh";

// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 1.; 

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 0; 
    u(1) = 0; 
    
    if(x(0) < 1e-8){
        u(0) = 0.25 * x(1) * (1-x(1)); 
    } 
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
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