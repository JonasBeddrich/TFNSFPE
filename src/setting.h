/*
    Test setting 
*/

using namespace std;
using namespace mfem;

const std::string scenario = "Test"; 
const int dim = 2;

// #define alpha_05
// #define alpha_07
#define alpha_09
// #define alpha_1

const int n_modes = 10; // number of modes for the rational approximation
const int N = 10; // degree of Hermite decomposition 
const int vector_size = N*N;

const double alpha = 0.5; 

double t_final = 2;
double dt = 0.05;
int n_x = 16; 

const double nu = 0.59; 
double xi = 1000; 
double chi = 1000; 
const double eps = 0.; 

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

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 16 * x(0) * x(0) * (1-x(0)) * (1-x(0)); 
    u(1) = 0; 
}

void u_IC(const Vector &x, double t, Vector &u){
   u(0) = 0; 
   u(1) = 0; 
}