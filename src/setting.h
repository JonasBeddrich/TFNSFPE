/*
    Test setting 
*/

using namespace std;
using namespace mfem;

const std::string scenario = "Test"; 
const int dim = 2;

const int m = 1; // number of modes for the rational approximation
const int N = 10; // degree of Hermite decomposition 
const int vector_size = N*N;

const double alpha = 0.5; 

double t_final = 10;
double dt = 0.1;
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

double du1dx2(const Vector &x){
    return 1;  
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

// Rational Approximation 

std::vector<double> get_lambdas(){
    std::vector<double> lambdas(m,0.); 
    lambdas[0] = 0; 
    return lambdas; 
}

std::vector<double> get_weights(){
    std::vector<double> weights(m,0.); 
    weights[0] = 1; 
    return weights; 
}

double get_w_infinity(){
    return 0.; 
}

std::vector<double> get_gammas(double dt){
    std::vector<double> gammas(get_lambdas()); 
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::multiplies<double>(), std::placeholders::_1, dt));
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::plus<double>(), std::placeholders::_1, 1));
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::divides<double>(), 1, std::placeholders::_1));
    return gammas; 
} 

double get_beta(double dt){
    double beta = 0; 
    std::vector<double> gammas = get_gammas(dt); 
    std::vector<double> weights = get_weights(); 
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::multiplies<double>(), std::placeholders::_1, dt));
    std::transform(gammas.begin(), gammas.end(), weights.begin(), gammas.begin(), std::multiplies<double>()); 
    std::for_each(gammas.begin(), gammas.end(), [&] (double d) {beta += d;}); 
    return beta + get_w_infinity(); 
}

double get_delta(double dt){
    double delta = 0; 
    std::vector<double> weights = get_weights(); 
    std::for_each(weights.begin(), weights.end(), [&] (double d) {delta += d;}); 
    return delta; 
}