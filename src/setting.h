// #define alpha_07
// #define alpha_08
// #define alpha_09
#define alpha_1

#define Experiment4
// #define Experiment5
// #define Experiment6

/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

using namespace std;
using namespace mfem;

const int dim = 2;
const int n_modes = 20; // number of modes for the rational approximation
const int N = 10; // degree of Hermite decomposition 
const int vector_size = N*N;
const double alpha = 0.5; 

/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

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

/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

#if defined(Experiment4)
const std::string scenario = "Exp4"; 
double t_final = 2;
const char *mesh_file = "../src/test.mesh";
const int n_refine = 3; 

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
const int n_refine = 0; 

// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 1.; 

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 0; 
    u(1) = 0; 

    if(x(0) < 1e-8 || x(0) > 1 - 1e-8){
        u(0) = 0.25 * x(1) * (1-x(1)); 
    } 
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
}

#endif 

/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

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
    y(10) = -0.0093679970435956;
    y(12) = 0.006726880736189749;
    y(14) = -0.004861638355237339;
    y(16) = 0.0035304458222095232;
    y(18) = -0.0025732324695826657;
    y(40) = -0.042202327319864355;
    y(42) = 0.02238116387229778;
    y(44) = -0.014536992359754286;
    y(46) = 0.009952798292147094;
    y(48) = -0.0069824927341656075;
    y(50) = 0.004968130676746484;
    y(52) = -0.003567467238594695;
    y(54) = 0.002578273086498702;
    y(56) = -0.0018723016361220524;
    y(58) = 0.0013646625966084816;
    y(80) = 0.027411215668371434;
    y(82) = -0.014536992359754286;
    y(84) = 0.009442053508625628;
    y(86) = -0.006464532119806319;
    y(88) = 0.004535262067145761;


    // y(90) = -0.0032268955315374235;
    // y(92) = 0.0023171379418437956;
    // y(94) = -0.001674637493101192;
    // y(96) = 0.0012160955853216496;
    // y(98) = -0.0008863743571930357;
    // y(120) = -0.018767176437757598;
    // y(122) = 0.009952798292147094;
    // y(124) = -0.006464532119806319;
    // y(126) = 0.004425962582168262;
    // y(128) = -0.0031050816729665523;
    // y(130) = 0.0022093043416695863;
    // y(132) = -0.0015864358994986752;
    // y(134) = 0.0011465459132692808;
    // y(136) = -0.0008326037302038944;
    // y(138) = 0.0006068590372859589;
    // y(160) = 0.0131663145651048;
    // y(162) = -0.0069824927341656075;
    // y(164) = 0.004535262067145761;
    // y(166) = -0.0031050816729665523;
    // y(168) = 0.0021784034584109418;
    // y(170) = -0.001549961233057453;
    // y(172) = 0.001112981175375484;
    // y(174) = -0.0008043716223111528;
    // y(176) = 0.0005841221057574263;
    // y(178) = -0.0004257484874234121;
    // y(200) = -0.0093679970435956;
    // y(202) = 0.004968130676746484;
    // y(204) = -0.0032268955315374235;
    // y(206) = 0.0022093043416695863;
    // y(208) = -0.001549961233057453;
    // y(210) = 0.001102816750820539;
    // y(212) = -0.0007918999890925135;
    // y(214) = 0.00057232044263429;
    // y(216) = -0.00041561016431562465;
    // y(218) = 0.0003029253593916452;
    // y(240) = 0.006726880736189749;
    // y(242) = -0.003567467238594695;
    // y(244) = 0.0023171379418437956;
    // y(246) = -0.0015864358994986752;
    // y(248) = 0.001112981175375484;
    // y(250) = -0.0007918999890925135;
    // y(252) = 0.0005686398871418406;
    // y(254) = -0.00041096632957588153;
    // y(256) = 0.0002984373281811338;
    // y(258) = -0.00021752171303129417;
    // y(280) = -0.004861638355237339;
    // y(282) = 0.002578273086498702;
    // y(284) = -0.001674637493101192;
    // y(286) = 0.0011465459132692808;
    // y(288) = -0.0008043716223111528;
    // y(290) = 0.00057232044263429;
    // y(292) = -0.00041096632957588153;
    // y(294) = 0.0002970127981946219;
    // y(296) = -0.0002156860539409193;
    // y(298) = 0.0001572068756148139;
    // y(320) = 0.0035304458222095232;
    // y(322) = -0.0018723016361220524;
    // y(324) = 0.0012160955853216496;
    // y(326) = -0.0008326037302038944;
    // y(328) = 0.0005841221057574263;
    // y(330) = -0.00041561016431562465;
    // y(332) = 0.0002984373281811338;
    // y(334) = -0.0002156860539409193;
    // y(336) = 0.00015662784279794556;
    // y(338) = -0.00011416117709352666;
    // y(360) = -0.0025732324695826657;
    // y(362) = 0.0013646625966084816;
    // y(364) = -0.0008863743571930357;
    // y(366) = 0.0006068590372859589;
    // y(368) = -0.0004257484874234121;
    // y(370) = 0.0003029253593916452;
    // y(372) = -0.00021752171303129417;
    // y(374) = 0.0001572068756148139;
    // y(376) = -0.00011416117709352666;
    // y(378) = 8.320854148640875e-05;
}