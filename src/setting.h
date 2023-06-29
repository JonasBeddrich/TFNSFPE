using namespace std;
using namespace mfem;

// #define Experiment1
// #define Experiment2
// #define Experiment3
#define Experiment4
// #define Experiment4_Medea
// #define Experiment5_pres_u
// #define Experiment5_pres_C
// #define Experiment6
// #define Experiment7

/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

#define normal 
// #define symmetric
// #define skew_symmetric

#if defined(normal) 
    const std::string symmetry = ""; 
#endif 

#if defined(symmetric) 
    const std::string symmetry = "_sym"; 
#endif 

#if defined(skew_symmetric) 
    const std::string symmetry = "_skew"; 
#endif 

/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
// for some tests

void div_T(const Vector &x, double t, Vector &u){
    u(0) = -10; 
    u(1) = -10; 
} 