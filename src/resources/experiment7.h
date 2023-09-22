using namespace std;
using namespace mfem;

#if defined(Experiment7)
const std::string scenario = "Exp7_convergence_study"; 

const double alpha = 0.5; 
const int n_modes = 20;             // rational approximation 
const int N=3;                     // spectral method 
const double a = 0.5;               // weighted hermite polynomials 

// MESH 
const int dim = 2; 
const int n_refine = 1; // ACHTUNG NUR DER 5 REF SCHRITT HAT HIER EINE 1 !!! 
const char *mesh_file = "../src/resources/meshes/half_channel_ref4.msh";

// MESH 
void load_mesh(Mesh* mesh_){ 
    *mesh_ = Mesh(mesh_file); 
}

int boundary_D(int idx){
    if (idx == 4){
        return 0; 
    }
    return 1; 
    // attr[0] = 0; // Halbkreis  
    // attr[1] = 0; // unten links des Halbkreises
    // attr[2] = 0; // links 
    // attr[3] = 0; // oben 
    // attr[4] = 0; // rechts 
    // attr[5] = 0; // unten rechts des Halbkreises
}

// PARAMETERS  
//TODO: put xi and chi here 
// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 1.; 
const double scale_T = 5; 

// SIMULATION 
double t_final = 1;
// double dt = 0.00015625;
// double dt = 0.0000390625;
double dt = 0.000009765625; 
int plot_frequency = 1024; 

// CALCULATE INITIAL PROBABILITY DENSITY 
#define calculate_initial_condition
double t_final_IC = 1; 
double dt_IC = 0.01; 
int plot_frequency_IC = 10; 

// VELOCITY FIELD 
#define calculate_initial_velocity_field
// #define prescribed_velocity
double dt_IC_vel = 0.001; 
int n_iter_IC_vel = 10; 

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 0; 
    u(1) = 0; 

    if(x(0) < 1e-8){
        u(0) = x(1) * (1-x(1)); 
    } 
    if(x(0) > 3 - 1e-8){
        u(0) = x(1) * (1-x(1)); 
    } 

    if(t<0){
        u(0) *= 0; 
    } else if (t < 1){
        u(0) *= 1 - cos(2 * M_PI * t); 
    } else {
        u(0) *= 0; 
    }
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
}

void n_bdr(const Vector &x, Vector&n){
    // if (x[0] < 1e-8){ 
    //     // left 
    //     n[0] = -1; 
    //     n[1] = 0; 
    // } else if (x[0] > 3 - 1e-8){  
    //     // right 
    //     n[0] = 1; 
    //     n[1] = 0; 
    // } else if (x[1] < 1e-8){
    //     // bottom 
    //     n[0] = 0; 
    //     n[1] = -1; 
    // } else if (x[1] > 1 - 1e-8){
    //     // top 
    //     n[0] = 0; 
    //     n[1] = 1;
    // } else {
        n[0] = 0; 
        n[1] = 0;
    // }
}

#endif 
