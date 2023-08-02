using namespace std;
using namespace mfem;

#if defined(Experiment7)
const std::string scenario = "test"; 

const double alpha = 0.8; 
const int n_modes = 20;             // rational approximation 
const int N=3;                      // spectral method 
const double a = 0.5;               // weighted hermite polynomials 

// MESH 
const int dim = 2; 
const int n_refine = 0; 
const char *mesh_file = "../src/resources/meshes/half_channel_ref2.msh";

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

// SIMULATION 
double t_final = 1;
double dt = 0.001; 
int plot_frequency = 10; 

// CALCULATE INITIAL PROBABILITY DENSITY 
#define calculate_initial_condition
double t_final_IC = 1; 
double dt_IC = 0.01; 
int plot_frequency_IC = 10; 

// VELOCITY FIELD 
#define calculate_initial_velocity_field
// #define prescribed_velocity
double dt_IC_vel = 0.001; 
int n_iter_IC_vel = 100; 

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 0; 
    u(1) = 0; 

    if(x(0) < 1e-8){
        u(0) = 0.25 * x(1) * (1-x(1)); 
    } 
    if(x(0) > 3 - 1e-8){
        u(0) = 0.25 * x(1) * (1-x(1)); 
    } 

    if(t<0){
        u(0) *= 1; 
    }else{
        u(0) *= 2 - cos(2 * M_PI * t); 
    }
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
}
#endif 