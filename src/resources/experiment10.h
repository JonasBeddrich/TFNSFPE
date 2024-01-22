using namespace std;
using namespace mfem;

#if defined(Experiment10)
const std::string scenario = "Exp10NSFP_Re1000";

const double alpha = 0.7; 
const int n_modes = 20;             // rational approximation 
const int N=3;                      // spectral method 
const double a = 0.5;               // weighted hermite polynomials 

// MESH 
const int dim = 2; 
const int n_refine = 0; 
const char *mesh_file = "../src/resources/meshes/vortex_street.msh";
    
// MESH 
void load_mesh(Mesh* mesh_){ 
    *mesh_ = Mesh(mesh_file); 
}

int boundary_D(int idx){
    if (idx == 3){
        return 0; 
    }
    return 1; 
    cout << idx << endl; 
}

// PARAMETERS  
//TODO: put xi and chi here 
// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.0001; 
const double eps = 1.; 
const double scale_T = 1; 

// SIMULATION 
double t_final = 10;
double dt = 0.0001; 
int n_plots = 1000; 
int plot_frequency = int(t_final / dt / n_plots); 

// CALCULATE INITIAL PROBABILITY DENSITY 
// #define calculate_initial_condition
double t_final_IC = 1; 
double dt_IC = 0.01; 
int plot_frequency_IC = 10; 

// VELOCITY FIELD 
// #define calculate_initial_velocity_field
// #define prescribed_velocity
double dt_IC_vel = 0.0001; 
int n_iter_IC_vel = 1; 

void u_BC(const Vector &x, double t, Vector&u){
    
    if(x(0) < 1e-8){
        u(0) = 6 * x(1) * (0.41-x(1)) / (0.41 * 0.41); 
        u(1) = 0 ; 
    } else if(x(0) > 2.2 - 1e-8){
        
    }
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
}


void n_bdr(const Vector &x, Vector&n){
    n[0] = 0; 
    n[1] = 0;
}


#endif 