using namespace std;
using namespace mfem;

#if defined(Experiment9)
const std::string scenario = "testingExp9parallel";

const double alpha = 1; 
const int n_modes = 20;             // rational approximation 
const int N=3;                      // spectral method 
const double a = 0.5;               // weighted hermite polynomials 

// MESH 
const int dim = 2; 
const int n_refine = 0; 
const char *mesh_file = "../src/resources/meshes/donut.mesh";
// const char *mesh_file = "../src/resources/meshes/donut2d_0.msh";
    
// MESH 
void load_mesh(Mesh* mesh_){ 
    *mesh_ = Mesh(mesh_file); 
}
int boundary_D(int idx){
    return 1; 
}

// PARAMETERS  
//TODO: put xi and chi here 
// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 1.; 
const double scale_T = 1; 

// SIMULATION 
double t_final = 0.01;
double dt = 0.001; 
int n_plots = 20; 
int plot_frequency = 1; // int(t_final / dt / n_plots); 

// CALCULATE INITIAL PROBABILITY DENSITY 
// #define calculate_initial_condition
double t_final_IC = 1; 
double dt_IC = 0.01; 
int plot_frequency_IC = 10; 

// VELOCITY FIELD 
#define calculate_initial_velocity_field
// #define prescribed_velocity
double dt_IC_vel = 0.001; 
int n_iter_IC_vel = 1; 

void u_BC(const Vector &x, double t, Vector&u){

    u(0) = x(1);
    u(1) = -x(0);

    double inner_circle = 1; 
    double outer_circle = -1; 

    if(x(0) * x(0) + x(1) * x(1) < 1.5){
        u *= inner_circle * sin(0.5 * M_PI * t); 
    } else {
        u *= outer_circle * sin(0.5 * M_PI * t);    
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