using namespace std;
using namespace mfem;

#if defined(Experiment2)

#define prescribed_velocity

const std::string scenario = "InflVel_non_diagonal"; 

const double alpha = 0.5; 
const int n_modes = 20;             // rational approximation 
const int N=20;                     // spectral method 
const double a = 0.5;               // weighted hermite polynomials 

// MESH 
const int dim = 2; 
const int n_refine = 1; 

void load_mesh(Mesh* mesh){
    *mesh = Mesh::MakeCartesian2D(8, 8, Element::Type::QUADRILATERAL);
}

// PARAMETERS  
// TODO: put xi and chi here 
// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.; 
const double eps = 0.01; 

// SIMULATION 
double t_final = 1.;
double dt = 0.001; 
int n_plots = 100; 
int plot_frequency = 100; // int(t_final / dt / n_plots); 

// CALCULATE INITIAL PROBABILITY DENSITY 
// #define calculate_initial_condition
double t_final_IC = 10; 
double dt_IC = 0.01; 
int plot_frequency_IC = 100; 

// VELOCITY FIELD 
double dt_IC_vel = 0.001; 
int n_iter_IC_vel = 1000; 

void u_BC(const Vector &x, double t, Vector&u){}; 

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0.5 * (x[0] + x[1]); 
    u(1) = - 0.5 * (x[0] + x[1]); 
}

void n_bdr(const Vector &x, Vector&n){
    n[0] = 0; 
    n[1] = 0;
    
    if (x[0] < 1e-8){ 
        // left 
        n[0] = -1; 
        n[1] = 0; 
    } else if (x[0] > 1 - 1e-8){  
        // right 
        n[0] = 1; 
        n[1] = 0; 
    } else if (x[1] < 1e-8){
        // bottom 
        n[0] = 0; 
        n[1] = -1; 
    } else if (x[1] > 1 - 1e-8){
        // top 
        n[0] = 0; 
        n[1] = 1;
    } 
}

#endif 
