using namespace std;
using namespace mfem;

#if defined(Experiment3)
const std::string scenario = "Exp3la"; 

const double alpha = 0.5; 
const int n_modes = 20;             // rational approximation 
const int N=10;                     // spectral method 
const double a = 0.5;               // weighted hermite polynomials 

// MESH 
const int dim = 2; 
const int n_refine = 0; 

void load_mesh(Mesh* mesh){ 
    *mesh = Mesh::MakeCartesian2D(4, 4, Element::Type::QUADRILATERAL);
}

// PARAMETERS  
//TODO: put xi and chi here 
// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0; 
const double eps = 0.001; 

// SIMULATION 
double t_final = 10;
double dt = 0.1; 
int n_plots = 100; 
int plot_frequency = int(t_final / dt / n_plots); 

// CALCULATE INITIAL PROBABILITY DENSITY 
// #define calculate_initial_condition
double t_final_IC = 10; 
double dt_IC = 0.01; 
int plot_frequency_IC = 10; 

// VELOCITY FIELD 
// #define calculate_initial_velocity_field
#define prescribed_velocity
double dt_IC_vel = 0.001; 
int n_iter_IC_vel = 1000; 

void u_BC(const Vector &x, double t, Vector&u){}; 

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0.5 * x(0); 
    u(1) = -0.5 * x(1);
}
#endif 