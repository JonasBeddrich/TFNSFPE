using namespace std;
using namespace mfem;

#if defined(Experiment4)
const std::string scenario = "Exp4_rev"; 

const double alpha = 1.0; 
const int n_modes = 20;             // rational approximation 
const int N=10;                      // spectral method 
const double a = 0.5;               // weighted hermite polynomials 

// MESH 
const int dim = 2; 
const int n_refine = 3; 

void load_mesh(Mesh* mesh){ 
    *mesh = Mesh::MakeCartesian2D(2, 2, Element::Type::QUADRILATERAL);
}



// PARAMETERS  
// TODO: put xi and chi here 
// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 1.; 

// SIMULATION 
double t_final = 0.1;
double dt = 0.0001; 
int plot_frequency = 1000; 

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