using namespace std;
using namespace mfem;

#if defined(Experiment8)
const std::string scenario = "Exp8"; 

const double alpha = 1.0; 
const int n_modes = 20;             // rational approximation 
const int N=3;                      // spectral method 
const double a = 0.5;               // weighted hermite polynomials 

// MESH 
const int dim = 2; 
const int n_refine = 3; 

// MESH 
void load_mesh(Mesh* mesh){ 
    *mesh = Mesh::MakeCartesian2D(2, 2, Element::Type::QUADRILATERAL);
}

int boundary_D(int idx){
    if (idx == 1){
        return 1; 
    }
    return 0; 
}

// PARAMETERS  
//TODO: put xi and chi here 
// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 1.; 

// SIMULATION 
double t_final = 2.5;
double dt = 0.0001; 
int plot_frequency = 100; 

// CALCULATE INITIAL PROBABILITY DENSITY 
#define calculate_initial_condition
double t_final_IC = 1; 
double dt_IC = 0.01; 
int plot_frequency_IC = 10; 

// VELOCITY FIELD 
#define calculate_initial_velocity_field
// #define prescribed_velocity
double dt_IC_vel = 0.001; 
int n_iter_IC_vel = 2; 

void u_BC(const Vector &x, double t, Vector&u){
    
    u(1) = 0; 
    u(0) = x(1) * x(1) * (1-x(1)) * (1-x(1)); 
    
    if (t<0){
        u(0) = 0; 
    } else if (t < 2) {
        u(0) *= sin(0.5 * M_PI * t); 
    } else {
        u(0) = 0; 
    }
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
}
#endif 