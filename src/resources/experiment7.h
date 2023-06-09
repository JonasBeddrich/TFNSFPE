#if defined(Experiment7)

const std::string scenario = "Exp7"; 

// SIMULATION 
double t_final = 5;
double dt = 0.01; 
int plot_frequency = 10; 

// CALCULATE INITIAL VELOCITY FIELD 
#define calculate_initial_velocity_field
// #define prescribed_velocity
double dt_IC_vel = 0.001; 
int n_iter_IC_vel = 100; 

// CALCULATE INITIAL PROBABILITY DENSITY 
#define calculate_initial_condition
double t_final_IC = 1; 
double dt_IC = 0.001; 
int plot_frequency_IC = 100; 

// MESH 
const int dim = 2; 
const char *mesh_file = "../src/resources/meshes/half_channel_ref1.msh";
const int n_refine = 0; 


// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 1.; 

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 0; 
    u(1) = 0; 

    if(x(0) < 1e-8 || x(0) > 3 - 1e-8){
        u(0) = 0.25 * x(1) * (1-x(1)); 
    } 
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
}

#endif 