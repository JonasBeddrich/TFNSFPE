#if defined(Experiment1)
const std::string scenario = "Exp1"; 
double t_final = 1;
const int n_refine = 0; 
bool prescribed_velocity = true; 

// xi and chi are set depening on psi ... thus not defined here 
const double eps = 1; 
const double nu = 0.; 

void u_BC(const Vector &x, double t, Vector &u){}; 
void u_IC(const Vector &x, double t, Vector &u){
    u(0) = x(1); 
    u(1) = 0; 
}
#endif 

// #if defined(Experiment2)
// const std::string scenario = "Exp2_uneven"; 
// double t_final = 60;
// const int n_refine = 0; 
// bool prescribed_velocity = true; 

// // xi and chi are set depening on psi ... thus not defined here 
// const double eps = 0.01; 
// const double nu = 0.; 

// void u_BC(const Vector &x, double t, Vector &u){}; 
// void u_IC(const Vector &x, double t, Vector &u){
//     u(0) = x(1) * (1-x(1)); 
//     u(1) = 0; 
// }
// #endif 

// #if defined(Experiment3)
// const std::string scenario = "testing"; 
// double t_final = 10;
// const int n_refine = 0; 
// bool prescribed_velocity = true; 

// // xi and chi are set depening on psi ... thus not defined here 
// const double eps = 0.01; 
// const double nu = 0.; 

// void u_BC(const Vector &x, double t, Vector &u){}; 
// void u_IC(const Vector &x, double t, Vector &u){
//     u(0) = 0.5 * x(0); 
//     u(1) = - 0.5 * x(1);
// }
// #endif 

// #if defined(Experiment4)
// const std::string scenario = "Exp4_no16_testing"; 
// double t_final = 10;
// const char *mesh_file = "../src/test.mesh";
// const int n_refine = 4; 
// bool prescribed_velocity = false;


// void load_mesh(Mesh* mesh){ 
//     *mesh = Mesh::MakeCartesian2D(2, 2, Element::Type::QUADRILATERAL);
// }

// // xi and chi are set depening on psi ... thus not defined here 
// const double nu = 0.59; 
// const double eps = 1.; // this is changed in comparison to the Mizerova and She paper 

// // void u_BC(const Vector &x, double t, Vector&u){
// //     if (x(1) > 0.5){ // top boundary 
// //         u(0) = 16 * x(0) * x(0) * (1-x(0)) * (1-x(0)); 
// //     } else { // bottom boundary 
// //         u(0) = 0; 
// //     }
// //     u(1) = 0; 
// // }

// void u_BC(const Vector &x, double t, Vector&u){
//     // if (x(1) > 0.999999){ // top boundary 
//     u(0) = x(0) * x(0) * (1-x(0)) * (1-x(0)); 
//     u(1) = - x(0) * x(1) * (4 * x(0) * x(0) - 6 * x(0) + 2); 
//     // } else { // other boundaries 
//     //     u(0) = 0; 
//     //     u(1) = 0; 
//     // }
// } 

// void u_IC(const Vector &x, double t, Vector &u){
//     // u(0) = x(0) * x(0) * (1-x(0)) * (1-x(0)); 
//     // u(1) = - x(0) * x(1) * (4 * x(0) * x(0) - 6 * x(0) + 2); 
//     u(0) = 0; 
//     u(1) = 0; 
// }
// #endif 

#if defined(Experiment4_Medea)
const std::string scenario = "Exp4_Medea_100"; 
double t_final = 10;
const char *mesh_file = "../src/test.mesh";
const int n_refine = 2; 
bool prescribed_velocity = false;

// xi and chi are set depening on psi ... thus not defined here 
const double nu = 1.; 
const double eps = 1.; 

void u_BC(const Vector &x, double t, Vector&u){
    if (x(1) > 0.999999){ // top boundary 
        u(0) = 0.1;
    } else { // other boundaries 
        u(0) = 0; 
    }
    u(1) = 0; 
} 

void u_IC(const Vector &x, double t, Vector &u){
    // u(0) = x(0) * x(0) * (1-x(0)) * (1-x(0)); 
    // u(1) = - x(0) * x(1) * (4 * x(0) * x(0) - 6 * x(0) + 2); 
    u(0) = 0; 
    u(1) = 0; 
}
#endif 

#if defined(Experiment5_pres_u)
const std::string scenario = "Exp5_pres_u"; 
double t_final = 2;
const int n_refine = 0;  
bool prescribed_velocity = true; 

// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.5; 
const double eps = 0.; 

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 0; 
    u(1) = 0; 
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = x(1) * (1-x(1)); 
    u(1) = 0; 
}
#endif 

#if defined(Experiment5_pres_C)
const std::string scenario = "Exp5_pres_C"; 
double t_final = 2;
const int n_refine = 0;  
bool prescribed_velocity = false; 

// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.5; 
const double eps = 0.; 

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 0; 
    u(1) = 0; 
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = x(1) * (1-x(1)); 
    u(1) = 0; 
}
#endif 

#if defined(Experiment6)
const std::string scenario = "Exp6"; 
double t_final = 4;
const char *mesh_file = "../src/square-disc.mesh";
const int n_refine = 0; 
bool prescribed_velocity = false; 

// xi and chi are set depening on psi ... thus not defined here 
const double nu = 0.59; 
const double eps = 1.; 

void u_BC(const Vector &x, double t, Vector&u){
    u(0) = 0; 
    u(1) = 0; 

    if(x(0) < 1e-8 || x(0) > 1 - 1e-8){
        u(0) = 0.25 * x(1) * (1-x(1)); 
    } 
}

void u_IC(const Vector &x, double t, Vector &u){
    u(0) = 0; 
    u(1) = 0; 
}

#endif 


// mesh = new Mesh(mesh_file);  

    // // *mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL);

    // #if defined(Experiment1)
    //     *mesh = Mesh::MakeCartesian2D(8, 8, Element::Type::QUADRILATERAL);
    // #endif

    // #if defined(Experiment2)
    //     *mesh = Mesh::MakeCartesian2D(8, 8, Element::Type::QUADRILATERAL);
    // #endif

    // #if defined(Experiment3)
    //     *mesh = Mesh::MakeCartesian2D(8, 8, Element::Type::QUADRILATERAL);
    // #endif

    // // #if defined(Experiment4)
    // //     *mesh = Mesh::MakeCartesian2D(2, 2, Element::Type::QUADRILATERAL);
    // // #endif

    // #if defined(Experiment4_Medea)
    //     *mesh = Mesh::MakeCartesian2D(2, 2, Element::Type::QUADRILATERAL);
    // #endif

    // #if defined(Experiment5_pres_u)
    //     Mesh square = Mesh::MakeCartesian2D(10, 10, Element::QUADRILATERAL);
    //     square.Save("Exp5_non_periodic.mesh");
    //     Vector left_right_translation({1.0, 0.0});
    //     std::vector<Vector> translations = {left_right_translation};
    //     *mesh = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
    //     mesh->RemoveInternalBoundaries();
    //     mesh->Save("Exp5_periodic.mesh");
    // #endif

    // #if defined(Experiment5_pres_C)
    //     Mesh square = Mesh::MakeCartesian2D(10, 10, Element::QUADRILATERAL);
    //     square.Save("Exp5_non_periodic.mesh");
    //     Vector left_right_translation({1.0, 0.0});
    //     std::vector<Vector> translations = {left_right_translation};
    //     *mesh = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
    //     mesh->RemoveInternalBoundaries();
    //     mesh->Save("Exp5_periodic.mesh");
    // #endif

    // #if defined(Experiment6)
    //     mesh = new Mesh(mesh_file);
    // #endif

    // #if defined(Experiment7)
    //     mesh = new Mesh(mesh_file);
    // #endif