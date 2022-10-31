#include "mfem.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional> 

#include <setting.h>
#include <rational_approximation.h> 
#include <finite_difference_scheme.h>
#include <PhysicalSpaceSolver.h>
#include <ConfigurationSpaceSolver.h>

#include <../../mfem-4.5/miniapps/navier/navier_solver.hpp>
#include <../../mfem-4.5/miniapps/navier/navier_solver.cpp>

using namespace std;
using namespace mfem;
using namespace navier; 

int main(int argc, char *argv[]){    
    Mpi::Init(); 
    Hypre::Init(); 

    // Create a triangulation of the unitsquare (2 * n_x ^ 2 elements) 
    Mesh *mesh = new Mesh(); 
    *mesh = Mesh::MakeCartesian2D(n_x, n_x, Element::Type::TRIANGLE);
    mesh->UniformRefinement(); 

    // Define the finite element and the finite element space
    H1_FECollection fec(1, dim); 
    FiniteElementSpace fespace(mesh, &fec); // scalar 
    FiniteElementSpace vfespace(mesh, &fec, vector_size, Ordering::byNODES); // vector 
    const int n_dof = fespace.GetVSize(); 

    // setup navier stokes solver ... needed earlier than the others  ... 
    auto *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh); 
    // delete mesh; 
    NavierSolver flowsolver(pmesh, 2, 1.0 / 100000.0);

    // Create offset vector for phi vector
    // i.e. indices where the next phi_ij starts  
    Array<int> block_offsets(vector_size+1); 
    block_offsets[0]=0; 
    for (int i = 1; i < vector_size + 1; i++){
        block_offsets[i] =  n_dof; 
    }
    block_offsets.PartialSum(); 

    // Create blockvetor phi and phi^0 
    BlockVector phi_block(block_offsets); 
    BlockVector phi0_block(block_offsets);
    BlockVector F_R_block(block_offsets); 
    BlockVector F_x_block(block_offsets); 

    phi_block = 0.; 
    phi0_block = 0.; 
    F_R_block = 0.; 
    F_x_block = 0.; 
    
    // velocity initial conditions 
    ParGridFunction *u_ic = flowsolver.GetCurrentVelocity(); 
    VectorFunctionCoefficient u_IC_coeff(pmesh->Dimension(), u_IC); 
    u_ic->ProjectCoefficient(u_IC_coeff); 

    // Load initial conditions and fill phi^0 
    GridFunction phi(&vfespace, phi_block.GetData());  
    GridFunction phi0(&vfespace, phi0_block.GetData());  
    VectorFunctionCoefficient phi_initial(vector_size, phi0_function);  
    phi.ProjectCoefficient(phi_initial); 
    phi0.ProjectCoefficient(phi_initial); 

     // Create references to the single phis 
    std::vector<GridFunction> phis(vector_size);
    for (int i = 0; i < vector_size; i++ ){
        phis[i].MakeRef(&fespace, phi_block.GetBlock(i),0);
    }

    // Create vector of blockvectors for modes 
    std::vector<BlockVector> phi_modes(n_modes,BlockVector(block_offsets)); 
    for (int i = 0; i < n_modes; i++){
        phi_modes[i] = 0.; 
    }

    // velocity boundary conditions  
    Array<int> attr(pmesh->bdr_attributes.Max());
    attr = 1;
    flowsolver.AddVelDirichletBC(u_BC, attr); 

    // ****************************************************************
    // Rational Approximation

    std::vector<double> lambdas = get_lambdas(); 
    std::vector<double> weights = get_weights(); 
    double w_inf = get_w_infinity(); 

    std::vector<double> gammas = get_gammas(dt); 
    double beta = get_beta(dt); 
    double delta = get_delta(dt); 

    for (auto i :weights){
        cout << i << " "; 
    }
    cout << endl; 

    for (auto i :gammas){
        cout << i << " "; 
    }
    cout << endl; 
   
    // ****************************************************************
    // Setup of the simulation 

    double t = 0.0;
    double tmp = 0.0; 
    bool done = false;

    Vector tmp_vector(n_dof); 
    BlockVector tmp_block_vector(block_offsets); 

    // Define two ODE solvers 
    ODESolver *ode_solver_css = NULL;
    ODESolver *ode_solver_pss = NULL;
    ode_solver_css = new BackwardEulerSolver;
    ode_solver_pss = new BackwardEulerSolver;

    // creating a dummy u to base the FD scheme from a gridfunction 
    FiniteElementSpace v2dfespace(mesh, &fec, 2, Ordering::byNODES); // vector 
    cout << "ready to load u" << endl;  
    ParGridFunction *u_gf_NS = flowsolver.GetCurrentVelocity(); 

    GridFunction u_gf(&v2dfespace); 
    GridFunction u1_gf(&fespace); 
    GridFunction u2_gf(&fespace); 

    Vector e1(2); 
    Vector e2(2);
    e1[0] = 1.;  
    e2[1] = 1.; 
    VectorConstantCoefficient e1_coeff(e1); 
    VectorConstantCoefficient e2_coeff(e2); 

    VectorGridFunctionCoefficient u_coeff(u_gf_NS);
    InnerProductCoefficient u1_coeff(u_coeff, e1_coeff); 
    InnerProductCoefficient u2_coeff(u_coeff, e2_coeff); 

    u1_gf.ProjectCoefficient(u1_coeff); 
    u2_gf.ProjectCoefficient(u2_coeff); 
    
    // Configuration space solver 
    CSS css(fespace, vector_size, block_offsets, u1_gf, u2_gf);
    css.SetTime(t);
    ode_solver_css->Init(css);

    // Physical space solver  
    PSS pss(fespace, phi0_block, phi_modes); 
    pss.SetTime(t); 
    ode_solver_pss -> Init(pss); 

    // Navier Stokes 
    flowsolver.Setup(dt); 

    BilinearForm m(&fespace);
    m.AddDomainIntegrator(new MassIntegrator); 
    m.Assemble(); 

    BiCGSTABSolver m_solver; 
    m_solver.iterative_mode=false; 
    m_solver.SetRelTol(1e-12);
    m_solver.SetMaxIter(1000);
    m_solver.SetPrintLevel(0);
    m_solver.SetOperator(m);

    // ****************************************************************
    // Output 

    // Prepare paraview output and save initial conditions  
    ParaViewDataCollection *pd = NULL;
    pd = new ParaViewDataCollection(scenario 
        + "_" + tf_degree
        + "_N=" + to_string(N) 
        + "_m=" + to_string(n_modes)
        + "_nx=" + to_string(n_x) 
        + "_dt=" + to_string(dt), mesh);
    pd->SetPrefixPath("ParaView");
    for (int i = 0; i < vector_size; i++ ){
        pd->RegisterField("phi " + std::to_string(i), &phis[i]);
    } 
    pd->RegisterField("velocity", u_gf_NS); 
    pd->SetLevelsOfDetail(1);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();

    // ****************************************************************
    // Time loop 
    
    for (int ti = 0; !done; ){

        cout << "t: " << t << "s / " << t_final << "s - dt: " << dt << endl;  
        
        // configuration space solver 
        tmp = t; 
        ode_solver_css->Step(phi_block, t, dt);
        t = tmp; 

        // calculate FR 
        css.Mult(phi_block, F_R_block); 

        // calculate M^-1 FR 
        for(int i= 0; i < vector_size; i++){
            m_solver.Mult(F_R_block.GetBlock(i), tmp_vector);
            F_R_block.GetBlock(i) = tmp_vector; 
        }

        // update modes for the configuration space 
        for (int l = 0; l < n_modes; l++){
            F_R_block *= dt * weights[l]; 
            phi_modes[l] += F_R_block; 
        }

        // physical space solver 
        tmp = t; 
        for(int i = 0; i < vector_size; i++){
            pss.set_current_block(i); 
            ode_solver_pss -> Step(phi_block.GetBlock(i),t,dt); 
        }   
        t = tmp; 

        // calculate Fx 
        for(int i= 0; i < vector_size; i++){
            pss.Mult(phi_block.GetBlock(i), tmp_vector); 
            m_solver.Mult(tmp_vector, F_x_block.GetBlock(i)); 
        }
        
        // update modes for the physical space 
        for (int l = 0; l < n_modes; l++){
            phi_modes[l].Add(dt * weights[l], F_x_block); 
            phi_modes[l] *= gammas[l]; 
        }

        flowsolver.Step(t, dt, ti); 

        // advance iteration counter and save output 
        ti++;
        pd->SetCycle(ti);
        pd->SetTime(t);
        pd->Save();

        done = (t >= t_final - 1e-8 * dt);
    }
    
    cout << "t: " << t << "s / " << t_final << "s" << endl;  
    
    // Free the used memory.
    delete ode_solver_css;
    delete pd; 
    return 0;
}

