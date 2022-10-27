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

using namespace std;
using namespace mfem;

int main(int argc, char *argv[]){    
    // Create a triangulation of the unitsquare (2 * n_x ^ 2 elements) 
    Mesh *mesh = new Mesh(); 
    *mesh = Mesh::MakeCartesian2D(n_x, n_x, Element::Type::TRIANGLE);

    // Define the finite element and the finite element space
    H1_FECollection fec(1, dim); 
    FiniteElementSpace fespace(mesh, &fec); // scalar 
    FiniteElementSpace vfespace(mesh, &fec, vector_size, Ordering::byNODES); // vector 
    const int n_dof = fespace.GetVSize(); 

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
    
    // Load initial conditions and fill phi^0 
    GridFunction phi(&vfespace, phi_block.GetData());  
    GridFunction phi0(&vfespace, phi0_block.GetData());  
    VectorFunctionCoefficient phi_initial(vector_size, phi0_function);  
    phi.ProjectCoefficient(phi_initial); 
    phi0.ProjectCoefficient(phi_initial); 

    // Create vector of blockvectors for modes 
    std::vector<BlockVector> phi_modes(n_modes,BlockVector(block_offsets)); 
    for (int i = 0; i < n_modes; i++){
        phi_modes[i] = 0.; 
    }

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
    // Output 

    // Create references to the single phis 
    std::vector<GridFunction> phis(vector_size);
    for (int i = 0; i < vector_size; i++ ){
        phis[i].MakeRef(&fespace, phi_block.GetBlock(i),0);
    }

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
    pd->SetLevelsOfDetail(1);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();
    
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

    // Configuration space solver 
    CSS css(fespace, vector_size, block_offsets);
    css.SetTime(t);
    ode_solver_css->Init(css);

    // Physical space solver  
    PSS pss(fespace, phi0_block, phi_modes); 
    pss.SetTime(t); 
    ode_solver_pss -> Init(pss); 

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
    // Time loop 
    
    for (int ti = 0; !done; ){

        cout << "t: " << t << "s / " << t_final << "s - dt: " << dt << endl;  
        
        // configuration space solver 
        ode_solver_css->Step(phi_block, t, dt);
        
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

        // int blocknumber = 42; 
        // for (int i = 0; i < phi_modes[0].GetBlock(blocknumber).Size(); i++){
        //     cout << phi_modes[0].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[1].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[2].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[3].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[4].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[5].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[6].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[7].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[2].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[8].GetBlock(blocknumber)[i] <<  " ";
        //     cout << phi_modes[9].GetBlock(blocknumber)[i] << endl; 
        // }

        // advance iteration counter and save output 
        ti++;
        pd->SetCycle(ti);
        pd->SetTime(t);
        pd->Save();

        done = (t >= t_final - 1e-8 * dt);
        // if(ti >= 20){
        //     return 0; 
        // }
    }
    
    cout << "t: " << t << "s / " << t_final << "s" << endl;  
    
    // Free the used memory.
    delete ode_solver_css;
    delete pd; 
    return 0;
}