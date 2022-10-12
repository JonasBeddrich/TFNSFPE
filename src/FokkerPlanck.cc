#include "mfem.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional> 

#include <experiment_settings.h>
#include <finite_difference_scheme.h>
#include <PhysicalSpaceSolver.h>
#include <ConfigurationSpaceSolver.h>

using namespace std;
using namespace mfem;

// Size of the vector phi 
const int vector_size = N*N;
const int dim = 2;  

int main(int argc, char *argv[]){    
    // Create a triangulation of the unitsquare (2 * n_x ^ 2 elements) 
    Mesh *mesh = new Mesh(); 
    *mesh = Mesh::MakeCartesian2D(n_x, n_x, Element::Type::TRIANGLE);

    // Define the finite element and the finite element space
    H1_FECollection fec(1, dim); 
    FiniteElementSpace fespace(mesh, &fec); // scalar 
    FiniteElementSpace vfespace(mesh, &fec, vector_size, Ordering::byNODES); // vector 

    // Create offset vector for phi vector
    // i.e. indices where the next phi_ij starts  
    Array<int> block_offsets(vector_size+1); 
    block_offsets[0]=0; 
    for (int i = 1; i < vector_size + 1; i++){
        block_offsets[i] =  fespace.GetVSize(); 
    }
    block_offsets.PartialSum(); 

    // Create blockvetor for all phis 
    BlockVector phi_block(block_offsets); 
    phi_block = 1.0; 

    // Initial conditions
    GridFunction phi(&vfespace, phi_block.GetData());  
    VectorFunctionCoefficient phi0(vector_size, phi0_function);  
    phi.ProjectCoefficient(phi0); 

    // ****************************************************************
    // Configuration Space Solver 

    // Create CoefficientFactory objects for the finite difference scheme 
    std::vector<CoefficientFactory> CoefficientVector(vector_size, CoefficientFactory()); 
    fill_CoefficientVector(CoefficientVector); 

    // Create matrix A and A_entries 
    std::vector<std::vector<FunctionCoefficient>> A(vector_size, std::vector<FunctionCoefficient>(vector_size, FunctionCoefficient(zero))); 
    std::vector<std::vector<bool>> A_entries(vector_size, std::vector<bool>(vector_size)); 
    // Fill matrix A and A_entries 
    fill_A(A, A_entries, CoefficientVector); 
  
    // Create blockoprators 
    BlockOperator MmdtA_BO(block_offsets); 
    BlockOperator A_BO(block_offsets); 
    
    // The sparse matrices will be stored here 
    std::vector<std::vector<SparseMatrix>> MmdtA_SpMat(vector_size, std::vector<SparseMatrix>(vector_size)); 
    std::vector<std::vector<SparseMatrix>> A_SpMat(vector_size, std::vector<SparseMatrix>(vector_size)); 

    // Create mass matrix 
    BilinearForm m(&fespace);
    m.AddDomainIntegrator(new MassIntegrator); 
    m.Assemble(); 

    // Create zero matrix 
    BilinearForm m0(&fespace); 
    ConstantCoefficient m0_coef(0.);
    m0.AddDomainIntegrator(new MassIntegrator(m0_coef)); 
    m0.Assemble(); 

    // Fill the block operators with A and (Id - dt A) 
    for(int i = 0; i < vector_size; i++){
        for(int j = 0; j < vector_size; j++){
            if(A_entries[i][j]){
                BilinearForm a(&fespace);  
                a.AddDomainIntegrator(new MassIntegrator(A[i][j])); 
                a.Assemble();
                // store the block matrices of (Id - dt A)
                if(i == j){
                    MmdtA_SpMat[i][j] = m.SpMat(); 
                    MmdtA_SpMat[i][j].Add(-dt, a.SpMat());
                }
                else{
                    MmdtA_SpMat[i][j] = m0.SpMat(); 
                    MmdtA_SpMat[i][j].Add(-dt, a.SpMat());
                }      
            // store the block matrices of A
            A_SpMat[i][j] = a.SpMat();  
            MmdtA_BO.SetBlock(i, j, &MmdtA_SpMat[i][j]); 
            A_BO.SetBlock(i, j, &A_SpMat[i][j]); 
            }
        }
    }

    // ****************************************************************
    // Physical Space Solver 

    BilinearForm q(&fespace); 
    ConstantCoefficient eps_coeff(-1.0 * eps);
    VectorFunctionCoefficient u_coeff(dim, u, new ConstantCoefficient(-1.0));
    q.AddDomainIntegrator(new DiffusionIntegrator(eps_coeff)); 
    q.AddDomainIntegrator(new ConvectionIntegrator(u_coeff)); 
    q.Assemble(); 

    SparseMatrix Tmat = m.SpMat();
    Tmat.Add(-dt, q.SpMat()); 

    // ****************************************************************
    // Output 

    // Create references to the single phis 
    std::vector<GridFunction> phis(vector_size);
    for (int i = 0; i < vector_size; i++ ){
        phis[i].MakeRef(&fespace, phi_block.GetBlock(i),0);
    }

    // Prepare paraview output and save initial conditions  
    ParaViewDataCollection *pd = NULL;
    pd = new ParaViewDataCollection(scenario + "_N=" + to_string(N) + "_nx=" + to_string(n_x) + "_dt=" + to_string(dt), mesh);
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

    // Define two ODE solvers 
    ODESolver *ode_solver_css = NULL;
    ODESolver *ode_solver_pss = NULL;
    ode_solver_css = new BackwardEulerSolver;
    ode_solver_pss = new BackwardEulerSolver;

    // Physical space solver  
    PSS pss(q, Tmat); 
    pss.SetTime(t); 
    ode_solver_pss -> Init(pss); 

    // Configuration space solver 
    CSS css(A_BO, MmdtA_BO, vector_size);
    css.SetTime(t);
    ode_solver_css->Init(css);

    // ****************************************************************
    // Setup of the simulation 
    
    for (int ti = 0; !done; ){
      cout << "t: " << t << "s / " << t_final << "s - dt: " << dt << endl;  
      
      cout << "CSS "; 
      // apply css to the phi vector 
      ode_solver_css->Step(phi, t, dt);
      
      // advance iteration counter and save output 
      ti++;
      pd->SetCycle(ti);
      pd->SetTime(t);
      pd->Save();

      cout << "PSS " << endl; 
      tmp = t; 
      // apply pss to componentwise 
      for(int i = 0; i < vector_size; i++){
          ode_solver_pss -> Step(phi_block.GetBlock(i),t,dt); 
      }   
      t = tmp + dt; 

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