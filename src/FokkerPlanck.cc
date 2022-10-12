#include "mfem.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional> 

#include <experiment_settings.h>
#include <finite_difference_scheme.h>

using namespace std;
using namespace mfem;

//***********************************************************************************

// Physical Space Solver 
class PSS : public TimeDependentOperator {
private: 
  Operator &M, &K; 
  GMRESSolver T_solver;

  mutable Vector z; 

public: 
  PSS(Operator &M, Operator &K); 
  virtual void Mult(const Vector &x, Vector &y) const;
  virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);
  virtual ~PSS();
}; 

// Configuration Space Solver 
class CSS : public TimeDependentOperator {
private:
  Operator &M, &K, &A02, &A11, &A20; 
  int vector_size, N, n_dof; 
  SparseMatrix T; 
  BiCGSTABSolver T_solver;
  BlockDiagonalPreconditioner T_prec;
  mutable Vector z, tmp;

public:
  CSS(BlockOperator &M_, BlockOperator &K_, int vector_size_, Operator &A02_, Operator &A11_, Operator &A20_, Array<int> &block_offsets_partial);
  virtual void Mult(const Vector &x, Vector &y) const;
  virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);
  virtual ~CSS();
};

//***********************************************************************************

int main(int argc, char *argv[]){
  // 1. Options 
  const int vector_size = N*N; 

  // 2. Create mesh 
  Mesh *mesh = new Mesh(); 
  *mesh = Mesh::MakeCartesian2D(n_x, n_x, Element::Type::TRIANGLE);
  int dim = mesh -> Dimension();

  // 3. Define the ODE solver 
  ODESolver *ode_solver_css = NULL;
  ODESolver *ode_solver_pss = NULL;
  ode_solver_css = new BackwardEulerSolver;
  ode_solver_pss = new BackwardEulerSolver;

  // 5. Define the FECollection and FiniteElementSpace 
  H1_FECollection fec(1, dim); 
  FiniteElementSpace fespace(mesh, &fec); 
  FiniteElementSpace vfespace(mesh, &fec, vector_size, Ordering::byNODES); 

  // create offset vector for all phis 
  Array<int> block_offsets(vector_size+1); 
  block_offsets[0]=0; 
  for (int i = 1; i < vector_size + 1; i++){
    block_offsets[i] =  fespace.GetVSize(); 
  }
  block_offsets.PartialSum(); 

  // create offset vector for all phis except phi_0 = phi_00 
  int vector_size_partial = vector_size - 1; 
  Array<int> block_offsets_partial(block_offsets); 
  block_offsets_partial.DeleteLast(); 

  // create blockvetor for all phis 
  BlockVector phi_block(block_offsets); 
  phi_block = 1.0; 

  // 5.9 initial conditions
  GridFunction phi(&vfespace, phi_block.GetData());  
  VectorFunctionCoefficient phi0(vector_size, phi0_function);  
  phi.ProjectCoefficient(phi0); 

  // 6. Set up the bilinear and linear forms 
  // 6.1 Configuration Space Solver  

  std::vector<std::vector<FunctionCoefficient>> A(vector_size_partial, std::vector<FunctionCoefficient>(vector_size_partial, FunctionCoefficient(zero))); 
  std::vector<std::vector<bool>> A_entries(vector_size_partial, std::vector<bool>(vector_size_partial)); 
  std::vector<CoefficientFactory> CoefficientVector(vector_size, CoefficientFactory()); 

  fill_CoefficientVector(CoefficientVector); 
  fill_A(A, A_entries, CoefficientVector); 
  
  BlockOperator MmdtA_BO(block_offsets_partial); 
  BlockOperator A_BO(block_offsets_partial); 

  std::vector<std::vector<SparseMatrix>> MmdtA_SpMat(vector_size_partial, std::vector<SparseMatrix>(vector_size_partial)); 
  std::vector<std::vector<SparseMatrix>> A_SpMat(vector_size_partial, std::vector<SparseMatrix>(vector_size_partial)); 

  BilinearForm m(&fespace);
  m.AddDomainIntegrator(new MassIntegrator); 
  m.Assemble(); 

  BilinearForm m0(&fespace); 
  ConstantCoefficient m0_coef(0.);
  m0.AddDomainIntegrator(new MassIntegrator(m0_coef)); 
  m0.Assemble(); 

  for(int i = 0; i < vector_size_partial; i++){
    for(int j = 0; j < vector_size_partial; j++){
      if(A_entries[i][j]){
        BilinearForm a(&fespace);  
        a.AddDomainIntegrator(new MassIntegrator(A[i][j])); 
        a.Assemble();

        if(i == j){
          MmdtA_SpMat[i][j] = m.SpMat(); 
          MmdtA_SpMat[i][j].Add(-dt, a.SpMat());
        }
        else{
          MmdtA_SpMat[i][j] = m0.SpMat(); 
          MmdtA_SpMat[i][j].Add(-dt, a.SpMat());
        }
      
        A_SpMat[i][j] = a.SpMat();  

        MmdtA_BO.SetBlock(i, j, &MmdtA_SpMat[i][j]); 
        A_BO.SetBlock(i, j, &A_SpMat[i][j]); 
      }
    }
  }

  // Special treatment for the phi_00 since the matrix contains a row of zeros 
  auto function02 = std::bind(&CoefficientFactory::Eval_z_km2, CoefficientVector[2], std::placeholders::_1, std::placeholders::_2); 
  auto function11 = std::bind(&CoefficientFactory::Eval_zm1_km1, CoefficientVector[1 + N], std::placeholders::_1, std::placeholders::_2); 
  auto function20 = std::bind(&CoefficientFactory::Eval_zm2_k, CoefficientVector[2 * N], std::placeholders::_1, std::placeholders::_2); 
  
  FunctionCoefficient a02(function02); 
  FunctionCoefficient a11(function11); 
  FunctionCoefficient a20(function20); 

  BilinearForm A02(&fespace); 
  BilinearForm A11(&fespace); 
  BilinearForm A20(&fespace); 

  A02.AddDomainIntegrator(new MassIntegrator(a02)); 
  A11.AddDomainIntegrator(new MassIntegrator(a11)); 
  A20.AddDomainIntegrator(new MassIntegrator(a20)); 

  A02.Assemble(); 
  A11.Assemble(); 
  A20.Assemble(); 

  // 6.2 Physical Space Solver   
  BilinearForm q(&fespace); 
  ConstantCoefficient eps_coeff(-1.0 * eps);
  VectorFunctionCoefficient u_coeff(dim, u, new ConstantCoefficient(-1.0));
  q.AddDomainIntegrator(new DiffusionIntegrator(eps_coeff)); 
  q.AddDomainIntegrator(new ConvectionIntegrator(u_coeff)); 
  q.Assemble(); 

  SparseMatrix Tmat = m.SpMat();
  Tmat.Add(-dt, q.SpMat()); 

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
  
  // Setup of the simulation 
  double t = 0.0;
  double tmp = 0.0; 
  bool done = false;

  // Time-dependent evolution operator describing the ODE 
  PSS pss(q, Tmat); 
  pss.SetTime(t); 
  ode_solver_pss -> Init(pss); 

  CSS css(A_BO, MmdtA_BO, vector_size, A02, A11, A20, block_offsets_partial);
  css.SetTime(t);
  ode_solver_css->Init(css);

  for (int ti = 0; !done; ){

    cout << "t: " << t << "s / " << t_final << "s - dt: " << dt << endl;  
  
    cout << "CSS "; 
    ode_solver_css->Step(phi, t, dt);
    
    ti++;
    pd->SetCycle(ti);
    pd->SetTime(t);
    pd->Save();

    cout << "PSS " << endl; 
    tmp = t; 
    for(int i = 0; i < vector_size; i++){
      ode_solver_pss -> Step(phi_block.GetBlock(i),t,dt); 
    }   
    t = tmp + dt; 

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

//***********************************************************************************

PSS::PSS(Operator &M_, Operator &K_): 
  TimeDependentOperator(M_.Height()), 
  M(M_), 
  K(K_), 
  z(M_.Height()){

  T_solver.iterative_mode = false;
  T_solver.SetRelTol(1e-12);
  T_solver.SetMaxIter(1000);
  T_solver.SetPrintLevel(0);
  T_solver.SetOperator(K); 
}

void PSS::Mult(const Vector &u, Vector &k) const{}

void PSS::ImplicitSolve(const double dt, const Vector &u, Vector &k){
  M.Mult(u,z);
  T_solver.Mult(z,k); 
}

PSS::~PSS(){}

//***********************************************************************************

CSS::CSS(BlockOperator &M_, BlockOperator &K_, int vector_size_, Operator &A02_, Operator &A11_, Operator &A20_, Array<int> &block_offsets_partial): 
  TimeDependentOperator(M_.Height()), 
  M(M_), 
  K(K_),
  A02(A02_), 
  A11(A11_), 
  A20(A20_), 
  vector_size(vector_size_),
  N(sqrt(vector_size)), 
  n_dof(M_.Height() / (vector_size - 1)), 
  z(M_.Height()), 
  tmp(n_dof),
  T_prec(block_offsets_partial){

  T_solver.iterative_mode = false;
  T_solver.SetRelTol(1e-12);
  T_solver.SetMaxIter(1000);
  T_solver.SetPrintLevel(0);
  // T_solver.SetPreconditioner(T_prec);
  T_solver.SetOperator(K); 
}

void CSS::Mult(const Vector &u, Vector &k) const{}

void CSS::ImplicitSolve(const double dt, const Vector &phi_old, Vector &dphi_dt){
  // Compute: dphi_dt = f(phi_old + dt dphi_dt, t)
  
  Vector phi0_old(phi_old.GetData() + 0, n_dof); 
  Vector phir_old(phi_old.GetData() + n_dof, phi_old.Size() - n_dof); 
  Vector phi0_up(dphi_dt.GetData() + 0, n_dof); 
  Vector phir_up(dphi_dt.GetData() + n_dof, dphi_dt.Size() - n_dof); 

  // z doesn't contain phi_0 thus every entry has to be taken - 1 * n_dof    
  Vector phi02(z.GetData() + (2     - 1)  * n_dof, n_dof); 
  Vector phi11(z.GetData() + (N + 1 - 1)  * n_dof, n_dof); 
  Vector phi20(z.GetData() + (2 * N - 1)  * n_dof, n_dof); 

  M.Mult(phir_old,z);

  A02.Mult(phi0_old, tmp); 
  phi02.Add(1, tmp); 

  A11.Mult(phi0_old, tmp);   
  phi11.Add(1, tmp); 

  A20.Mult(phi0_old, tmp); 
  phi20.Add(1, tmp);  

  T_solver.Mult(z,phir_up); 

  phi0_up = 0.; 
}

CSS::~CSS(){}