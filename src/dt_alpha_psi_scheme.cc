#include "mfem.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>


#include <setting.h>
#include <resources/experiment2.h>
#include <resources/experiment3.h>
#include <resources/experiment4.h>
#include <resources/experiment7.h>
#include <resources/experiment8.h>
#include <resources/experiments.h>
#include <resources/initialConditions.h>
#include <resources/initialConditionsAnalyticalSolution.h>
#include <resources/rational_approximation.h>
#include <operators/finite_difference_scheme.h>
#include <operators/PhysicalSpaceSolver.h>
#include <operators/ConfigurationSpaceSolver.h>

#include <../../mfem-4.5/miniapps/navier/navier_solver.hpp>
#include <../../mfem-4.5/miniapps/navier/navier_solver.cpp>

#include <chrono>

using namespace std;
using namespace std::chrono;
using namespace mfem;
using namespace navier;

int main(int argc, char *argv[]){
    auto start = high_resolution_clock::now();
    Mpi::Init(argc, argv);
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    bool verbose = myid == 0; 
    Hypre::Init();

    Mesh *mesh = new Mesh();
    load_mesh(mesh); 

    auto *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    for (int i = 0; i < n_refine; i++){
        pmesh->UniformRefinement();
    }
    delete mesh;

    ParaViewDataCollection *pd = NULL;
    pd = new ParaViewDataCollection("dt_a_PSI_" 
        + scenario
        + symmetry 
        + "_alpha=" + std::to_string(alpha)
        + "_N=" + to_string(N)
        + "_m=" + to_string(n_modes)
        + "_dt=" + to_string(dt), pmesh);
    pd->SetPrefixPath("ParaView");
    pd->SetLevelsOfDetail(1);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);

    // Define the finite element and the finite element space
    const int vector_size = N*N; 
    H1_FECollection fec(1, dim);
    ParFiniteElementSpace fespace(pmesh, &fec); // scalar
    ParFiniteElementSpace v2dfespace(pmesh, &fec, 2, Ordering::byNODES); // 2D vector
    ParFiniteElementSpace vfespace(pmesh, &fec, vector_size, Ordering::byNODES); // phi vector
    
    const int n_dof = fespace.GetVSize();
    const int n_true_dof = fespace.GetTrueVSize();
    
    // Create offset vector for phi vector
    // i.e. indices where the next phi_ij starts
    Array<int> block_offsets(vector_size+1);
    block_offsets[0]=0;
    for (int i = 1; i < vector_size + 1; i++){
        block_offsets[i] =  n_dof;
    }
    block_offsets.PartialSum();
    
    Array<int> true_block_offsets(vector_size+1);
    true_block_offsets[0]=0;
    for (int i = 1; i < vector_size + 1; i++){
        true_block_offsets[i] =  n_true_dof;
    }
    true_block_offsets.PartialSum();
    
    // Blockvector for eta, psi, psi0, eta0 is 0 anyway 
    BlockVector phi_eta_block(block_offsets);
    BlockVector phi_psi_block(block_offsets);
    BlockVector phi_psi0_block(block_offsets);

    BlockVector phi_Fpsi0_block(block_offsets);
    phi_eta_block = 0.;
    phi_psi_block = 0.;
    phi_psi0_block = 0; 

    BlockVector phi_eta_true_block(true_block_offsets); 
    BlockVector phi_psi_true_block(true_block_offsets); 
    BlockVector phi_psi0_true_block(true_block_offsets); 
    phi_eta_true_block = 0; 
    phi_psi_true_block = 0; 
    phi_psi0_true_block = 0; 

    // Blockvectors for F(phi) and so on ... better than calling everything tmp 
    BlockVector phi_Fpsi0_true_block(true_block_offsets);
    BlockVector phi_FR_true_block(true_block_offsets);
    BlockVector phi_Fx_true_block(true_block_offsets);
    phi_Fpsi0_true_block = 0;
    phi_FR_true_block = 0.;
    phi_Fx_true_block = 0.;

    // Needed to project the initial condition 
    ParGridFunction phi_psi0(&vfespace, phi_psi0_block.GetData());
    ParGridFunction phi_psi(&vfespace, phi_psi_block.GetData());
    NavierSolver flowsolver(pmesh, 2, nu);
    
    // Velocity initial conditions
    ParGridFunction *u_ic = flowsolver.GetCurrentVelocity();
    VectorFunctionCoefficient u_IC_coeff(pmesh->Dimension(), u_IC);
    u_ic->ProjectCoefficient(u_IC_coeff);
    
    // Load initial conditions and fill psiM
    VectorFunctionCoefficient phi_psiM_coeff(vector_size, phi_psiM_function);
    // phi_psiM.ProjectCoefficient(phi_psiM_coeff);
    phi_psi0.ProjectCoefficient(phi_psiM_coeff);
    phi_psi.ProjectCoefficient(phi_psiM_coeff);

    // Create references to the single phis
    std::vector<ParGridFunction> phis_psi0(vector_size); 
    std::vector<ParGridFunction> phis_eta(vector_size);
    std::vector<ParGridFunction> phis_psi(vector_size);
    std::vector<ParGridFunction> phis_Fpsi0(vector_size);

    for (int i = 0; i < vector_size; i++ ){
        phis_psi0[i].MakeRef(&fespace, phi_psi0_block.GetBlock(i),0);
        phis_eta[i].MakeRef(&fespace, phi_eta_block.GetBlock(i),0);
        phis_psi[i].MakeRef(&fespace, phi_psi_block.GetBlock(i),0);
        phis_Fpsi0[i].MakeRef(&fespace, phi_Fpsi0_block.GetBlock(i),0);
    }

    for (int i = 0; i < vector_size; i++ ){
        pd->RegisterField("phi_eta " + std::to_string(i), &phis_eta[i]);
        pd->RegisterField("phi_psi " + std::to_string(i), &phis_psi[i]);
        // pd->RegisterField("phi_Fpsi_0 " + std::to_string(i), &phis_Fpsi0[i]);f
    }

    for (int i = 0; i < vector_size; i++ ){
        phis_psi0[i].GetTrueDofs(phi_psi0_true_block.GetBlock(i)); 
    }

    // Create vector of blockvectors for modes
    std::vector<BlockVector> phi_eta_true_modes(n_modes,BlockVector(true_block_offsets));
    std::vector<BlockVector> phi_psi_true_modes(n_modes,BlockVector(true_block_offsets));
    for (int i = 0; i < n_modes; i++){
        phi_eta_true_modes[i] = 0.;
        phi_psi_true_modes[i] = 0.;
    }

 #if !defined(prescribed_velocity)
    // velocity boundary conditions
    Array<int> attr(pmesh->bdr_attributes.Max());
    attr = 1;
    for(int i = 0; i < pmesh->bdr_attributes.Max(); i++){
        attr[i] = boundary_D(i); 
    }
    flowsolver.AddVelDirichletBC(u_BC, attr);
#endif

    // Calculating T
    ParGridFunction *phi00 = &phis_psi[0];
    ParGridFunction *phi02 = &phis_psi[2];
    ParGridFunction *phi11 = &phis_psi[1+N];
    ParGridFunction *phi20 = &phis_psi[2*N];

    GradientGridFunctionCoefficient grad_phi00_coeff(phi00);
    GradientGridFunctionCoefficient grad_phi02_coeff(phi02);
    GradientGridFunctionCoefficient grad_phi11_coeff(phi11);
    GradientGridFunctionCoefficient grad_phi20_coeff(phi20);

    DenseMatrix swap(dim);
    swap = 0.;
    swap.Elem(0,1) = 25.132741228718345;
    swap.Elem(1,0) = 25.132741228718345;
    MatrixConstantCoefficient swap_coeff(swap);

    DenseMatrix elim_c1(dim);
    elim_c1 = 0.;
    elim_c1.Elem(1,1) = 35.54306350526693;
    MatrixConstantCoefficient elim_c1_coeff(elim_c1);

    DenseMatrix elim_c2(dim);
    elim_c2 = 0.;
    elim_c2.Elem(0,0) = 35.54306350526693;
    MatrixConstantCoefficient elim_c2_coeff(elim_c2);

    ScalarVectorProductCoefficient T_1(25.132741228718345, grad_phi00_coeff);
    MatrixVectorProductCoefficient T_2(swap_coeff, grad_phi11_coeff);
    MatrixVectorProductCoefficient T_3(elim_c1_coeff, grad_phi02_coeff);
    MatrixVectorProductCoefficient T_4(elim_c2_coeff, grad_phi20_coeff);

    VectorSumCoefficient T_12(T_1, T_2);
    VectorSumCoefficient T_34(T_3, T_4);

    // #if !defined(Experiment5_pres_C)
    double scale_T = 100; 
    VectorCoefficient *T = new VectorSumCoefficient(T_12,T_34,scale_T,scale_T);
    // #endif 

    // // #if defined(Experiment5_pres_C)
    // VectorFunctionCoefficient *T = new VectorFunctionCoefficient(dim, div_T); 
    // // #endif

    Array<int> domain_attr(pmesh->attributes.Max());
    domain_attr = 1;
    flowsolver.AddAccelTerm(T, domain_attr);

    // Calculating chi and xi coefficients
    GridFunctionCoefficient phi00_coeff(phi00);
    GridFunctionCoefficient phi02_coeff(phi02);
    GridFunctionCoefficient phi11_coeff(phi11);
    GridFunctionCoefficient phi20_coeff(phi20);

    SumCoefficient phi02n20_coeff(phi02_coeff, phi20_coeff);
    SumCoefficient trace_C_coeff(phi02n20_coeff, phi00_coeff, 35.54306350526693, 2 * 25.132741228718345);

    // ****************************************************************
    // Set parameters xi and chi 

    #if defined(Experiment1)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
    #endif

    #if defined(Experiment2)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
    #endif

    #if defined(Experiment3)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
    #endif

    #if defined(Experiment4)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
        // ProductCoefficient xi_coeff(trace_C_coeff, trace_C_coeff);
        // ProductCoefficient chi_coeff(trace_C_coeff,trace_C_coeff);
    #endif

     #if defined(Experiment4_Medea)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
    #endif

    #if defined(Experiment5_pres_u)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
    #endif

    #if defined(Experiment5_pres_C)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
    #endif

    #if defined(Experiment6)
        ProductCoefficient xi_coeff(trace_C_coeff, trace_C_coeff);
        ProductCoefficient chi_coeff(1.0, trace_C_coeff);
    #endif

    #if defined(Experiment7)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
        // ProductCoefficient xi_coeff(trace_C_coeff, trace_C_coeff);
        // ProductCoefficient chi_coeff(1.0, trace_C_coeff);
    #endif

    #if defined(Experiment8)
        ConstantCoefficient one_coeff(1.0);
        ProductCoefficient xi_coeff(1.0, one_coeff);
        ProductCoefficient chi_coeff(1.0, one_coeff);
        // ProductCoefficient xi_coeff(trace_C_coeff, trace_C_coeff);
        // ProductCoefficient chi_coeff(1.0, trace_C_coeff);
    #endif

    // ****************************************************************
    // Stuff for the output

    ParGridFunction *chi_gf = new ParGridFunction(&fespace);
    chi_gf->ProjectCoefficient(chi_coeff);

    ParGridFunction *xi_gf = new ParGridFunction(&fespace);
    xi_gf->ProjectCoefficient(xi_coeff);

    ParGridFunction *T_gf = new ParGridFunction(&v2dfespace);
    T_gf->ProjectCoefficient(*T);

    ParGridFunction *C11_gf = new ParGridFunction(&fespace);
    SumCoefficient *C11_coeff = new SumCoefficient(phi00_coeff, phi20_coeff, 25.132741228718345, 35.54306350526693);
    C11_gf->ProjectCoefficient(*C11_coeff);

    ParGridFunction *C12_gf = new ParGridFunction(&fespace);
    ProductCoefficient *C12_coeff = new ProductCoefficient(25.132741228718345, phi11_coeff);
    C12_gf->ProjectCoefficient(*C12_coeff);

    ParGridFunction *C22_gf = new ParGridFunction(&fespace);
    SumCoefficient *C22_coeff = new SumCoefficient(phi00_coeff, phi02_coeff, 25.132741228718345, 35.54306350526693);
    C22_gf->ProjectCoefficient(*C22_coeff);

    // ****************************************************************
    // Rational Approximation for the TFNSFP system  

    std::vector<double> lambdas = get_lambdas(alpha);
    std::vector<double> weights = get_weights(alpha);
    double w_inf = get_w_infinity(alpha);

    std::vector<double> gammas = get_gammas(alpha, dt);
    double beta = get_beta(alpha, dt);

    for(auto i: gammas){
        cout << i << " "; 
    }
    cout << endl; 

    for(auto i: weights){
        cout << i << " "; 
    }
    cout << endl; 

    for(auto i: lambdas){
        cout << i << " "; 
    }
    cout << endl; 

    cout << w_inf << endl; 
   
    // ****************************************************************
    // Setup of the simulation

    Vector tmp_vector(n_true_dof);
    BlockVector tmp_block_vector(true_block_offsets);
    BlockVector tmp2_block_vector(true_block_offsets);

    // Get velocity and pressure from Navier-Stokes solver 
    ParGridFunction *u_gf_NS = flowsolver.GetCurrentVelocity();
    ParGridFunction *p_gf_NS = flowsolver.GetCurrentPressure();
    VectorGridFunctionCoefficient u_coeff(u_gf_NS);

    // Divergernce of u 
    DivergenceGridFunctionCoefficient div_u_coeff(u_gf_NS);
    ParGridFunction *div_u_gf = new ParGridFunction(&fespace);
    div_u_gf->ProjectCoefficient(div_u_coeff);

    // ****************************************************************
    // Use a prescribed velocity 
    VectorFunctionCoefficient u_prescribed_coeff(dim, u_IC);
    #if defined(prescribed_velocity)
        u_gf_NS->ProjectCoefficient(u_prescribed_coeff);
    #endif 

    // ****************************************************************
    // Visualization of the gradient of u 
    const FiniteElementCollection *u_fec = u_gf_NS->ParFESpace()->FEColl(); 
    ParFiniteElementSpace u_fes(pmesh, u_fec, 1); 

    ParGridFunction *dxu1_gf = new ParGridFunction(&u_fes); 
    ParGridFunction *dyu1_gf = new ParGridFunction(&u_fes); 
    ParGridFunction *dxu2_gf = new ParGridFunction(&u_fes); 
    ParGridFunction *dyu2_gf = new ParGridFunction(&u_fes); 
    
    u_gf_NS->GetDerivative(1,0,*dxu1_gf); 
    u_gf_NS->GetDerivative(1,1,*dyu1_gf); 
    u_gf_NS->GetDerivative(2,0,*dxu2_gf); 
    u_gf_NS->GetDerivative(2,1,*dyu2_gf); 
    
    GridFunctionCoefficient d1u1_coeff(dxu1_gf); 
    GridFunctionCoefficient d2u1_coeff(dyu1_gf); 
    GridFunctionCoefficient d1u2_coeff(dxu2_gf); 
    GridFunctionCoefficient d2u2_coeff(dyu2_gf); 

    // Navier Stokes
    flowsolver.Setup(dt_IC_vel);

    ParBilinearForm m(&fespace);
    m.AddDomainIntegrator(new MassIntegrator);
    m.Assemble();
    m.Finalize(); 
    HypreParMatrix *m_HPM = m.ParallelAssemble(); 

    BiCGSTABSolver m_solver(MPI_COMM_WORLD);
    m_solver.iterative_mode=false;
    m_solver.SetRelTol(1e-12);
    m_solver.SetMaxIter(1000);
    m_solver.SetPrintLevel(0);
    m_solver.SetOperator(*m_HPM);

    // ****************************************************************
    // Output

    // Prepare paraview output and save initial conditions
    pd->RegisterField("_velocity", u_gf_NS);
    pd->RegisterField("_pressure", p_gf_NS);

    pd->RegisterField("_div_u", div_u_gf);
    pd->RegisterField("dxu1", dxu1_gf);
    pd->RegisterField("dyu1", dyu1_gf);
    pd->RegisterField("dxu2", dxu2_gf);
    pd->RegisterField("dyu2", dyu2_gf);
    
    pd->RegisterField("_chi", chi_gf);
    pd->RegisterField("_T", T_gf);

    pd->RegisterField("_C11", C11_gf);
    pd->RegisterField("_C12", C12_gf);
    pd->RegisterField("_C22", C22_gf);
    
    // ****************************************************************
    // Calculate initial conditions - velocity field  

#if defined(calculate_initial_velocity_field)

    if(verbose){
        cout << "" << endl; 
        cout << ">>> COMPUTING INITIAL VELOCITY FIELD <<<" << endl; 
    }

    for (int i = 0; i < n_iter_IC_vel; i++){
        if(verbose){
            cout << endl; 
            cout << "Initial velocity field - Iteration: " << i << endl; 
        }
        double j=-10; 
        flowsolver.Step(j, dt_IC_vel, 0);
    }
    
    if(verbose){
        cout << endl; 
        cout << "Recalculate operators for initial velocity field" << endl; 
    }

#endif 
    
    // ****************************************************************
    // Calculate initial conditions - probability density  

    int ti_total = 0; 
    phi_psi_true_block = phi_psi0_true_block; 

#if defined(calculate_initial_condition)

    if(verbose){
        cout << "" << endl; 
        cout << ">>> COMPUTING INITIAL CONDITION FOR PSI <<<" << endl; 
        cout << "" << endl; 
    }
    
    // Operators for the integer-order problem (Id - dt * F)
    CSS css_IC(fespace, vector_size, true_block_offsets, u_gf_NS, chi_coeff, xi_coeff, dt_IC);
    PSS pss_IC(fespace, u_coeff, dt_IC);

    // css_IC.set_theta(dt_IC); 
    // pss_IC.set_beta(dt_IC); 
    
    int ti_IC = 0; 
    double t_IC = 0.0; 
    bool done_IC = false; 

    for (int ti_IC = 0; !done_IC; ){
        
        for(int i = 0; i < vector_size; i++){
            phis_eta[i].Distribute(phi_eta_true_block.GetBlock(i)); 
            phis_psi[i].Distribute(phi_psi_true_block.GetBlock(i)); 
            phis_Fpsi0[i].Distribute(phi_Fpsi0_true_block.GetBlock(i)); 
        }

        if(ti_IC % plot_frequency_IC == 0){
            pd->SetCycle(ti_total);
            pd->SetTime(-(t_final_IC - t_IC));
            pd->Save();
        }

        if(verbose){
            cout << "t: " << t_IC << "s / " << t_final_IC << "s - dt: " << dt_IC << endl;
        }        

        css_IC.solve_Id_minus_theta_FR(phi_psi0_true_block, tmp_block_vector);                

        for(int i = 0; i < vector_size; i++){
            pss_IC.solve_Id_minus_beta_Fx(tmp_block_vector.GetBlock(i), phi_psi0_true_block.GetBlock(i)); 
        }

        phi_psi_true_block = phi_psi0_true_block; 

        t_IC += dt_IC; 
        ti_IC++;
        ti_total++;         
        
        done_IC = (t_final_IC <= t_IC + 1e-8 * dt_IC);
    }

    // if(verbose){
    //     cout << endl; 
    //     cout << "Store initial condition for psi, "; 
    // }

    // std::ofstream outstream;
    // outstream.open("Initial_Condition"+scenario + symmetry + "_N=" + to_string(N));
    // outstream << phi_psi0_block.Size() << endl; 
    // phi_psi0_block.Print(outstream);
    // outstream.close();  

#endif 

    // cout << "Load Initial Condition" << endl; 
    // std::ifstream inputstream;     
    // inputstream.open("Initial_Condition"+scenario + symmetry + "_N=" + to_string(N)); 
    // phi_psi0_block.Load(inputstream); 
    // inputstream.close(); 
    // phi_psi_block = phi_psi0_block; 

    // ****************************************************************
    // Time loop

    double t = 0.0;
    bool done = false; 

    if(verbose){
        cout << ">>> OUTPUT INITIAL CONDITIONS <<<" << endl; 
        cout << "" << endl; 
    }

    phi_psi_true_block = phi_psi0_true_block; 

    pd->SetCycle(ti_total);
    pd->SetTime(t);
    pd->Save();

    if(verbose){
        cout << ">>> START TIME LOOP <<<" << endl; 
        cout << "" << endl; 
    }
    std::ofstream outstream;

    for (int ti = 0; !done; ){
        
        if(verbose){
            cout << "t: " << t << "s / " << t_final << "s - dt: " << dt << endl;
        }

        #if !defined(prescribed_velocity)
            flowsolver.Step(t, dt, ti);
            t -= dt; // set back to have t^n = t and t^n+1 = t + dt 
        #endif 
        
        // Calculate the derivatives of u for the output 
        u_gf_NS->GetDerivative(1,0,*dxu1_gf); 
        u_gf_NS->GetDerivative(1,1,*dyu1_gf); 
        u_gf_NS->GetDerivative(2,0,*dxu2_gf); 
        u_gf_NS->GetDerivative(2,1,*dyu2_gf); 

         // Project coefficients on grid functions for the output
        dxu1_gf->ProjectCoefficient(d1u1_coeff); 
        dyu1_gf->ProjectCoefficient(d2u1_coeff); 
        dxu2_gf->ProjectCoefficient(d1u2_coeff); 
        dyu2_gf->ProjectCoefficient(d2u2_coeff); 
        div_u_gf->ProjectCoefficient(div_u_coeff);
        
        // ORIGINAL VERSION 
        // CSS css(fespace, vector_size, true_block_offsets, u_gf_NS, chi_coeff, xi_coeff, get_theta(alpha, dt));
        // PSS pss(fespace, u_coeff, get_beta(alpha, dt));
        
        // CORRECTED VERSION 
        CSS css(fespace, vector_size, true_block_offsets, u_gf_NS, chi_coeff, xi_coeff, get_beta(alpha, dt));
        PSS pss(fespace, u_coeff, get_beta(alpha, dt)); 

        // cout << "theta: " << get_theta(alpha,dt); 
        // cout << " / beta: " << get_beta(alpha,dt) << endl; 
        
        // ****************************************************************
        // Solve the configuration space problem 
        
        // summation over modes and right-hand side 
        tmp_block_vector = phi_psi0_true_block; 
        
        for (int l = 0; l < n_modes; l++){
            // ORIGINAL VERSION 
            // tmp_block_vector.Add(1., phi_psi_true_modes[l]); 
            
            // CORRECTED VERSION 
            tmp_block_vector.Add(gammas[l], phi_psi_true_modes[l]); 
        }

        // solve for eta^n+0.5 
        css.solve_Id_minus_theta_FR(tmp_block_vector, phi_psi_true_block); 
        
        // calculate M^-1 FR eta 
        css.apply_FR(phi_psi_true_block, tmp_block_vector); 
        for(int i= 0; i < vector_size; i++){
            m_solver.Mult(tmp_block_vector.GetBlock(i), phi_FR_true_block.GetBlock(i));
        }

        // update modes 
        for (int l = 0; l < n_modes; l++){
            phi_psi_true_modes[l].Add(dt * weights[l], phi_FR_true_block);
        }

        // ****************************************************************
        // Solve the physical space problem  

        // summation over modes and right-hand side 
        tmp_block_vector = phi_psi0_true_block;

        for (int l = 0; l < n_modes; l++){
            tmp_block_vector.Add(gammas[l], phi_psi_true_modes[l]); 
        }

        // solve for eta^n+1 and store it 
        for(int i = 0; i < vector_size; i++){
            pss.solve_Id_minus_beta_Fx(tmp_block_vector.GetBlock(i), phi_psi_true_block.GetBlock(i)); 
        }

        // calculate M^⁻1 Fx eta   
        for(int i = 0; i < vector_size; i++){
            pss.apply_Fx(phi_psi_true_block.GetBlock(i), tmp_block_vector.GetBlock(i)); 
            m_solver.Mult(tmp_block_vector.GetBlock(i), phi_Fx_true_block.GetBlock(i));
        }
        
        // update modes 
        for (int l = 0; l < n_modes; l++){
            phi_psi_true_modes[l].Add(dt * weights[l], phi_Fx_true_block);
            phi_psi_true_modes[l] *= gammas[l];  
        }

        
        // ****************************************************************
        // Update time step and check if done 
        
        t += dt;         
        done = (t >= t_final - 1e-8 * dt);

        // ****************************************************************
        // Load output 

        for(int i = 0; i < vector_size; i++){
            phis_eta[i].Distribute(phi_eta_true_block.GetBlock(i)); 
            phis_psi[i].Distribute(phi_psi_true_block.GetBlock(i)); 
            phis_Fpsi0[i].Distribute(phi_Fpsi0_true_block.GetBlock(i)); 
        }

        // Calculate the derivatives of u for the output 
        u_gf_NS->GetDerivative(1,0,*dxu1_gf); 
        u_gf_NS->GetDerivative(1,1,*dyu1_gf); 
        u_gf_NS->GetDerivative(2,0,*dxu2_gf); 
        u_gf_NS->GetDerivative(2,1,*dyu2_gf); 

        // Project coefficients on grid functions for the output
        dxu1_gf->ProjectCoefficient(d1u1_coeff); 
        dyu1_gf->ProjectCoefficient(d2u1_coeff); 
        dxu2_gf->ProjectCoefficient(d1u2_coeff); 
        dyu2_gf->ProjectCoefficient(d2u2_coeff); 
        div_u_gf->ProjectCoefficient(div_u_coeff);
        T_gf->ProjectCoefficient(*T);
        chi_gf->ProjectCoefficient(chi_coeff);
        xi_gf->ProjectCoefficient(xi_coeff);
        C11_gf->ProjectCoefficient(*C11_coeff);
        C12_gf->ProjectCoefficient(*C12_coeff);
        C22_gf->ProjectCoefficient(*C22_coeff);

        div_u_gf->Distribute(div_u_gf->GetTrueDofs());         
        T_gf->Distribute(T_gf->GetTrueDofs());         

        // Advance iteration counters and write output
        ti++;
        ti_total++; 
        if(ti % plot_frequency == 0){
            pd->SetCycle(ti_total);
            pd->SetTime(t);
            pd->Save();
        }        
    }
    
    if(verbose){
        cout << "t: " << t << "s / " << t_final << "s" << endl;
        cout << "" << endl; 
        cout << ">>> SIMULATION FINISHED <<<" << endl; 
        cout << "" << endl; 
    }

    // Free the used memory.
    // delete pd;
    // delete m_HPM; 

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    
    if(verbose){
        cout << "The simulation took " << duration.count()/1000000. << " seconds."<< endl;
    }

    return 0;
}