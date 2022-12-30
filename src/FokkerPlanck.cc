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

    Mesh *mesh = new Mesh();
    *mesh = Mesh::MakeCartesian2D(4, 4, Element::QUADRILATERAL);

    #if defined(Experiment1)
        *mesh = Mesh::MakeCartesian2D(8, 8, Element::Type::QUADRILATERAL);
    #endif

    #if defined(Experiment2)
        *mesh = Mesh::MakeCartesian2D(64, 64, Element::Type::QUADRILATERAL);
    #endif

    #if defined(Experiment3)
        *mesh = Mesh::MakeCartesian2D(4, 4, Element::Type::QUADRILATERAL);
    #endif

    #if defined(Experiment4)
        *mesh = Mesh::MakeCartesian2D(2, 2, Element::Type::QUADRILATERAL);
    #endif

    #if defined(Experiment5_pres_u)
        Mesh square = Mesh::MakeCartesian2D(10, 10, Element::QUADRILATERAL);
        square.Save("Exp5_non_periodic.mesh");
        Vector left_right_translation({1.0, 0.0});
        std::vector<Vector> translations = {left_right_translation};
        *mesh = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
        mesh->RemoveInternalBoundaries();
        mesh->Save("Exp5_periodic.mesh");
    #endif

    #if defined(Experiment5_pres_C)
        Mesh square = Mesh::MakeCartesian2D(10, 10, Element::QUADRILATERAL);
        square.Save("Exp5_non_periodic.mesh");
        Vector left_right_translation({1.0, 0.0});
        std::vector<Vector> translations = {left_right_translation};
        *mesh = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
        mesh->RemoveInternalBoundaries();
        mesh->Save("Exp5_periodic.mesh");
    #endif

    #if defined(Experiment6)
        mesh = new Mesh(mesh_file);
    #endif

    auto *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    for (int i = 0; i < n_refine; i++){
        pmesh->UniformRefinement();
    }
    delete mesh;

    // Define the finite element and the finite element space
    H1_FECollection fec(1, dim);
    FiniteElementSpace fespace(pmesh, &fec); // scalar
    FiniteElementSpace vfespace(pmesh, &fec, vector_size, Ordering::byNODES); // phi vector
    FiniteElementSpace v2dfespace(pmesh, &fec, 2, Ordering::byNODES); // 2D vector
    const int n_dof = fespace.GetVSize();

    NavierSolver flowsolver(pmesh, 2, nu);

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

    // Calculating T
    GridFunction *phi00 = &phis[0];
    GridFunction *phi02 = &phis[2];
    GridFunction *phi11 = &phis[1+N];
    GridFunction *phi20 = &phis[2*N];

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

    #if !defined(Experiment5_pres_C)
    VectorCoefficient *T = new VectorSumCoefficient(T_12, T_34);
    #endif 

    #if defined(Experiment5_pres_C)
        VectorFunctionCoefficient *T = new VectorFunctionCoefficient(dim, div_T); 
    #endif

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
        ProductCoefficient xi_coeff(trace_C_coeff, trace_C_coeff);
        ProductCoefficient chi_coeff(trace_C_coeff,trace_C_coeff);
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

    // ****************************************************************
    // Stuff for the output

    GridFunction *chi_gf = new GridFunction(&fespace);
    chi_gf->ProjectCoefficient(chi_coeff);

    GridFunction *xi_gf = new GridFunction(&fespace);
    xi_gf->ProjectCoefficient(xi_coeff);

    GridFunction *T_gf = new GridFunction(&v2dfespace);
    T_gf->ProjectCoefficient(*T);

    GridFunction *C11_gf = new GridFunction(&fespace);
    SumCoefficient *C11_coeff = new SumCoefficient(phi00_coeff, phi20_coeff, 25.132741228718345, 35.54306350526693);
    C11_gf->ProjectCoefficient(*C11_coeff);

    GridFunction *C12_gf = new GridFunction(&fespace);
    ProductCoefficient *C12_coeff = new ProductCoefficient(25.132741228718345, phi11_coeff);
    C12_gf->ProjectCoefficient(*C12_coeff);

    GridFunction *C22_gf = new GridFunction(&fespace);
    SumCoefficient *C22_coeff = new SumCoefficient(phi00_coeff, phi02_coeff, 25.132741228718345, 35.54306350526693);
    C22_gf->ProjectCoefficient(*C22_coeff);

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
    ParGridFunction *u_gf_NS = flowsolver.GetCurrentVelocity();
    ParGridFunction *p_gf_NS = flowsolver.GetCurrentPressure();

    // Calculating the divergence of u
    DivergenceGridFunctionCoefficient div_u_coeff(u_gf_NS);
    GridFunction *div_u_gf = new GridFunction(&fespace);
    div_u_gf->ProjectCoefficient(div_u_coeff);

    VectorGridFunctionCoefficient u_coeff(u_gf_NS);

    // ****************************************************************
    // Use a prescribed velocity 

    VectorFunctionCoefficient u_prescribed_coeff(dim, u_IC);
    if(prescribed_velocity){
        u_gf_NS->ProjectCoefficient(u_prescribed_coeff);
    }

    // ****************************************************************

    // to visualize the gradient of u
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
    

    // Configuration space solver
    CSS css(fespace, phi0_block, phi_modes, vector_size, block_offsets, u_gf_NS, chi_coeff, xi_coeff);
    css.SetTime(t);
    ode_solver_css->Init(css);

    // Physical space solver
    PSS pss(fespace, phi0_block, phi_modes, u_coeff);
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
        + "_dt=" + to_string(dt), pmesh);
    pd->SetPrefixPath("ParaView");

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
    // Time loop

    cout << "delta" << endl; 
    cout << delta << endl; 

    for (int ti = 0; !done; ){
        cout << "t: " << t << "s / " << t_final << "s - dt: " << dt << endl;

        // configuration space solver
        tmp = t;
        ode_solver_css->Step(phi_block, t, dt);

        if(!prescribed_velocity){
            t = tmp;
        }

        // calculate FR
        css.Mult(phi_block, F_R_block);

        // calculate M^-1 FR
        for(int i= 0; i < vector_size; i++){
            m_solver.Mult(F_R_block.GetBlock(i), tmp_vector);
            F_R_block.GetBlock(i) = tmp_vector;
        }

        // update modes for the configuration space
        for (int l = 0; l < n_modes; l++){
            phi_modes[l].Add(dt * weights[l], F_R_block);
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

        if(!prescribed_velocity){
            flowsolver.Step(t, dt, ti);
        }

        // recalculate the operators for the new u
        css.calculate_operators();
        pss.calculate_operators();

        // just for output
        div_u_gf->ProjectCoefficient(div_u_coeff);
        T_gf->ProjectCoefficient(*T);
        chi_gf->ProjectCoefficient(chi_coeff);
        xi_gf->ProjectCoefficient(xi_coeff);
        C11_gf->ProjectCoefficient(*C11_coeff);
        C12_gf->ProjectCoefficient(*C12_coeff);
        C22_gf->ProjectCoefficient(*C22_coeff);

        // advance iteration counter and save output
        ti++;
        if(ti % plot_frequency == 0){
            pd->SetCycle(ti);
            pd->SetTime(t);
            pd->Save();
        }

        done = (t >= t_final - 1e-8 * dt);
    }

    cout << "t: " << t << "s / " << t_final << "s" << endl;

    // Free the used memory.
    delete ode_solver_css;
    delete pd;
    return 0;
}