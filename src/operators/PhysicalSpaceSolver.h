/*
    This time dependent operator deals with the physical space part 
    d/dt psi = u * grad(psi) + eps laplace(psi)
*/

// Physical Space Solver 
class PSS {

private: 
    int current_block; 
    double beta; 
    std::vector<double> gammas; 

    double h; 
    double (coth_ptr)(double); 

    MatrixCoefficient* uuT_coeff; 
    Coefficient* uTu_coeff; 

    ScalarVectorProductCoefficient* tau_u_coeff; 
    ScalarMatrixProductCoefficient* tau_uuT_coeff; 
    ScalarMatrixProductCoefficient* beta_tau_uuT_coeff; 

    VectorFunctionCoefficient* n_coeff; 
    InnerProductCoefficient* u_n_coeff; 

    ConstantCoefficient* negative_coeff; 
    ConstantCoefficient* h_coeff; 
    ConstantCoefficient* h_o_2_coeff; 
    ConstantCoefficient* h_by_2eps_coeff; 

    ProductCoefficient* Pe1; 
    ProductCoefficient* Pe2; 

    ProductCoefficient* beta_eps_coeff; 
    ConstantCoefficient* eps_coeff; 
    ConstantCoefficient* m_1o_beta_coeff; 
    VectorGridFunctionCoefficient &u_coeff; 

    ConstantCoefficient* two_coeff; 

    ParFiniteElementSpace &fespace; 
    ParBilinearForm* m; 
    ParBilinearForm* Fx;
    ParBilinearForm* M_m_beta_Fx;
    

    ParBilinearForm* RHS;
    HypreParMatrix* RHS_HPM; 

    ParBilinearForm* Stab_LHS; 
    HypreParMatrix* Stab_LHS_HPM; 

    HypreParMatrix* m_HPM;  
    HypreParMatrix* Fx_HPM;  
    HypreParMatrix* M_m_beta_Fx_HPM;  
    
    BiCGSTABSolver pss_solver;
    
    mutable Vector z, tmp; 

public: 
    PSS(ParFiniteElementSpace &fespace_, VectorGridFunctionCoefficient &u_coeff_, double beta_):  
        fespace(fespace_),                    
        z(fespace_.GetVSize()), 
        tmp(fespace_.GetVSize()),
        u_coeff(u_coeff_), 
        beta(beta_), 
        pss_solver(MPI_COMM_WORLD){
        
        // mass matrix 
        m = new ParBilinearForm(&fespace);
        m->AddDomainIntegrator(new MassIntegrator); 
        m->Assemble(); 
        m->Finalize(); 
        m_HPM = m->ParallelAssemble(); 
 
        eps_coeff = new ConstantCoefficient(-1.0 * eps); 

        negative_coeff = new ConstantCoefficient(1.0); 
        n_coeff = new VectorFunctionCoefficient(2, n_bdr,negative_coeff);
        u_n_coeff = new InnerProductCoefficient(*n_coeff, u_coeff); 
        
        // TODO 
        h = 0.25; 
        h_coeff = new ConstantCoefficient(h); 
        h_o_2_coeff = new ConstantCoefficient(h/2); 
        h_by_2eps_coeff = new ConstantCoefficient(h/(2*eps)); 

        uuT_coeff = new OuterProductCoefficient(u_coeff, u_coeff);
        uTu_coeff = new InnerProductCoefficient(u_coeff, u_coeff); 
        
        Vector e1(2); 
        Vector e2(2);
        e1[0] = 1;  
        e2[1] = 1; 

        VectorConstantCoefficient e1_coeff(e1); 
        VectorConstantCoefficient e2_coeff(e2);  
        
        InnerProductCoefficient u1_coeff(u_coeff, e1_coeff); 
        InnerProductCoefficient u2_coeff(u_coeff, e2_coeff);
        
        Pe1 = new ProductCoefficient(u1_coeff, *h_by_2eps_coeff); 
        Pe2 = new ProductCoefficient(u1_coeff, *h_by_2eps_coeff); 
        
        TransformedCoefficient coth_Pe1(Pe1, &coth); 
        TransformedCoefficient coth_Pe2(Pe2, &coth); 

        RatioCoefficient m1_o_Pe1(-1., *Pe1); 
        RatioCoefficient m1_o_Pe2(-1., *Pe2); 

        SumCoefficient f_Pe1_coeff(coth_Pe1, m1_o_Pe1); 
        SumCoefficient f_Pe2_coeff(coth_Pe2, m1_o_Pe2);

        ProductCoefficient f_Pe1_u1_coeff(f_Pe1_coeff, u1_coeff); 
        ProductCoefficient f_Pe2_u2_coeff(f_Pe2_coeff, u2_coeff); 

        ProductCoefficient h_o_2_f_Pe1_u1_coeff(f_Pe1_u1_coeff, *h_o_2_coeff); 
        ProductCoefficient h_o_2_f_Pe2_u2_coeff(f_Pe2_u2_coeff, *h_o_2_coeff); 

        SumCoefficient v_bar_coeff(h_o_2_f_Pe1_u1_coeff, h_o_2_f_Pe2_u2_coeff);
        RatioCoefficient tau_coeff(v_bar_coeff, *uTu_coeff); 

        tau_u_coeff = new ScalarVectorProductCoefficient(tau_coeff, u_coeff); 

        tau_uuT_coeff = new ScalarMatrixProductCoefficient(tau_coeff, *uuT_coeff); 
        beta_tau_uuT_coeff = new ScalarMatrixProductCoefficient(beta, *tau_uuT_coeff); 

        beta_eps_coeff = new ProductCoefficient(beta, *eps_coeff); 

        two_coeff = new ConstantCoefficient(2.); 

        calculate_operators(); 
    }

    void apply_Fx(const Vector &a, Vector &b) const{
        // cout << Fx->NumCols() << " " << Fx_HPM->NumCols() << " "; 
        // cout << a.Size() << " " << b.Size() << endl; 
        Fx_HPM->Mult(a,b);
    }

    // This method solves (Id - beta Fx) x = y 
    void solve_Id_minus_beta_Fx(const Vector &y, Vector &x){
        // calculates x = (Id - beta Fx)^-1 M*y 
        
        // z = M*y 
        m_HPM->Mult(y,z); 
        // RHS_HPM->Mult(y,z); 
        
        // x = (Id - beta Fx)^-1 z  
        pss_solver.Mult(z, x); 

        // for (int l = 0; l < n_modes; l++){    
        //     m->Mult(phi_modes[l].GetBlock(current_block), tmp);    
        //     z.Add(gammas[l], tmp); 
        // } 

        // // This function solves (Id-beta Fx)^-1  
        // // accumulate phi_0 and modes in z
        // m->Mult(phi0.GetBlock(current_block), z);
        // for (int l = 0; l < n_modes; l++){    
        //     m->Mult(phi_modes[l].GetBlock(current_block), tmp);    
        //     z.Add(gammas[l], tmp); 
        // }

        // // solver for phi^n+1
        // pss_solver.Mult(z, tmp); 

        // // calculate update 
        // phi_up = tmp; 
        // phi_up.Add(-1.0, phi_old); 
        // phi_up /= dt; 
    }
    
    void set_current_block(int i){
        current_block = i; 
    }

    void calculate_operators(){
        Fx = new ParBilinearForm(&fespace); 

        Fx->AddDomainIntegrator(new DiffusionIntegrator(*eps_coeff)); 
        Fx->AddDomainIntegrator(new ConvectionIntegrator(u_coeff,-1.)); 
        // boundary part ... let's see 
        Fx->AddBoundaryIntegrator(new MassIntegrator(*u_n_coeff)); 

        Fx->Assemble();
        Fx->Finalize(); 
        Fx_HPM = Fx->ParallelAssemble(); 
        
        m_1o_beta_coeff = new ConstantCoefficient(-1.0 / beta); 
        M_m_beta_Fx = new ParBilinearForm(&fespace); 

        M_m_beta_Fx->AddDomainIntegrator(new DiffusionIntegrator(*eps_coeff)); 
        M_m_beta_Fx->AddDomainIntegrator(new ConvectionIntegrator(u_coeff,-1.)); 
        // boundary part ... let's see 
        M_m_beta_Fx->AddBoundaryIntegrator(new MassIntegrator(*u_n_coeff)); 

        M_m_beta_Fx->AddDomainIntegrator(new MassIntegrator(*m_1o_beta_coeff)); 
        M_m_beta_Fx->Assemble();
        M_m_beta_Fx->Finalize(); 
        M_m_beta_Fx_HPM = M_m_beta_Fx->ParallelAssemble(); 
        (*M_m_beta_Fx_HPM) *= (-1.0 * beta); 

        // set up solver 
        pss_solver.iterative_mode = false;
        pss_solver.SetRelTol(1e-12);
        pss_solver.SetMaxIter(1000);
        pss_solver.SetPrintLevel(0);
        pss_solver.SetOperator(*M_m_beta_Fx_HPM); 

        // RHS = new ParBilinearForm(&fespace); 
        // RHS->AddDomainIntegrator(new ConservativeConvectionIntegrator(*tau_u_coeff, -1)); 
        // RHS->AddDomainIntegrator(new MassIntegrator()); 
        // RHS->Assemble();
        // RHS->Finalize(); 
        // RHS_HPM = RHS->ParallelAssemble(); 

        // Stab_LHS = new ParBilinearForm(&fespace); 
        // Stab_LHS->AddDomainIntegrator(new ConservativeConvectionIntegrator(*tau_u_coeff, -1. )); 
        // Stab_LHS->AddDomainIntegrator(new DiffusionIntegrator(*beta_tau_uuT_coeff)); 
        // Stab_LHS->AddDomainIntegrator(new MassIntegrator()); 
        // Stab_LHS->AddDomainIntegrator(new ConvectionIntegrator(u_coeff, beta)); 
        // Stab_LHS->AddDomainIntegrator(new DiffusionIntegrator()); 
        // Stab_LHS->Assemble();
        // Stab_LHS->Finalize(); 
        // Stab_LHS_HPM = Stab_LHS->ParallelAssemble(); 

        // // set up solver 
        // pss_solver.iterative_mode = false;
        // pss_solver.SetRelTol(1e-12);
        // pss_solver.SetMaxIter(1000);
        // pss_solver.SetPrintLevel(0);
        // pss_solver.SetOperator(*Stab_LHS_HPM); 
    }

    void set_beta(double beta_){
        beta = beta_; 
        calculate_operators(); 
    }

    void reset_beta(){
        beta = get_beta(alpha, dt); 
        calculate_operators(); 
    }

    ~PSS(){
        delete m; 
        delete m_HPM; 
        delete eps_coeff; 
        delete Fx; 
        delete M_m_beta_Fx; 
        delete Fx_HPM; 
        delete M_m_beta_Fx_HPM; 
        delete n_coeff; 
        delete u_n_coeff; 
    }
}; 