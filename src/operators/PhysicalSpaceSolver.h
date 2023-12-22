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
        z(fespace_.GetTrueVSize()), 
        tmp(fespace_.GetTrueVSize()),
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
        
        calculate_operators(); 
    }

    void apply_Fx(const Vector &a, Vector &b) const{
        Fx_HPM->Mult(a,b);
    }

    void solve_Id_minus_beta_Fx(const Vector &y, Vector &x){

        m_HPM->Mult(y,z);   
        pss_solver.Mult(z, x); 

    }
    
    void set_current_block(int i){
        current_block = i; 
    }

    void calculate_operators(){

        Fx = new ParBilinearForm(&fespace); 
        Fx->AddDomainIntegrator(new DiffusionIntegrator(*eps_coeff)); 
        Fx->AddDomainIntegrator(new ConvectionIntegrator(u_coeff,-1.)); 
        Fx->AddBoundaryIntegrator(new MassIntegrator(*u_n_coeff)); 

        Fx->Assemble();
        Fx->Finalize(); 
        Fx_HPM = Fx->ParallelAssemble(); 
        
        m_1o_beta_coeff = new ConstantCoefficient(-1.0 / beta); 
        M_m_beta_Fx = new ParBilinearForm(&fespace); 

        M_m_beta_Fx->AddDomainIntegrator(new DiffusionIntegrator(*eps_coeff)); 
        M_m_beta_Fx->AddDomainIntegrator(new ConvectionIntegrator(u_coeff,-1.)); 
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