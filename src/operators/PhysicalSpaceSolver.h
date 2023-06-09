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

    ConstantCoefficient* eps_coeff; 
    ConstantCoefficient* m_1o_beta_coeff; 
    VectorGridFunctionCoefficient &u_coeff; 

    ParFiniteElementSpace &fespace; 
    ParBilinearForm* m; 
    ParBilinearForm* Fx;
    ParBilinearForm* M_m_beta_Fx;
    
    HypreParMatrix* m_HPM;  
    HypreParMatrix* Fx_HPM;  
    HypreParMatrix* M_m_beta_Fx_HPM;  
    
    BlockVector &phi0; 
    std::vector<BlockVector> &phi_modes; 
    
    GMRESSolver pss_solver;
    
    mutable Vector z, tmp; 

public: 
    PSS(ParFiniteElementSpace &fespace_, BlockVector &phi0_, std::vector<BlockVector> &phi_modes_, VectorGridFunctionCoefficient &u_coeff_):  
        fespace(fespace_),                    
        phi0(phi0_),
        phi_modes(phi_modes_), 
        z(fespace_.GetVSize()), 
        tmp(fespace_.GetVSize()),
        u_coeff(u_coeff_), 
        pss_solver(MPI_COMM_WORLD){
        
        // weights for rational approximation 
        beta = get_beta(alpha, dt); 
        gammas = get_gammas(alpha, dt); 
        
        // mass matrix 
        m = new ParBilinearForm(&fespace);
        m->AddDomainIntegrator(new MassIntegrator); 
        m->Assemble(); 
        m->Finalize(); 
        m_HPM = m->ParallelAssemble(); 
 
        // coefficients for the physical space operator 
        // eps_coeff = new ConstantCoefficient(eps); 
        eps_coeff = new ConstantCoefficient(-1.0 * eps); 
        m_1o_beta_coeff = new ConstantCoefficient(-1.0 / beta); 
        
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
        // physical space operator        
        Fx = new ParBilinearForm(&fespace); 
        Fx->AddDomainIntegrator(new DiffusionIntegrator(*eps_coeff)); 
        Fx->AddDomainIntegrator(new ConvectionIntegrator(u_coeff)); 
        Fx->Assemble();
        Fx->Finalize(); 
        Fx_HPM = Fx->ParallelAssemble(); 

        M_m_beta_Fx = new ParBilinearForm(&fespace); 
        M_m_beta_Fx->AddDomainIntegrator(new DiffusionIntegrator(*eps_coeff)); 
        M_m_beta_Fx->AddDomainIntegrator(new ConvectionIntegrator(u_coeff)); 
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

    ~PSS(){}
}; 