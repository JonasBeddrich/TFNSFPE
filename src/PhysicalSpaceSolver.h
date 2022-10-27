/*
    This time dependent operator deals with the physical space part 
    d/dt psi = u * grad(psi) + eps laplace(psi)
*/

// Physical Space Solver 
class PSS : public TimeDependentOperator {

private: 
    int current_block; 
    double beta; 
    std::vector<double> gammas; 

    ConstantCoefficient* eps_coeff; 
    VectorFunctionCoefficient* u_coeff; 

    FiniteElementSpace &fespace; 
    BilinearForm* m; 
    BilinearForm* Fx; 
    SparseMatrix* M_m_beta_Fx; 
    
    BlockVector &phi0; 
    std::vector<BlockVector> &phi_modes; 
    
    GMRESSolver pss_solver;
    
    mutable Vector z, tmp; 

public: 
    PSS(FiniteElementSpace &fespace_, BlockVector &phi0_, std::vector<BlockVector> &phi_modes_): 
        TimeDependentOperator(fespace_.GetVSize()), 
        fespace(fespace_),                    
        phi0(phi0_),
        phi_modes(phi_modes_), 
        z(fespace_.GetVSize()), 
        tmp(fespace_.GetVSize()){
        
        // weights for rational approximation 
        beta = get_beta(dt); 
        gammas = get_gammas(dt); 
        
        // mass matrix 
        m = new BilinearForm(&fespace);
        m->AddDomainIntegrator(new MassIntegrator); 
        m->Assemble(); 
 
        // coefficients for the physical space operator 
        eps_coeff = new ConstantCoefficient(-1.0 * eps); 
        u_coeff = new VectorFunctionCoefficient(dim, u, new ConstantCoefficient(-1.0)); 

        // physical space operator        
        Fx = new BilinearForm(&fespace); 
        Fx->AddDomainIntegrator(new DiffusionIntegrator(*eps_coeff)); 
        Fx->AddDomainIntegrator(new ConvectionIntegrator(*u_coeff)); 
        Fx->Assemble(); 
        
        // left side operator 
        M_m_beta_Fx = new SparseMatrix(m->SpMat());
        M_m_beta_Fx->Add(- beta, Fx->SpMat()); 

        // set up solver 
        pss_solver.iterative_mode = false;
        pss_solver.SetRelTol(1e-12);
        pss_solver.SetMaxIter(1000);
        pss_solver.SetPrintLevel(0);
        pss_solver.SetOperator(*M_m_beta_Fx); 
    }

    void Mult(const Vector &a, Vector &b) const{
        Fx->Mult(a,b);
    }

    void ImplicitSolve(const double dt, const Vector &phi_old, Vector &phi_up){
    
        // accumulate phi_0 and modes in z
        m->Mult(phi0.GetBlock(current_block), z);
        for (int l = 0; l < n_modes; l++){    
            m->Mult(phi_modes[l].GetBlock(current_block), tmp);    
            z.Add(gammas[l], tmp); 
        }

        // solver for phi^n+1
        pss_solver.Mult(z, tmp); 

        // calculate update 
        phi_up = tmp; 
        phi_up.Add(-1.0, phi_old); 
        phi_up /= dt; 

        // phi_up = 0.; 
    }
    
    void set_current_block(int i){
        current_block = i; 
    }

    ~PSS(){}
}; 
