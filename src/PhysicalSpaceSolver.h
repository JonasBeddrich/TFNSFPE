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

    double *t; 

    ConstantCoefficient* eps_coeff; 
    VectorGridFunctionCoefficient &u_coeff; 

    ParFiniteElementSpace &fespace; 
    ParBilinearForm* m; 
    ParBilinearForm* Fx; 
    SparseMatrix* M_m_beta_Fx; 
    
    BlockVector &phi0; 
    std::vector<BlockVector> &phi_modes; 
    
    GMRESSolver pss_solver;
    
    mutable Vector z, tmp; 

public: 
    PSS(ParFiniteElementSpace &fespace_, BlockVector &phi0_, std::vector<BlockVector> &phi_modes_, VectorGridFunctionCoefficient &u_coeff_, double *t_):  
        fespace(fespace_),                    
        phi0(phi0_),
        phi_modes(phi_modes_), 
        z(fespace_.GetVSize()), 
        tmp(fespace_.GetVSize()),
        u_coeff(u_coeff_), 
        t(t_){
        
        // weights for rational approximation 
        beta = get_theta(dt); 
        gammas = get_gammas(dt); 
        
        // mass matrix 
        m = new ParBilinearForm(&fespace);
        m->AddDomainIntegrator(new MassIntegrator); 
        m->Assemble(); 
 
        // coefficients for the physical space operator 
        // eps_coeff = new ConstantCoefficient(eps); 
        eps_coeff = new ConstantCoefficient(-1.0 * eps); 
        
        calculate_operators(); 
    }

    void apply_Fx(const Vector &a, Vector &b) const{
        Fx->Mult(a,b);
    }

    // This method solves (Id - beta Fx) x = y 
    void solve_Id_minus_beta_Fx(const Vector &y, Vector &x){
        // calculates x = (Id - beta Fx)^-1 M*y 
        
        // z = M*y 
        m->Mult(y,z); 
        
        // x = (Id - beta Fx)^-1 z  
        pss_solver.Mult(z, x); 

        // for (int l = 0; l < n_modes; l++){    
        //     m->Mult(phi_modes[l].GetBlock(current_block), tmp);    
        //     z.Add(gammas[l], tmp); 
        // } 

        // // This function solves (Id-beta Fx)^-1  
        // // cout << "t: " << *t << endl; 
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

    ~PSS(){}
}; 
