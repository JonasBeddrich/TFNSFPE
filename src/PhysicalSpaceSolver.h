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
    ProductCoefficient* beta_eps_coeff; 
    VectorGridFunctionCoefficient &u_coeff; 

    ParFiniteElementSpace &fespace; 
    Array<int> ess_tdof_list;

    ParBilinearForm* M; 
    ParBilinearForm* Fx; 
    ParBilinearForm* M_m_beta_Fx;

    HypreParMatrix M_mat;  
    HypreParMatrix Fx_mat;  
    HypreParMatrix M_m_beta_Fx_mat;  
    
    BlockVector &phi0; 
    std::vector<BlockVector> &phi_modes; 
    
    GMRESSolver pss_solver;
    
    mutable Vector z, tmp; 

public: 
    PSS(ParFiniteElementSpace &fespace_, BlockVector &phi0_, std::vector<BlockVector> &phi_modes_, VectorGridFunctionCoefficient &u_coeff_, double *t_):  
        fespace(fespace_),                    
        phi0(phi0_),
        phi_modes(phi_modes_), 
        z(fespace_.TrueVSize()), 
        tmp(fespace_.TrueVSize()),
        u_coeff(u_coeff_), 
        t(t_){
        
        Array<int> ess_bdr(fespace.GetParMesh()->bdr_attributes.Max());
        ess_bdr = 1;
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

        // weights for rational approximation 
        beta = get_beta(alpha, dt); 
        gammas = get_gammas(alpha, dt); 
        
        // mass matrix 
        M = new ParBilinearForm(&fespace);
        M->AddDomainIntegrator(new MassIntegrator); 
        M->Assemble(); 
        M->FormSystemMatrix(ess_tdof_list, M_mat);  

        // constant diffusion coefficient 
        eps_coeff = new ConstantCoefficient(-1.0 * eps); 
        
        calculate_operators(); 
    }

    void apply_Fx(const Vector &a, Vector &b) const{
        Fx_mat.Mult(a,b);
    }

    // This method solves (Id - beta Fx) x = y 
    void solve_Id_minus_beta_Fx(const Vector &y, Vector &x){
        // calculates x = (Id - beta Fx)^-1 M*y 

        // z = M*y 
        M_mat.Mult(y,z); 
        
        // x = (Id - beta Fx)^-1 z  
        pss_solver.Mult(z, x); 
    }
    
    void set_current_block(int i){
        current_block = i; 
    }

    void calculate_operators(){
        beta_eps_coeff = new ProductCoefficient(-beta, *eps_coeff); 

        // rhs        
        Fx = new ParBilinearForm(&fespace); 
        Fx->AddDomainIntegrator(new DiffusionIntegrator(*eps_coeff)); 
        Fx->AddDomainIntegrator(new ConvectionIntegrator(u_coeff)); 
        Fx->Assemble();
        Fx->FormSystemMatrix(ess_tdof_list, Fx_mat);  

        // lhs 
        M_m_beta_Fx = new ParBilinearForm(&fespace); 
        M_m_beta_Fx->AddDomainIntegrator(new MassIntegrator); 
        M_m_beta_Fx->AddDomainIntegrator(new ConvectionIntegrator(u_coeff, -beta)); 
        M_m_beta_Fx->AddDomainIntegrator(new DiffusionIntegrator(*beta_eps_coeff)); 
        M_m_beta_Fx->Assemble(); 
        M_m_beta_Fx->FormSystemMatrix(ess_tdof_list, M_m_beta_Fx_mat); 

        // solver 
        pss_solver.iterative_mode = false;
        pss_solver.SetRelTol(1e-12);
        pss_solver.SetMaxIter(1000);
        pss_solver.SetPrintLevel(0);
        pss_solver.SetOperator(M_m_beta_Fx_mat); 
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
