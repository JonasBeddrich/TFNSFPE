/*
    This time dependent operator deals with the physical space part 
    d/dt psi = u * grad(psi) + eps laplace(psi)
*/

// Physical Space Solver 
class PSS : public TimeDependentOperator {

private: 
    Operator &M, &Fx, &Id_m_beta_Fx; 
    GMRESSolver pss_solver;
    
    BlockVector &phi0; 
    std::vector<BlockVector> &phi_modes; 
    int current_block; 

    double beta; 
    std::vector<double> gammas; 

    mutable Vector z, tmp; 

public: 
    PSS(Operator &M_, Operator &Fx_, Operator &Id_m_beta_Fx_, BlockVector &phi0_, std::vector<BlockVector> &phi_modes_): 
        TimeDependentOperator(Fx_.Height()), 
        M(M_), 
        Fx(Fx_),                    
        Id_m_beta_Fx(Id_m_beta_Fx_),      
        phi0(phi0_),
        phi_modes(phi_modes_), 
        z(Fx_.Height()), 
        tmp(Fx_.Height()){
        
        beta = get_beta(dt); 
        gammas = get_gammas(dt); 

        pss_solver.iterative_mode = false;
        pss_solver.SetRelTol(1e-12);
        pss_solver.SetMaxIter(1000);
        pss_solver.SetPrintLevel(0);
        pss_solver.SetOperator(Id_m_beta_Fx); 
    }

    void Mult(const Vector &a, Vector &b) const{
        Fx.Mult(a,b);
    }

    void ImplicitSolve(const double dt, const Vector &phi_old, Vector &phi_up){
        // let's violate some principles ... 

        M.Mult(phi0.GetBlock(current_block), z);
        for (int l = 0; l < n_modes; l++){    
            M.Mult(phi_modes[l].GetBlock(current_block), tmp);    
            z.Add(gammas[l], tmp); 
        }

        pss_solver.Mult(z, tmp); 

        phi_up = (tmp - phi_old) / dt; 

        phi_up = 0.; 
    }
    
    void set_current_block(int i){
        current_block = i; 
    }

    ~PSS(){}
}; 
