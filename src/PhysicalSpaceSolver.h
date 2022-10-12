/*
    This time dependent operator deals with the physical space part 
    d/dt psi = u * grad(psi) + eps laplace(psi)
*/

// Physical Space Solver 
class PSS : public TimeDependentOperator {

private: 
    Operator &M, &Fx, &Id_m_dtFx; 
    GMRESSolver pss_solver;
    
    BlockVector phi0; 
    int current_block; 

    mutable Vector z, tmp; 

public: 
    PSS(Operator &M_, Operator &Fx_, Operator &Id_m_dtFx_, BlockVector phi0_): 
        TimeDependentOperator(Fx_.Height()), 
        M(M_), 
        Fx(Fx_),                    
        Id_m_dtFx(Id_m_dtFx_),      
        phi0(phi0_),
        z(Fx_.Height()), 
        tmp(Fx_.Height()){
        
        pss_solver.iterative_mode = false;
        pss_solver.SetRelTol(1e-12);
        pss_solver.SetMaxIter(1000);
        pss_solver.SetPrintLevel(0);
        pss_solver.SetOperator(Id_m_dtFx); 
    }

    void Mult(const Vector &a, Vector &b) const{}

    void ImplicitSolve(const double dt, const Vector &phi_old, Vector &phi_up){
        M.Mult(phi0.GetBlock(current_block), tmp);
        Fx.Mult(phi_old,z);
        
        // ready 
        // z += tmp; 
        
        pss_solver.Mult(z,phi_up); 
    }
    
    void set_current_block(int i){
        current_block = i; 
    }

    ~PSS(){}
}; 
