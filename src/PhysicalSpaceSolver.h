/*
    This time dependent operator deals with the physical space part 
    d/dt psi = u * grad(psi) + eps laplace(psi)
*/

// Physical Space Solver 
class PSS : public TimeDependentOperator {

private: 
    Operator &M, &K; 
    GMRESSolver T_solver;

    mutable Vector z; 

public: 
    PSS(Operator &M_, Operator &K_): 
        TimeDependentOperator(M_.Height()), 
        M(M_),          // mass matrix 
        K(K_),          // rhs matrix
        z(M_.Height()){
        
        T_solver.iterative_mode = false;
        T_solver.SetRelTol(1e-12);
        T_solver.SetMaxIter(1000);
        T_solver.SetPrintLevel(0);
        T_solver.SetOperator(K); 
    }

    void Mult(const Vector &u, Vector &k) const{}

    void ImplicitSolve(const double dt, const Vector &u, Vector &k){
        M.Mult(u,z);
        T_solver.Mult(z,k); 
    }
    
    ~PSS(){}
}; 
