/*
    This time dependent operator deals with the configuration space
    d/dt psi = - div(grad(u) R psi) + div(xi R psi) + chi laplace(psi)

    The finite difference scheme d/dt phi = A phi 
    has a row full of zeros for phi_00 and thus leaves it constant 
    This requires some special treatment throughout everything. 
*/

class CSS : public TimeDependentOperator {

private:
    Operator &M, &K; 

    int vector_size, N, n_dof; 

    SparseMatrix T; 
    BiCGSTABSolver T_solver;
    
    mutable Vector z, tmp;

public:

    CSS(BlockOperator &M_, BlockOperator &K_, int vector_size_): 
        TimeDependentOperator(M_.Height()), 
        M(M_),              // A 
        K(K_),              // Id - dt A, where A is the matrix arising from the finite difference scheme  
        vector_size(vector_size_),
        N(sqrt(vector_size)), 
        n_dof(M_.Height() / (vector_size)), 
        z(M_.Height()), 
        tmp(n_dof){

        T_solver.iterative_mode = false;
        T_solver.SetRelTol(1e-12);
        T_solver.SetMaxIter(1000);
        T_solver.SetPrintLevel(0);
        T_solver.SetOperator(K); 
    }

    void Mult(const Vector &u, Vector &k) const{}

    void ImplicitSolve(const double dt, const Vector &phi_old, Vector &dphi_dt){
        // A * phi^n 
        M.Mult(phi_old,z);

        // calculate phi^n+1
        T_solver.Mult(z,dphi_dt); 
    }
    
    ~CSS(){}
};
