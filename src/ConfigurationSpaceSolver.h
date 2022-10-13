/*
    This time dependent operator deals with the configuration space
    d/dt psi = - div(grad(u) R psi) + div(xi R psi) + chi laplace(psi)

    The finite difference scheme d/dt phi = A phi 
    has a row full of zeros for phi_00 and thus leaves it constant 
    This requires some special treatment throughout everything. 
*/

class CSS : public TimeDependentOperator {

private:
    Operator &A_BO, &Id_p_delta_A_BO, &Id_m_delta_dt_A_BO; 

    int vector_size, N, n_dof; 

    BiCGSTABSolver css_solver;
    
    mutable Vector z, tmp;

public:

    CSS(BlockOperator &A_, BlockOperator &M_, BlockOperator &K_, int vector_size_): 
        TimeDependentOperator(M_.Height()), 
        A_BO(A_),
        Id_p_delta_A_BO(M_),
        Id_m_delta_dt_A_BO(K_),
        vector_size(vector_size_),
        N(sqrt(vector_size)), 
        n_dof(M_.Height() / (vector_size)), 
        z(M_.Height()), 
        tmp(n_dof){

        css_solver.iterative_mode = false;
        css_solver.SetRelTol(1e-12);
        css_solver.SetMaxIter(1000);
        css_solver.SetPrintLevel(0);
        css_solver.SetOperator(Id_m_delta_dt_A_BO); 
    }

    void Mult(const Vector &u, Vector &k) const{
        A_BO.Mult(u,k); 
    }

    void ImplicitSolve(const double dt, const Vector &phi_old, Vector &dphi_dt){
        // A * phi^n 
        Id_p_delta_A_BO.Mult(phi_old,z);

        // calculate phi^n+1
        css_solver.Mult(z,dphi_dt); 
    }
    
    ~CSS(){}
};
