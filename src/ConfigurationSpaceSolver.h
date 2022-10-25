/*
    This time dependent operator deals with the configuration space
    d/dt psi = - div(grad(u) R psi) + div(xi R psi) + chi laplace(psi)

    The finite difference scheme d/dt phi = A phi 
    has a row full of zeros for phi_00 and thus leaves it constant 
    This requires some special treatment throughout everything. 
*/

class CSS : public TimeDependentOperator {

private:
    FiniteElementSpace &fespace; 
    
    BlockOperator *my_A; 
    BlockOperator *my_Id_m; 

    BilinearForm* m;
    BilinearForm* m0;
     
    std::vector<CoefficientFactory> *CoefficientVector;
    std::vector<std::vector<FunctionCoefficient>> *A; 
    std::vector<std::vector<bool>> *A_entries;  
    
    std::vector<std::vector<SparseMatrix>> *Id_m_delta_dt_A_SpMat; 
    std::vector<std::vector<SparseMatrix>> *A_SpMat; 
    
    double delta; 

    Operator &A_BO, &Id_m_delta_dt_A_BO; 

    int vector_size, N, n_dof; 

    BiCGSTABSolver css_solver;
    
    mutable Vector z, tmp;

public:

    CSS(FiniteElementSpace &fespace_, BlockOperator &A_, BlockOperator &K_, int vector_size_): 
        TimeDependentOperator(A_.Height()),
        fespace(fespace_),  
        A_BO(A_),
        Id_m_delta_dt_A_BO(K_),
        vector_size(vector_size_),
        N(sqrt(vector_size)), 
        n_dof(A_.Height() / (vector_size)), 
        z(A_.Height()), 
        tmp(n_dof){

        m = new BilinearForm(&fespace);
        m->AddDomainIntegrator(new MassIntegrator); 
        m->Assemble(); 

        m0 = new BilinearForm(&fespace);
        ConstantCoefficient m0_coef(0.);
        m0->AddDomainIntegrator(new MassIntegrator(m0_coef)); 
        m0->Assemble(); 

        CoefficientVector = new std::vector<CoefficientFactory>(vector_size, CoefficientFactory()); 
        fill_CoefficientVector(*CoefficientVector); 

        A = new std::vector<std::vector<FunctionCoefficient>>(vector_size, std::vector<FunctionCoefficient>(vector_size, FunctionCoefficient(zero))); 
        A_entries = new std::vector<std::vector<bool>>(vector_size, std::vector<bool>(vector_size)); 
        fill_A(*A, *A_entries, *CoefficientVector); 

        A_SpMat = new std::vector<std::vector<SparseMatrix>>(vector_size, std::vector<SparseMatrix>(vector_size,m0->SpMat())); 


        // A_SpMat(vector_size, std::vector<SparseMatrix*>(vector_size)); 
        // Id_m_delta_dt_A_SpMat = new std::vector<std::vector<SparseMatrix>>(vector_size, std::vector<SparseMatrix>(vector_size)); 

        // Fill the block operators with A and (Id - dt A) 
        for(int i = 0; i < vector_size; i++){
            for(int j = 0; j < vector_size; j++){
                if(A_entries->at(i)[j]){
                    BilinearForm a(&fespace);  
                    a.AddDomainIntegrator(new MassIntegrator(A->at(i)[j])); 
                    a.Assemble();

                    // store the block matrices of (Id - dt A)
                    if(i == j){
                        // Id_m_delta_dt_A_SpMat->at(i)[j] = m->SpMat(); 
                        // Id_m_delta_dt_A_SpMat->at(i)[j].Add(- delta * dt, a.SpMat());
                    } else{
                        // Id_m_delta_dt_A_SpMat->at(i)[j] = m0->SpMat(); 
                        // Id_m_delta_dt_A_SpMat->at(i)[j].Add(- delta * dt, a.SpMat());
                    }
                    // cout << i << " " << j << endl; 
                    A_SpMat->at(i)[j] = a.SpMat(); 
                    // store the block matrices of A
                    // my_Id_m->SetBlock(i, j, &(Id_m_delta_dt_A_SpMat->at(i)[j])); 
                    // cout << "me" << endl;  
                    my_A->SetBlock(i, j, A_SpMat->at(i)[j]); 
                }
            }
        } 

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
        A_BO.Mult(phi_old,z);

        // calculate phi^n+1
        css_solver.Mult(z,dphi_dt); 
    }
    
    ~CSS(){}
};

