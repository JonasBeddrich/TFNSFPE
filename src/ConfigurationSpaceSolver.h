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

    BilinearForm* a; 
    BilinearForm* m;
    BilinearForm* m0;
     
    std::vector<CoefficientFactory> *CoefficientVector;
    std::vector<std::vector<FunctionCoefficient>> *A; 
    std::vector<std::vector<bool>> *A_entries;

    std::vector<std::vector<SparseMatrix>> *A_SpMat;
    BlockOperator *A_BO; 
    
    std::vector<std::vector<SparseMatrix>> *L_SpMat;
    BlockOperator *L_BO; 

    BiCGSTABSolver css_solver;
    
    double delta; 
    int vector_size, N, n_dof; 

    mutable Vector z;


public:

    CSS(FiniteElementSpace &fespace_, int vector_size_, const Array<int> &offsets_): 
        TimeDependentOperator(fespace_.GetVSize() * vector_size_),
        z(fespace_.GetVSize() * vector_size_), 
        fespace(fespace_),  
        vector_size(vector_size_),
        n_dof(fespace_.GetVSize()), 
        N(sqrt(vector_size_)){

        delta = get_delta(dt); 

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
        L_SpMat = new std::vector<std::vector<SparseMatrix>>(vector_size, std::vector<SparseMatrix>(vector_size,m0->SpMat()));

        A_BO = new BlockOperator(offsets_); 
        L_BO = new BlockOperator(offsets_);

        for(int i = 0; i < vector_size; i++){
            for(int j = 0; j < vector_size; j++){
                if((*A_entries)[i][j]){
                    cout << i << " " << j << endl; 
                    BilinearForm a(&fespace);  
                    a.AddDomainIntegrator(new MassIntegrator((*A)[i][j])); 
                    a.Assemble();

                    (*A_SpMat)[i][j] = a.SpMat(); 

                    (*L_SpMat)[i][j] = a.SpMat();
                    (*L_SpMat)[i][j] *= - delta * dt; 
                    if(i == j){
                        (*L_SpMat)[i][j] += m->SpMat();
                    }

                    A_BO->SetBlock(i, j, &((*A_SpMat)[i][j])); 
                    L_BO->SetBlock(i, j, &((*L_SpMat)[i][j])); 
                }
            }
        } 

        for(int i = 0; i < vector_size; i++){
            for(int j = 0; j < vector_size; j++){
                if((*A_entries)[i][j]){
                    if(i == j){
                    } else {
                        cout << endl; 
                        cout << "A - MATRIX" << endl; 
                        cout << endl; 
                        cout << (*A_SpMat)[i][j] << endl; 
                        cout << endl; 
                        cout << "L - MATRIX" << endl; 
                        cout << endl; 
                        cout << (*L_SpMat)[i][j] << endl; 
                    }
                }
            }
        } 
        
        cout << "Delta times dt" << endl; 
        cout << - delta * dt << endl; 

        // cout << delta << " " << dt << endl; 

        // SparseMatrix* block = (SparseMatrix *) &(L_BO->GetBlock(13,11)); 
        // cout << *block << endl; 

        css_solver.iterative_mode = false;
        css_solver.SetRelTol(1e-12);
        css_solver.SetMaxIter(1000);
        css_solver.SetPrintLevel(0);
        css_solver.SetOperator(*L_BO);
    }

    void Mult(const Vector &phi, Vector &A_phi) const{
        A_BO->Mult(phi,A_phi); 
    }

    void ImplicitSolve(const double dt, const Vector &phi_old, Vector &dphi_dt){
        // A * phi^n 
        A_BO->Mult(phi_old,z);
        // calculate phi^n+1
        css_solver.Mult(z,dphi_dt); 
    }

    ~CSS(){}
};

