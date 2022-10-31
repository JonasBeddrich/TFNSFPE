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

    std::vector<std::vector<Coefficient*>> *coeff_matrix; 
    GridFunction &u1_gf; 
    GridFunction &u2_gf; 

    std::vector<std::vector<SparseMatrix>> *A_SpMat;
    BlockOperator *A_BO; 
    
    std::vector<std::vector<SparseMatrix>> *L_SpMat;
    BlockOperator *L_BO; 

    BiCGSTABSolver css_solver;
    
    double delta; 
    int vector_size, N, n_dof; 

    mutable Vector z, e1, e2;

    VectorConstantCoefficient *e1_coeff, *e2_coeff; 
    GradientGridFunctionCoefficient *u1_gf_coeff, *u2_gf_coeff;
    InnerProductCoefficient *d1u1_coeff, *d1u2_coeff, *d2u1_coeff, *d2u2_coeff; 
    SumCoefficient *A11_coeff, *A12_coeff, *A21_coeff, *A22_coeff; 

public:

    CSS(FiniteElementSpace &fespace_, int vector_size_, const Array<int> &offsets_, GridFunction &u1_gf, GridFunction &u2_gf): 
        TimeDependentOperator(fespace_.GetVSize() * vector_size_),
        z(fespace_.GetVSize() * vector_size_), 
        fespace(fespace_),  
        vector_size(vector_size_),
        n_dof(fespace_.GetVSize()), 
        N(sqrt(vector_size_)), 
        u1_gf(u1_gf), 
        u2_gf(u2_gf), 
        e1(2), 
        e2(2){

        setup_coefficients(); 

        // CoefficientVector = new std::vector<CoefficientFactory>(vector_size, CoefficientFactory()); 
        // fill_CoefficientVector(*CoefficientVector); 

        // A = new std::vector<std::vector<FunctionCoefficient>>(vector_size, std::vector<FunctionCoefficient>(vector_size, FunctionCoefficient(zero))); 
        A_entries = new std::vector<std::vector<bool>>(vector_size, std::vector<bool>(vector_size)); 
        // fill_A(*A, *A_entries, *CoefficientVector); 

        ConstantCoefficient *ptr_0_coeff = new ConstantCoefficient(0.); 
        coeff_matrix = new std::vector<std::vector<Coefficient*>>(vector_size, std::vector<Coefficient*>(vector_size, ptr_0_coeff)); 
        fill_coefficient_matrix(*coeff_matrix, *A_entries, A11_coeff, A12_coeff, A21_coeff, A22_coeff); 

        A_SpMat = new std::vector<std::vector<SparseMatrix>>(vector_size, std::vector<SparseMatrix>(vector_size,m0->SpMat()));
        L_SpMat = new std::vector<std::vector<SparseMatrix>>(vector_size, std::vector<SparseMatrix>(vector_size,m0->SpMat()));

        A_BO = new BlockOperator(offsets_); 
        L_BO = new BlockOperator(offsets_);

        for(int i = 0; i < vector_size; i++){
            for(int j = 0; j < vector_size; j++){
                
                if((*A_entries)[i][j]){
                    BilinearForm a(&fespace);  
                    a.AddDomainIntegrator(new MassIntegrator(*((*coeff_matrix)[i][j]))); 
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

    void setup_coefficients(){
        delta = get_delta(dt); 
        e1[0] = 1.;  
        e2[1] = 1.; 
        
        m = new BilinearForm(&fespace);
        m->AddDomainIntegrator(new MassIntegrator); 
        m->Assemble(); 

        m0 = new BilinearForm(&fespace);
        ConstantCoefficient m0_coef(0.);
        m0->AddDomainIntegrator(new MassIntegrator(m0_coef)); 
        m0->Assemble(); 

        e1_coeff = new VectorConstantCoefficient(e1); 
        e2_coeff = new VectorConstantCoefficient(e2); 

        u1_gf_coeff = new GradientGridFunctionCoefficient(&u1_gf);
        u2_gf_coeff = new GradientGridFunctionCoefficient(&u2_gf);

        d1u1_coeff = new InnerProductCoefficient(*u1_gf_coeff, *e1_coeff); 
        d2u1_coeff = new InnerProductCoefficient(*u1_gf_coeff, *e2_coeff); 
        d1u2_coeff = new InnerProductCoefficient(*u2_gf_coeff, *e1_coeff); 
        d2u2_coeff = new InnerProductCoefficient(*u2_gf_coeff, *e2_coeff); 

        A11_coeff = new SumCoefficient(xi, *d1u1_coeff, 1.0, -1.0); 
        A12_coeff = new SumCoefficient(0., *d2u1_coeff, 0.0, -1.0); 
        A21_coeff = new SumCoefficient(0., *d1u2_coeff, 0.0, -1.0); 
        A22_coeff = new SumCoefficient(xi, *d2u2_coeff, 1.0, -1.0); 
    }

    ~CSS(){}
};

