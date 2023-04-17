/*
    This time dependent operator deals with the configuration space
    d/dt psi = - div(grad(u) R psi) + div(xi R psi) + chi laplace(psi)

    The finite difference scheme d/dt phi = A phi 
    has a row full of zeros for phi_00 and thus leaves it constant 
    This requires some special treatment throughout everything. 
*/

class CSS {

private:

    double theta; 
    int vector_size, N, n_dof; 

    ParFiniteElementSpace &fespace; 
    const FiniteElementCollection *u_fec; 
    ParFiniteElementSpace *u_fes; 

    BilinearForm* m;
    BilinearForm* m0;

    ParGridFunction *u_gf_NS; 
    ParGridFunction *dxu1_gf, *dyu1_gf, *dxu2_gf, *dyu2_gf;
    
    BlockVector &phi0; 
    std::vector<BlockVector> &phi_modes; 
    
    std::vector<std::vector<bool>> *A_entries;
    std::vector<std::vector<Coefficient*>> *coeff_matrix; 
    std::vector<std::vector<SparseMatrix>> *A_SpMat;
    std::vector<std::vector<SparseMatrix>> *L_SpMat;
    
    BlockOperator *M_BO; 
    BlockOperator *A_BO; 
    BlockOperator *L_BO; 
    BlockOperator *A44;
    BlockOperator *A42;
     

    ProductCoefficient &chi_coeff, &xi_coeff; 
    GridFunctionCoefficient *d1u1_coeff, *d1u2_coeff, *d2u1_coeff, *d2u2_coeff; 
    SumCoefficient *A11_coeff, *A12_coeff, *A21_coeff, *A22_coeff, *du_sym_coeff, *skew_12_coeff, *skew_21_coeff; 
    ProductCoefficient *a_a_2_chi_coeff; 
    ConstantCoefficient *m0_coeff;

    BiCGSTABSolver css_solver;
    BiCGSTABSolver A44_solver;
    mutable Vector z, tmp, test;

public:

    CSS(ParFiniteElementSpace &fespace_, BlockVector &phi0_, std::vector<BlockVector> &phi_modes_, int vector_size_, const Array<int> &offsets_, ParGridFunction *u_gf_NS_, ProductCoefficient &chi_coeff_, ProductCoefficient &xi_coeff_): 
        z(fespace_.GetVSize() * vector_size_), 
        tmp(fespace_.GetVSize() * vector_size_), 
        test(fespace_.GetVSize() * vector_size_), 
        fespace(fespace_),  
        phi0(phi0_),
        phi_modes(phi_modes_), 
        vector_size(vector_size_),
        n_dof(fespace_.GetVSize()), 
        N(sqrt(vector_size_)), 
        u_gf_NS(u_gf_NS_), 
        chi_coeff(chi_coeff_), 
        xi_coeff(xi_coeff_){

        M_BO = new BlockOperator(offsets_); 
        A_BO = new BlockOperator(offsets_); 
        L_BO = new BlockOperator(offsets_);

        setup_coefficients(); 

        A_entries = new std::vector<std::vector<bool>>(vector_size, std::vector<bool>(vector_size)); 
        fill_A_entries(*A_entries); 

        ConstantCoefficient *ptr_0_coeff = new ConstantCoefficient(0.); 
        coeff_matrix = new std::vector<std::vector<Coefficient*>>(vector_size, std::vector<Coefficient*>(vector_size, ptr_0_coeff)); 
        
        A_SpMat = new std::vector<std::vector<SparseMatrix>>(vector_size, std::vector<SparseMatrix>(vector_size,m0->SpMat()));
        L_SpMat = new std::vector<std::vector<SparseMatrix>>(vector_size, std::vector<SparseMatrix>(vector_size,m0->SpMat()));        

        calculate_operators(); 
    }

    void apply_FR(const Vector &phi, Vector &A_phi) const{
        // setup_block_operators(); 
        A_BO->Mult(phi,A_phi); 
    }

    void solve_Id_minus_theta_FR(const BlockVector &y, BlockVector &x){

        // ****************************************************************
        // calculates x = (Id - theta FR)^-1 M*y 

        // z = M*y 
        M_BO->Mult(y, z); 

        // x = (Id - theta FR)^-1 z
        css_solver.Mult(z,x); 

        // ****************************************************************
        // The following code is concerned with solving A_4,4 phi = A_4,2 mu 

        for (int counter = 3; counter < N; counter = counter +2){
            // Setup 
            Array<int> part_offsets(counter+1);
            part_offsets[0] = 0; 
            for (int i = 1; i < counter + 1; i++){
                part_offsets[i] = n_dof;
            }
            part_offsets.PartialSum();

            BlockVector mu_k(part_offsets); 
            BlockVector phi_k(part_offsets); 
            BlockVector phi_km2(part_offsets);
            BlockVector A42_phi_km2(part_offsets); 

            A44 = new BlockOperator(part_offsets); 
            A42 = new BlockOperator(part_offsets); 
            
            // Load mu k 
            // 0,4 - 1,3 - 2,2 - 3,1 - 4,0 

            for(int i=0; i < counter; i++){
                int z = i; 
                int k = (counter-1) -i; 
                // cout << z << " " << k << " " << z*N + k << endl;  
                m->Mult(y.GetBlock(z*N + k), mu_k.GetBlock(i)); 
            }

            // Load phi km2 
            // 0,2 - 1,1 - 2,0 

            for(int i=0; i < counter-2; i++){
                int z = i; 
                int k = (counter-1) -2 -i; 
                // cout << z << " " << k << " " << z*N + k << endl;  
                m->Mult(x.GetBlock(z*N + k), phi_km2.GetBlock(i)); 
            }

            // Load A_4,2

            for(int i=0; i < counter; i++){ // loop rows 
                for(int j=0; j < counter-2; j++){ // loop columns, where the indices add up to z+k-2 
                    
                    int z_i = i; 
                    int k_i = (counter-1) -i; 
                    
                    int z_j = j;
                    int k_j = (counter-1) -2 -j;   

                    // cout << "i " << " " << z_i << " " << k_i << endl; 
                    // cout << "j " << " " << z_j << " " << k_j << endl; 

                    if((*A_entries)[z_i * N + k_i][z_j * N + k_j]){
                        // cout << "A entries is true: " << z_i << " " << k_i << " " << z_j << " " << k_j << endl; 
                        A42->SetBlock(i,j,&(L_BO->GetBlock(z_i * N + k_i,z_j * N + k_j))); 
                    } 
                }
            }

            // Compute A_4,2 phi 
            A42->Mult(phi_km2, A42_phi_km2); 

            // mu_k = mu_k + A_42 phi_km2
            mu_k += A42_phi_km2; 

            // Load A_4,4
            for(int i=0; i < counter; i++){ // loop rows 
                for(int j=0; j < counter; j++){ // loop columns 
                    
                    int z_i = i; 
                    int k_i = (counter-1) -i; 
                    
                    int z_j = j;  
                    int k_j = (counter-1) -j; 

                    if((*A_entries)[z_i * N + k_i][z_j * N + k_j]){
                        A44->SetBlock(i,j,&(L_BO->GetBlock(z_i * N + k_i,z_j * N + k_j))); 
                    } 
                }
            }

            // Solve A44 ... 
            A44_solver.iterative_mode = false;
            A44_solver.SetRelTol(1e-12);
            A44_solver.SetMaxIter(2000);
            A44_solver.SetPrintLevel(0);
            A44_solver.SetOperator(*A44); 

            A44_solver.Mult(mu_k, phi_k); 

            // Write phi km2 
            // 0,4 - 1,3 - 2,2 - 3,1 - 4,0 

            for(int i=0; i < counter; i++){
                int z = i; 
                int k = (counter-1) -i; 
                // cout << z << " " << k << endl; 

                for (int j = 0; j < n_dof; j++){
                    // cout << phi_k.GetBlock(i)[j] << " " << x.GetBlock(z*N + k)[j] << endl; 
                }
                x.GetBlock(z*N + k) = phi_k.GetBlock(i); 
                for (int j = 0; j < n_dof; j++){
                    // cout << phi_k.GetBlock(i)[j] << " " << x.GetBlock(z*N + k)[j] << endl; 
                }
                // cout << z << " " << k << " " << z*N + k << endl;  
                // phi_k.GetBlock(i) = y.GetBlock(z*N + k); 
            }
        }

        // for(int i=0; i < 3; i++){
        //     for(int j=0; j < 3; j++){
        //         if((*A_entries)[i][j]){
        //             partial->SetBlock(i,j,&(L_BO->GetBlock(i,j))); 
        //         }
        //     }
        // }

        // cout << "multiplication" << endl; 
        // partial->Mult(asdf, bsdf); 

        // phi0 = mu0 

        // for i = 2, 4, ... N 

        // tmp_offsets = ... 

        // Array<int> block_offsets(vector_size+1);
        // block_offsets[0]=0;
        // for (int i = 1; i < vector_size + 1; i++){
        //     block_offsets[i] =  n_dof;
        // }
        // block_offsets.PartialSum();

        // LBO = new BlockOperator( tmp_offsets )

    }

    void setup_coefficients(){

        theta = get_theta(alpha, dt); 
        
        m = new ParBilinearForm(&fespace);
        m->AddDomainIntegrator(new MassIntegrator); 
        m->Assemble(); 

        for(int i = 0; i < vector_size; i++){
            M_BO->SetBlock(i, i, m); 
        }

        m0 = new ParBilinearForm(&fespace);
        m0_coeff = new ConstantCoefficient(0.);
        m0->AddDomainIntegrator(new MassIntegrator(*m0_coeff)); 
        m0->Assemble(); 

        u_fec = u_gf_NS->ParFESpace()->FEColl(); 
        u_fes = new ParFiniteElementSpace(u_gf_NS->ParFESpace()->GetParMesh(), u_fec, 1); 

        dxu1_gf = new ParGridFunction(u_fes); 
        dyu1_gf = new ParGridFunction(u_fes); 
        dxu2_gf = new ParGridFunction(u_fes); 
        dyu2_gf = new ParGridFunction(u_fes); 
        
        u_gf_NS->GetDerivative(1,0,*dxu1_gf); 
        u_gf_NS->GetDerivative(1,1,*dyu1_gf); 
        u_gf_NS->GetDerivative(2,0,*dxu2_gf); 
        u_gf_NS->GetDerivative(2,1,*dyu2_gf); 

        d1u1_coeff = new GridFunctionCoefficient(dxu1_gf); 
        d2u1_coeff = new GridFunctionCoefficient(dyu1_gf); 
        d1u2_coeff = new GridFunctionCoefficient(dxu2_gf); 
        d2u2_coeff = new GridFunctionCoefficient(dyu2_gf); 

        A11_coeff = new SumCoefficient(xi_coeff, *d1u1_coeff, 1.0, -1.0); 
        A12_coeff = new SumCoefficient(0., *d2u1_coeff, 0.0, -1.0); 
        A21_coeff = new SumCoefficient(0., *d1u2_coeff, 0.0, -1.0); 
        A22_coeff = new SumCoefficient(xi_coeff, *d2u2_coeff, 1.0, -1.0); 
        
        #if defined(symmetric)
            du_sym_coeff = new SumCoefficient(*d2u1_coeff, *d1u2_coeff, 0.5, 0.5); 

            A12_coeff = new SumCoefficient(0., *du_sym_coeff, 0.0, -1.0); 
            A21_coeff = new SumCoefficient(0., *du_sym_coeff, 0.0, -1.0); 
        #endif

        #if defined(skew_symmetric)
            skew_12_coeff = new SumCoefficient(*d2u1_coeff, *d1u2_coeff, 0.5, -0.5); // d2u1-d1u2 (row 1 column 2)
            skew_21_coeff = new SumCoefficient(*d1u2_coeff, *d2u1_coeff, 0.5, -0.5); // d1u2-d2u1 (row 2 column 1)

            A11_coeff = new SumCoefficient(xi_coeff, *m0_coeff, 1.0, -1.0); 
            A12_coeff = new SumCoefficient(0., *skew_12_coeff, 0.0, -1.0); 
            A21_coeff = new SumCoefficient(0., *skew_21_coeff, 0.0, -1.0); 
            A22_coeff = new SumCoefficient(xi_coeff, *m0_coeff, 1.0, -1.0); 
        #endif

        a_a_2_chi_coeff = new ProductCoefficient(2 * a * a, chi_coeff); 
    }

    void calculate_operators(){

        fill_coefficient_matrix(*coeff_matrix, A11_coeff, A12_coeff, A21_coeff, A22_coeff, a_a_2_chi_coeff); 
    
        u_gf_NS->GetDerivative(1,0,*dxu1_gf); 
        u_gf_NS->GetDerivative(1,1,*dyu1_gf); 
        u_gf_NS->GetDerivative(2,0,*dxu2_gf); 
        u_gf_NS->GetDerivative(2,1,*dyu2_gf); 

        for(int i = 0; i < vector_size; i++){
            for(int j = 0; j < vector_size; j++){
                if((*A_entries)[i][j]){

                    ParBilinearForm a(&fespace);  
                    a.AddDomainIntegrator(new MassIntegrator(*((*coeff_matrix)[i][j]))); 
                    a.Assemble();
                    
                    (*A_SpMat)[i][j] = a.SpMat(); 

                    (*L_SpMat)[i][j] = a.SpMat();
                    (*L_SpMat)[i][j] *= - theta; 

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
        css_solver.SetMaxIter(2000);
        css_solver.SetPrintLevel(0);
        css_solver.SetOperator(*L_BO); 
    }

    void set_theta(double theta_){
        theta = theta_; 
        calculate_operators(); 
    }

    void reset_theta(){
        theta = get_theta(alpha, dt); 
        calculate_operators(); 
    }

    ~CSS(){}
};

