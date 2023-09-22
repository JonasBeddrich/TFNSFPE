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

    ParGridFunction *u_gf_NS; 
    ParGridFunction *dxu1_gf, *dyu1_gf, *dxu2_gf, *dyu2_gf;
    
    std::vector<std::vector<bool>> *A_entries;

    BlockOperator M_BO; 
    BlockOperator A_BO; 
    BlockOperator L_BO; 
    BlockOperator *A44;
    BlockOperator *A42;

    HypreParMatrix m_HPM; 
    HypreParMatrix m0_HPM;  

    std::vector<std::vector<HypreParMatrix*>> *A_ptrs; 
    std::vector<std::vector<HypreParMatrix*>> *L_ptrs; 

    ProductCoefficient &chi_coeff, &xi_coeff; 
    GridFunctionCoefficient *d1u1_coeff, *d1u2_coeff, *d2u1_coeff, *d2u2_coeff; 
    SumCoefficient *A11_coeff, *A12_coeff, *A21_coeff, *A22_coeff, *du_sym_coeff, *skew_12_coeff, *skew_21_coeff; 
    ProductCoefficient *a_a_2_chi_coeff; 
    ConstantCoefficient *m0_coeff; 

    BiCGSTABSolver css_solver;
    BiCGSTABSolver A44_solver;
    mutable Vector z; 

public:

    CSS(ParFiniteElementSpace &fespace_, int vector_size_, const Array<int> &offsets_, 
    ParGridFunction *u_gf_NS_, ProductCoefficient &chi_coeff_, ProductCoefficient &xi_coeff_, 
    double theta_): 
        z(fespace_.GetTrueVSize() * vector_size_), 
        fespace(fespace_),  
        vector_size(vector_size_),
        n_dof(fespace_.GetTrueVSize()), 
        N(sqrt(vector_size_)), 
        u_gf_NS(u_gf_NS_), 
        chi_coeff(chi_coeff_), 
        xi_coeff(xi_coeff_), 
        css_solver(MPI_COMM_WORLD), 
        M_BO(offsets_), 
        L_BO(offsets_), 
        A_BO(offsets_), 
        theta(theta_){

        setup_coefficients();   
        
        A_entries = new std::vector<std::vector<bool>>(vector_size, std::vector<bool>(vector_size)); 
        fill_A_entries(*A_entries); 

        HypreParMatrix *empty_HPM = new HypreParMatrix(); 
        A_ptrs = new std::vector<std::vector<HypreParMatrix*>>(vector_size, std::vector(vector_size, empty_HPM)); 
        L_ptrs = new std::vector<std::vector<HypreParMatrix*>>(vector_size, std::vector(vector_size, empty_HPM)); 
                    
        calculate_operators(); 
    }

    void apply_FR(const Vector &phi, Vector &A_phi) const{
        // setup_block_operators(); 
        A_BO.Mult(phi,A_phi); 
    }

    void solve_Id_minus_theta_FR2(const Vector &y, Vector &x){
        // calculates x = (Id - theta FR)^-1 M*y 
                
        // z = M*y 
        M_BO.Mult(y, z); 

        // x = (Id - theta FR)^-1 z
        css_solver.Mult(z,x); 
    }

   
    void setup_coefficients(){

        // theta = get_theta(alpha, dt); 

        ParBilinearForm m(&fespace);
        m.AddDomainIntegrator(new MassIntegrator); 
        m.Assemble(); 
        m.Finalize(); 
        m_HPM = *(m.ParallelAssemble()); 

        for(int i = 0; i < vector_size; i++){
            M_BO.SetBlock(i, i, &m_HPM); 
        }
    
        m0_coeff = new ConstantCoefficient(0);
        ParBilinearForm m0(&fespace);
        m0.AddDomainIntegrator(new MassIntegrator(*m0_coeff)); 
        m0.Assemble(); 
        m0.Finalize(); 
        m0_HPM = *(m0.ParallelAssemble()); 

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
    
        u_gf_NS->GetDerivative(1,0,*dxu1_gf); 
        u_gf_NS->GetDerivative(1,1,*dyu1_gf); 
        u_gf_NS->GetDerivative(2,0,*dxu2_gf); 
        u_gf_NS->GetDerivative(2,1,*dyu2_gf); 

        std::vector<std::vector<Coefficient*>> coeff_matrix(vector_size, std::vector<Coefficient*>(vector_size,NULL));  
        fill_coefficient_matrix(coeff_matrix, A11_coeff, A12_coeff, A21_coeff, A22_coeff, a_a_2_chi_coeff); 

        for(int i = 0; i < vector_size; i++){
            for(int j = 0; j < vector_size; j++){
                if((*A_entries)[i][j]){

                    ParBilinearForm a(&fespace);  
                    a.AddDomainIntegrator(new MassIntegrator(*(coeff_matrix[i][j]))); 
                    a.Assemble();
                    a.Finalize();

                    // ParBilinearForm l(&fespace); 
                    // double diag = 0; 
                    // if(i==j){
                    //     diag = 1; 
                    // }
                    // SumCoefficient l_coeff(diag, *(coeff_matrix[i][j]), 1, -theta); 
                    // l.AddDomainIntegrator(new MassIntegrator(l_coeff)); 
                    // l.Assemble();
                    // l.Finalize(); 
                                        
                    // Keeping track of the matrices 
                    (*A_ptrs)[i][j] = a.ParallelAssemble();  
                    // (*L_ptrs)[i][j] = l.ParallelAssemble();  

                    if(i == j){
                        (*L_ptrs)[i][j] = new HypreParMatrix(m_HPM); 
                        (*L_ptrs)[i][j]->Add(-theta, *(*A_ptrs)[i][j]);  
                    } else {
                        (*L_ptrs)[i][j] = new HypreParMatrix(*(*A_ptrs)[i][j]); 
                        *((*L_ptrs)[i][j]) *= -theta;    
                    }
                  
                    A_BO.SetBlock(i, j, (*A_ptrs)[i][j]); 
                    L_BO.SetBlock(i, j, (*L_ptrs)[i][j]); 

                    delete coeff_matrix[i][j];

                    // A_SpMat[i][j] = *(a.ParallelAssemble());  
                    // L_SpMat[i][j] = *(l.ParallelAssemble());  
                                        
                    // if(i == j){
                    //     (*L_SpMat)[i][j] = m_HPM; // not the problem
                    //     B_SpMat[i][j] = m_HPM; 

                    //     (*L_SpMat)[i][j].Add(-theta, *(a.ParallelAssemble())); 
                    // } else {
                    //     (*L_SpMat)[i][j] = *(a.ParallelAssemble()); 
                    //     (*L_SpMat)[i][j] *= -theta;   
                    // }
                }
            }
        }

        css_solver.iterative_mode = false;
        css_solver.SetRelTol(1e-16);
        css_solver.SetAbsTol(1e-12); 
        css_solver.SetMaxIter(10000);
        css_solver.SetPrintLevel(0);
        css_solver.SetOperator(L_BO); 
    }

    void set_theta(double theta_){
        theta = theta_; 
        calculate_operators(); 
    }

    void reset_theta(){
        theta = get_theta(alpha, dt); 
        calculate_operators(); 
    }

    ~CSS(){                

        for(int i = 0; i < vector_size; i++){
            for(int j = 0; j < vector_size; j++){
                if((*A_entries)[i][j]){
                    delete (*A_ptrs)[i][j];
                    delete (*L_ptrs)[i][j]; 
                    
                }
            }
        }

        delete A_ptrs; 
        delete L_ptrs; 

        delete A_entries; 

        delete u_fes; 
        delete dxu1_gf; 
        delete dyu1_gf; 
        delete dxu2_gf; 
        delete dyu2_gf; 

        delete d1u1_coeff; 
        delete d1u2_coeff; 
        delete d2u1_coeff; 
        delete d2u2_coeff; 
        
        delete A11_coeff; 
        delete A12_coeff; 
        delete A21_coeff; 
        delete A22_coeff; 
        
        #if defined(symmetric)
            delete du_sym_coeff;
        #endif

        #if defined(skew_symmetric)
            delete skew_12_coeff; 
            delete skew_21_coeff; 
        #endif

        delete a_a_2_chi_coeff; 
        delete m0_coeff; 
    }

     void solve_Id_minus_theta_FR(BlockVector &y, BlockVector &x){
        A44_solver.iterative_mode = false;
        A44_solver.SetRelTol(1e-16);
        A44_solver.SetMaxIter(2000);
        A44_solver.SetPrintLevel(0);
        
        ParBilinearForm m(&fespace);
        m.AddDomainIntegrator(new MassIntegrator); 
        m.Assemble(); 
        m.Finalize(); 
        HypreParMatrix *m_HPM = m.ParallelAssemble(); 

        // ****************************************************************
        // The following code is concerned with solving A_4,4 phi = A_4,2 mu 

        for (int counter = 1; counter < 2*N; counter++){ // loop over the diagonals 
        // for (int counter = 1; counter < 3; counter++){ 
            if(counter < N+1){
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

                BlockOperator A44(part_offsets);
                BlockOperator A42(part_offsets);

                // Load mu k 
                for(int z=0; z < counter; z++){
                    int k = (counter-1) - z; 
                    m_HPM->Mult(y.GetBlock(z*N + k), mu_k.GetBlock(z)); 
                }
                
                // Load phi km2 
                for(int z=0; z < counter-2; z++){
                    int k = (counter-1) -2 -z; 
                    phi_km2.GetBlock(z) = x.GetBlock(z*N + k); 
                }

                // Load A_4,2
                for(int z_i=0; z_i < counter; z_i++){ // loop rows 
                    for(int z_j=0; z_j < counter-2; z_j++){ // loop columns, where the indices add up to z+k-2 
                        int k_i = (counter-1) -z_i; 
                        int k_j = (counter-1) -2 -z_j;   
                        if((*A_entries)[z_i * N + k_i][z_j * N + k_j]){
                            // A42.SetBlock(z_i,z_j,&(L_BO.GetBlock(z_i * N + k_i,z_j * N + k_j))); 
                            A42.SetBlock(z_i,z_j,(*L_ptrs)[z_i * N + k_i][z_j * N + k_j]); 
                        } 
                    }
                }

                // Compute A_4,2 phi 
                A42.Mult(phi_km2, A42_phi_km2); 

                mu_k -= A42_phi_km2; 

                // Load A_4,4
                for(int z_i=0; z_i < counter; z_i++){ // loop rows 
                    for(int z_j=0; z_j < counter; z_j++){ // loop columns 
                        int k_i = (counter-1) -z_i; 
                        int k_j = (counter-1) -z_j; 

                        if((*A_entries)[z_i * N + k_i][z_j * N + k_j]){
                            // A44.SetBlock(z_i,z_j,&(L_BO.GetBlock(z_i * N + k_i,z_j * N + k_j))); 
                            A44.SetBlock(z_i, z_j, (*L_ptrs)[z_i * N + k_i][z_j * N + k_j]); 
                            // cout << Mpi::WorldRank() << "-" << n_dof << "-" << A44.RowOffsets()[1] <<endl;
                        } 
                    }
                }

                // Solve A44 ... 
                A44_solver.SetOperator(A44); 
                A44_solver.Mult(mu_k, phi_k); 

                // Write phi km2 
                for(int z=0; z < counter; z++){
                    int k = (counter-1) -z; 
                    x.GetBlock(z*N + k) = phi_k.GetBlock(z); 
                }
            } 
                
            else if (counter == N+1){
                // cout << "else if" << endl; 
                // Setup 
                Array<int> part_offsets(N);
                part_offsets[0] = 0; 
                for (int i = 1; i < N; i++){
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
                for(int i=0; i < N-1; i++){
                    int z = i +1; 
                    int k = N -z; 
                    // cout << N << " " << i << " " << z << " " << k << " " << z*N + k << endl;  
                    m_HPM->Mult(y.GetBlock(z*N + k), mu_k.GetBlock(i)); 
                }

                // Load phi km2 
                for(int i=0; i < N-1; i++){
                    int k = N -2 -i; 
                    int z = N -2 -k; 
                    // cout << i << " " << z << " " << k << " " << z*N + k << endl;  
                    phi_km2.GetBlock(z) = x.GetBlock(z*N + k); 
                }

                // Load A_k,km2
                for(int i=0; i < N-1; i++){ // loop rows 
                    for(int j=0; j < N-1; j++){ // loop columns 
                        int z_i = i +1; 
                        int k_i = N -z_i;  
                        
                        int k_j = N -2 -j; 
                        int z_j = N -2 -k_j; 

                        // cout << "i " << i << " " << z_i << " " << k_i << endl; 
                        // cout << "j " << j << " " << z_j << " " << k_j << endl; 
                        // cout << endl; 
                        if((*A_entries)[z_i * N + k_i][z_j * N + k_j]){
                            // cout << "A entries is true: " << z_i << " " << k_i << " " << z_j << " " << k_j << endl; 
                            A42->SetBlock(i,j,&(L_BO.GetBlock(z_i * N + k_i,z_j * N + k_j))); 
                        } 
                    }
                }

                // Compute A_4,2 phi 
                A42->Mult(phi_km2, A42_phi_km2); 
                mu_k -= A42_phi_km2; 
                // TODO: I think this should be a minus 
                
                // Load A_4,4
                for(int i=0; i < N-1; i++){ // loop rows 
                    for(int j=0; j < N-1; j++){ // loop columns 
                        int z_i = i +1; 
                        int k_i = N -z_i; 
                        
                        int z_j = j +1; 
                        int k_j = N -z_j; 
                        
                        // cout << "i " << i << " " << z_i << " " << k_i << endl; 
                        // cout << "j " << j << " " << z_j << " " << k_j << endl; 
                        // cout << endl; 

                        if((*A_entries)[z_i * N + k_i][z_j * N + k_j]){
                            A44->SetBlock(i,j,&(L_BO.GetBlock(z_i * N + k_i,z_j * N + k_j))); 
                        } 
                    }
                }

                // Solve A44 ... 
                A44_solver.SetOperator(*A44); 
                A44_solver.Mult(mu_k, phi_k); 

                // Write phi_k 
                for(int i=0; i < N-1; i++){
                    int z = i +1; 
                    int k = N -z; 
                    // for (int j = 0; j < n_dof; j++){
                        // cout << phi_k.GetBlock(i)[j] << " " << x.GetBlock(z*N + k)[j] << endl; 
                    // }
                    x.GetBlock(z*N + k) = phi_k.GetBlock(i); 
                    // for (int j = 0; j < n_dof; j++){
                    //     cout << phi_k.GetBlock(i)[j] << " " << x.GetBlock(z*N + k)[j] << endl; 
                    // }
                    // cout << z << " " << k << " " << z*N + k << endl;  
                }

                // cout << "phi_k" << endl; 
                // for(int i = 0; i < phi_k.Size() / 9; i++){
                //     cout << i; 
                //     for(int j = 0; j < 9; j++){
                //     cout << " - " << phi_k[i*9+j]; 
                //     }
                //     cout << endl; 
                // }                
            } 

            
            
            else { // counter > N+1 
                Array<int> short_offsets(2*N -counter +1);
                short_offsets[0] = 0; 
                for (int i = 1; i < 2*N -counter +1; i++){
                    short_offsets[i] = n_dof;
                }
                short_offsets.PartialSum();

                Array<int> long_offsets(2*N -counter +2 +1);
                long_offsets[0] = 0; 
                for (int i = 1; i < 2*N -counter +2 +1; i++){
                    long_offsets[i] = n_dof;
                }
                long_offsets.PartialSum();

                BlockVector mu_k(short_offsets); 
                BlockVector phi_k(short_offsets); 
                BlockVector phi_km2(long_offsets); 
                BlockVector A42_phi_km2(long_offsets); 

                A44 = new BlockOperator(short_offsets); 
                A42 = new BlockOperator(long_offsets); 
                
                // Load mu k 
                for(int i=0; i < 2*N -counter; i++){
                    int z = counter - N + i;  
                    int k = N -i -1; 
                    // cout << N << " " << i << " " << z << " " << k << " " << z*N + k << endl;  
                    m_HPM->Mult(y.GetBlock(z*N+k), mu_k.GetBlock(i)); 
                }
                
                // Load phi km2 
                for(int i=0; i < 2*N -counter +2; i++){
                    int z = counter -N -2 +i;                      
                    int k = N -i -1; 
                    // cout << i << " " << z << " " << k << " " << z*N + k << endl;  
                    phi_km2.GetBlock(i) = x.GetBlock(z*N + k); 
                }
                
                // Load A_k,km2
                for(int i=0; i < 2*N -counter; i++){ // loop rows 
                    for(int j=0; j < 2*N -counter +2; j++){ // loop columns             
                        int z_i = counter - N + i;  
                        int k_i = N -i -1; 
                        
                        int z_j = counter -N -2 +j;                      
                        int k_j = N -j -1;

                        // cout << "i " << i << ": " << z_i << " " << k_i << endl; 
                        // cout << "j " << j << ": " << z_j << " " << k_j << endl; 
                        // cout << (*A_entries)[z_i * N + k_i][z_j * N + k_j] << endl; 

                        if((*A_entries)[z_i * N + k_i][z_j * N + k_j]){
                            // cout << "A entries is true: " << z_i << " " << k_i << " " << z_j << " " << k_j << " " << endl; 
                            A42->SetBlock(i,j,&(L_BO.GetBlock(z_i * N + k_i,z_j * N + k_j))); 
                        } 
                    }
                }
                // cout << "mu_k" << endl; 
                // for(int i = 0; i < mu_k.Size() / 9; i++){
                //     cout << i; 
                //     for(int j = 0; j < 9; j++){
                //     cout << " - " << mu_k[i*9+j]; 
                //     }
                //     cout << endl; 
                // }      

                // Compute A_4,2 phi 
                A42->Mult(phi_km2, A42_phi_km2); 
                mu_k -= A42_phi_km2; 
                // TODO: I think this should be a minus 

                // cout << "A42_phi_km2" << endl; 
                // for(int i = 0; i < A42_phi_km2.Size() / 9; i++){
                //     cout << i; 
                //     for(int j = 0; j < 9; j++){
                //     cout << " - " << A42_phi_km2[i*9+j]; 
                //     }
                //     cout << endl; 
                // }         

                // cout << "mu_k-A42_phi_km2" << endl; 
                // for(int i = 0; i < mu_k.Size() / 9; i++){
                //     cout << i; 
                //     for(int j = 0; j < 9; j++){
                //     cout << " - " << mu_k[i*9+j]; 
                //     }
                //     cout << endl; 
                // }      


                // Load A_4,4
                for(int i=0; i < 2*N - counter; i++){ // loop rows 
                    for(int j=0; j < 2*N - counter; j++){ // loop columns 
                        int z_i = counter - N + i;  
                        int k_i = N -i -1; 
                        int z_j = counter - N + j;  
                        int k_j = N -j -1; 
                        
                        // cout << "i " << i << " " << z_i << " " << k_i << endl; 
                        // cout << "j " << j << " " << z_j << " " << k_j << endl; 
                        // cout << endl; 

                        if((*A_entries)[z_i * N + k_i][z_j * N + k_j]){
                            // cout << "A entries is true: " << z_i << " " << k_i << " " << z_j << " " << k_j << " " << endl; 
                            A44->SetBlock(i,j,&(L_BO.GetBlock(z_i * N + k_i,z_j * N + k_j))); 
                        }
                    }
                }
                 
                // Solve A44 ... 
                A44_solver.SetOperator(*A44); 
                A44_solver.Mult(mu_k, phi_k); 

                // Write phi_k 
                for(int i=0; i < 2*N - counter; i++){
                    int z = counter - N + i;  
                    int k = N -i -1; 
                    // for (int j = 0; j < n_dof; j++){
                        // cout << phi_k.GetBlock(i)[j] << " " << x.GetBlock(z*N + k)[j] << endl; 
                    // }
                    x.GetBlock(z*N + k) = phi_k.GetBlock(i); 
                    // for (int j = 0; j < n_dof; j++){
                    //     cout << phi_k.GetBlock(i)[j] << " " << x.GetBlock(z*N + k)[j] << endl; 
                    // }
                    // cout << z << " " << k << " " << z*N + k << endl;  
                }

                // cout << "phi_k" << endl; 
                // for(int i = 0; i < phi_k.Size() / 9; i++){
                //     cout << i; 
                //     for(int j = 0; j < 9; j++){
                //     cout << " - " << phi_k[i*9+j]; 
                //     }
                //     cout << endl; 
                // }           
            }

        }

        // cout << endl; 
        // cout << "x" << endl;
        // for(int i = 0; i < x.Size()/9; i++){
        //     cout << "i: " << i; 
        //     for(int j = 0; j < 9; j++){
        //         cout << " - " << x[i*9+j] ; 
        //     }
        //     cout << endl; 
        // }
        // cout << endl; 
        delete m_HPM; 
    }

};

