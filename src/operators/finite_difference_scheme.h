/*
    The configuration space operator is transformed to a finite difference scheme 
    by using a spectral decomposition based on the orthogonal basis of 
    weighted Hermite polynomials

    The finite difference scheme can be written as 
    d/dt phi = A phi 

    This approach is outlined in: 
    H. Mizerová, B. She, A conservative scheme for the Fokker-Planck equation with applications to
    viscoelastic polymeric fluids, J. Comput. Phys. (2018), https://doi.org/10.1016/j.jcp.2018.08.015

    The following classes and functions provide the functionality to create the coefficient matrix A 
*/

void fill_A_entries(std::vector<std::vector<bool>> &A_entries){
    int vector_size = A_entries.size(); 
    int N = sqrt(vector_size); 

    for(int i = 0; i < vector_size; i++){
        for(int j = 0; j < vector_size; j++){
            A_entries[i][j] = false; 
        }
    }

    for(int i = 0; i < vector_size; i++){
        int z = i / N; 
        int k = i % N; 

        // phi z,k 
        A_entries[i][i] = true; 
      
        // phi z-1,k-1 
        if(z-1 >= 0 && k-1 >=0){
            A_entries[i][(z-1) * N + k-1] = true; 
        }

        // phi z-1 k+1 
        if(z-1 >= 0 && k+1 < N){
            A_entries[i][(z-1) * N + k+1] = true; 
        }

        // phi z+1 k-1 
        if(z+1 < N && k-1 >= 0){
            A_entries[i][(z+1) * N + k-1] = true; 
        }
      
        // phi z-2 k 
        if(z-2 >= 0){
            A_entries[i][(z-2) * N + k] = true; 
        }

        // phi z   k-2
        if(k-2 >= 0){
            A_entries[i][(z) * N + k - 2 ] = true; 
        }
    } 
}

void fill_coefficient_matrix(std::vector<std::vector<Coefficient*>> &coeff_matrix, 
        SumCoefficient* A11_coeff, SumCoefficient* A12_coeff, SumCoefficient* A21_coeff, SumCoefficient* A22_coeff, 
        ProductCoefficient* alpha_alpha_2_chi_coeff){

    int vector_size = coeff_matrix.size(); 
    int N = sqrt(vector_size); 

    for(int i = 0; i < vector_size; i++){
        int z = i / N; 
        int k = i % N; 

        // phi z,k 
        coeff_matrix[i][i] = new SumCoefficient(*A11_coeff, *A22_coeff, (double) - z, (double) - k);
      
        // phi z-1,k-1 
        if(z-1 >= 0 && k-1 >=0){
            coeff_matrix[i][(z-1) * N + k-1] = new SumCoefficient(*A12_coeff, *A21_coeff, - sqrt((double) z * k), - sqrt((double) z * k));  
        }

        // phi z-1 k+1 
        if(z-1 >= 0 && k+1 < N){
            coeff_matrix[i][(z-1) * N + k+1] = new ProductCoefficient(- sqrt((double) z * (k+1)), *A12_coeff); 
        }

        // phi z+1 k-1 
        if(z+1 < N && k-1 >= 0){
            coeff_matrix[i][(z+1) * N + k-1] = new ProductCoefficient(- sqrt((double) (z+1) * k), *A21_coeff);
        }
      
        // phi z-2 k 
        if(z-2 >= 0){
            coeff_matrix[i][(z - 2) * N + k] = new SumCoefficient(*alpha_alpha_2_chi_coeff, *A11_coeff, sqrt((double) z * (z-1)), - sqrt((double) z * (z-1))); 
        }

        // phi z   k-2
        if(k-2 >= 0){
            coeff_matrix[i][(z) * N + k - 2] = new SumCoefficient(*alpha_alpha_2_chi_coeff, *A22_coeff, sqrt((double) k * (k-1)), - sqrt((double) k * (k-1)));  
        }
    } 
}