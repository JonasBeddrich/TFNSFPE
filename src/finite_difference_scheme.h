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

// The Coefficient Factory provides how phi_zk depends on other phis 
// It represents one row in the matrix A 
class CoefficientFactory {

private: 
    double alpha, xi, chi, eps; 
    int z,k; 

    // the following notation is taken from Mizerová & She 2018 
    double A11(const Vector &x){
        return - du1dx1(x) + xi; 
    }

    double A12(const Vector &x){
        return - du1dx2(x); 
    }

    double A21(const Vector &x){
        return - du2dx1(x); 
    }

    double A22(const Vector &x){
        return - du2dx2(x) + xi; 
    }

public:
    CoefficientFactory(){}
    
    CoefficientFactory(double alpha_, double xi_, double chi_, double eps_, int z_, int k_):
        alpha(alpha_), 
        xi(xi_), 
        chi(chi_), 
        eps(eps_), 
        z(z_), 
        k(k_){}


    // z,k 
    double Eval_z_k(const Vector& x, double t){ 
        return (- A11(x) * z - A22(x) * k); 
    }

    // z-1,k-1
    double Eval_zm1_km1(const Vector& x, double t){
        return (- A12(x) - A21(x)) * sqrt(z * k); 
    }

    // z-1,k+1
    double Eval_zm1_kp1(const Vector& x, double t){
        return (- A12(x)) * sqrt(z * (k+1)); 
    }

    // z+1,k-1
    double Eval_zp1_km1(const Vector& x, double t){
        return (- A21(x)) * sqrt((z+1) * k); 
    }

    // z-2,k
    double Eval_zm2_k(const Vector& x, double t){
        return (2 * alpha * alpha * chi - A11(x)) * sqrt(z * (z-1)); 
    }

    // z,k-2
    double Eval_z_km2(const Vector& x, double t){
        return (2 * alpha * alpha * chi - A22(x)) * sqrt(k * (k-1)); 
    } 
}; 

// This function fills a vector with CoefficientFactory objects 
void fill_CoefficientVector(std::vector<CoefficientFactory> &CoefficientVector){
    int vector_size = CoefficientVector.size(); 
    int N = sqrt(vector_size); 

    for (int i = 0; i < vector_size; i++){
        int z = i / N; 
        int k = i % N; 
        CoefficientVector[i] = CoefficientFactory(alpha, xi, chi, eps, z, k); 
    }  
}


// This function fills the matrix A with FunctionCoefficients 
// based on the CoefficientFactory vector 
// It also fills A_entries with ones where A has an entry 
void fill_A(std::vector<std::vector<FunctionCoefficient>> &A, std::vector<std::vector<bool>> & A_entries, std::vector<CoefficientFactory> &CoefficientVector){
    
    int vector_size = A.size(); 
    int N = sqrt(vector_size); 

    for(int i = 0; i < vector_size; i++){
        for(int j = 0; j < vector_size; j++){
            A_entries[i][j] = false; 
        }
    }

    // The iteration starts with 1 since we ignore phi_00
    // Thus indices are considered -1  
    for(int i = 0; i < vector_size; i++){
        int z = i / N; 
        int k = i % N; 
      
        // phi z   k 
        auto eval_z_k = std::bind(&CoefficientFactory::Eval_z_k, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
        A[i][i] = FunctionCoefficient(eval_z_k); 
        A_entries[i][i] = true; 
      
        // phi z-1 k-1 
        if(z-1 >= 0 && k-1 >=0){
            auto eval_zm1_km1 = std::bind(&CoefficientFactory::Eval_zm1_km1, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
            A[i][(z-1) * N + k-1] = FunctionCoefficient(eval_zm1_km1); 
            A_entries[i][(z-1) * N + k-1] = true; 
        }

        // phi z-1 k+1 
        if(z-1 >= 0 && k+1 < N){
            auto eval_zm1_kp1 = std::bind(&CoefficientFactory::Eval_zm1_kp1, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
            A[i][(z-1) * N + k+1] = FunctionCoefficient(eval_zm1_kp1); 
            A_entries[i][(z-1) * N + k+1] = true; 
        }

        // phi z+1 k-1 
        if(z+1 < N && k-1 >= 0){
            auto eval_zp1_km1 = std::bind(&CoefficientFactory::Eval_zp1_km1, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
            A[i][(z+1) * N + k-1] = FunctionCoefficient(eval_zp1_km1); 
            A_entries[i][(z+1) * N + k-1] = true; 
        }
      
        // phi z-2 k 
        if(z-2 >= 0){
            auto eval_zm2_k = std::bind(&CoefficientFactory::Eval_zm2_k, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
            A[i][(z-2) * N + k] = FunctionCoefficient(eval_zm2_k); 
            A_entries[i][(z-2) * N + k] = true; 
        }

        // phi z   k-2
        if(k-2 >= 0){
            auto eval_z_km2 = std::bind(&CoefficientFactory::Eval_z_km2, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
            A[i][(z) * N + k - 2 ] = FunctionCoefficient(eval_z_km2); 
            A_entries[i][(z) * N + k - 2 ] = true; 
        }
    } 
}

void fill_coefficient_matrix(std::vector<std::vector<Coefficient*>> &coeff_matrix, std::vector<std::vector<bool>> &A_entries, 
        SumCoefficient* A11_coeff, SumCoefficient* A12_coeff, SumCoefficient* A21_coeff, SumCoefficient* A22_coeff, 
        ProductCoefficient* alpha_alpha_2_xi_coeff){
    int vector_size = coeff_matrix.size(); 
    int N = sqrt(vector_size); 

    for(int i = 0; i < vector_size; i++){
        for(int j = 0; j < vector_size; j++){
            A_entries[i][j] = false; 
        }
    }

    // The iteration starts with 1 since we ignore phi_00
    // Thus indices are considered -1  
    for(int i = 0; i < vector_size; i++){
        int z = i / N; 
        int k = i % N; 

        // phi z,k 
        coeff_matrix[i][i] = new SumCoefficient(*A11_coeff, *A22_coeff, (double) -z, (double) -k);  
        A_entries[i][i] = true; 
      
        // phi z-1,k-1 
        if(z-1 >= 0 && k-1 >=0){
            coeff_matrix[i][(z-1) * N + k-1] = new SumCoefficient(*A12_coeff, *A21_coeff, -1.0, -sqrt((double) z * k));  
            A_entries[i][(z-1) * N + k-1] = true; 
        }

        // phi z-1 k+1 
        if(z-1 >= 0 && k+1 < N){
            coeff_matrix[i][(z-1) * N + k+1] = new ProductCoefficient(- sqrt((double) z * (k+1)), *A12_coeff); 
            A_entries[i][(z-1) * N + k+1] = true; 
        }

        // phi z+1 k-1 
        if(z+1 < N && k-1 >= 0){
            coeff_matrix[i][(z+1) * N + k-1] = new ProductCoefficient(- sqrt((double) (z+1) * k), *A21_coeff);  
            A_entries[i][(z+1) * N + k-1] = true; 
        }
      
        // phi z-2 k 
        if(z-2 >= 0){
            coeff_matrix[i][(z-2) * N + k] = new SumCoefficient(*alpha_alpha_2_xi_coeff, *A11_coeff, sqrt((double) z * (z-1)), - sqrt((double) z * (z-1))); 
            A_entries[i][(z-2) * N + k] = true; 
        }

        // phi z   k-2
        if(k-2 >= 0){
            coeff_matrix[i][(z) * N + k - 2 ] = new SumCoefficient(*alpha_alpha_2_xi_coeff, *A22_coeff, sqrt((double) k * (k-1)), - sqrt((double) k * (k-1)));  
            A_entries[i][(z) * N + k - 2 ] = true; 
        }
    } 
}