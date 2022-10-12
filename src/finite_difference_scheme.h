/*

    This header files includes everything for the finite difference scheme of the configuration solver

*/

// The coefficient factory provides the functions to relate phi_z,k to all other relevant phis 
class CoefficientFactory {
  private: 
    double alpha, xi, chi, eps; 
    int z,k; 
  
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

/* 
    Create a coefficient factory vector of length N  
    i.e. one for every phi_z,k
*/
void fill_CoefficientVector(std::vector<CoefficientFactory> &CoefficientVector){
  int vector_size = CoefficientVector.size(); 
  int N = sqrt(vector_size); 

  for (int i = 0; i < vector_size; i++){
    int z = i / N; 
    int k = i % N; 
    CoefficientVector[i] = CoefficientFactory(alpha, xi, chi, eps, z, k); 
  }  
}

/*
    This function fills a matrix A with the FunctionCoefficients
    A_entries is true when there is an entry in A 
*/
void fill_A(std::vector<std::vector<FunctionCoefficient>> &A, std::vector<std::vector<bool>> & A_entries, std::vector<CoefficientFactory> &CoefficientVector){
  int vector_size_partial = A.size(); 
  int vector_size = vector_size_partial + 1; 
  int N = sqrt(vector_size); 

  for(int i = 0; i < vector_size_partial; i++){
    for(int j = 0; j < vector_size_partial; j++){
      A_entries[i][j] = false; 
    }
  }

  // The iteration starts with 1 since we ignore phi_00
  // Thus indices are considered -1  
  for(int i = 1; i < vector_size; i++){
    int z = i / N; 
    int k = i % N; 
    
    // phi z   k 
    auto eval_z_k = std::bind(&CoefficientFactory::Eval_z_k, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
    A[i-1][i-1] = FunctionCoefficient(eval_z_k); 
    A_entries[i-1][i-1] = true; 
    
    // phi z-1 k-1 
    // the case z=k=1 has to be ignored 
    if(z-1 >= 0 && k-1 >=0 && (z != 1 || k != 1)){
      auto eval_zm1_km1 = std::bind(&CoefficientFactory::Eval_zm1_km1, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
      A[i-1][(z-1) * N + k-1-1] = FunctionCoefficient(eval_zm1_km1); 
      A_entries[i-1][(z-1) * N + k-1-1] = true; 
    }

    // phi z-1 k+1 
    if(z-1 >= 0 && k+1 < N){
      auto eval_zm1_kp1 = std::bind(&CoefficientFactory::Eval_zm1_kp1, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
      A[i-1][(z-1) * N + k+1-1] = FunctionCoefficient(eval_zm1_kp1); 
      A_entries[i-1][(z-1) * N + k+1-1] = true; 
    }

    // phi z+1 k-1 
    if(z+1 < N && k-1 >= 0){
      auto eval_zp1_km1 = std::bind(&CoefficientFactory::Eval_zp1_km1, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
      A[i-1][(z+1) * N + k-1-1] = FunctionCoefficient(eval_zp1_km1); 
      A_entries[i-1][(z+1) * N + k-1-1] = true; 
    }
    
    // phi z-2 k 
    // the case z=2, k=0 has to be ignored 
    if(z-2 >= 0 && (z != 2 || k != 0)){
      auto eval_zm2_k = std::bind(&CoefficientFactory::Eval_zm2_k, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
      A[i-1][(z-2) * N + k-1] = FunctionCoefficient(eval_zm2_k); 
      A_entries[i-1][(z-2) * N + k-1] = true; 
    }

    // phi z   k-2
    // the case z=0, k=2 has to be ignored 
    if(k-2 >= 0 && (z != 0 || k != 2)){
      auto eval_z_km2 = std::bind(&CoefficientFactory::Eval_z_km2, CoefficientVector[i], std::placeholders::_1, std::placeholders::_2); 
      A[i-1][(z) * N + k - 2 - 1] = FunctionCoefficient(eval_z_km2); 
      A_entries[i-1][(z) * N + k - 2 - 1] = true; 
    }
  } 
}