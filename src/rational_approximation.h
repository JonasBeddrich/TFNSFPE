using namespace std; 
using namespace mfem;

#if defined(alpha_05)
    const std::string tf_degree = "alpha=0.5"; 
    std::vector<double> get_weights(){
        std::vector<double> weights(n_modes,0.); 
        weights[0] = 0.26643085402672834;
        weights[1] = 0.36221216226480507;
        weights[2] = 0.5984688438847446;
        weights[3] = 1.0961683250350185;
        weights[4] = 2.1408864174910995;
        weights[5] = 4.294576741770443;
        weights[6] = 8.638825902824571;
        weights[7] = 17.70998953275846;
        weights[8] = 40.73670520684686;
        weights[9] = 151.12709891305988;
        return weights; 
    }

    std::vector<double> get_lambdas(){
        std::vector<double> lambdas(n_modes,0.); 
        lambdas[0] = 0.04140430416394466;
        lambdas[1] = 0.46642099889240357;
        lambdas[2] = 1.9988411809027893;
        lambdas[3] = 7.274221042614074;
        lambdas[4] = 26.412593428296983;
        lambdas[5] = 99.91058648562047;
        lambdas[6] = 390.35715107132563;
        lambdas[7] = 1562.7061487506387;
        lambdas[8] = 6755.098291449528;
        lambdas[9] = 42052.445913844596;
        return lambdas; 
    }

    double get_w_infinity(){
        return 0.00143131; 
    }
#endif 

#if defined(alpha_07)
    const std::string tf_degree = "alpha=0.7"; 
    std::vector<double> get_weights(){
        std::vector<double> weights(n_modes,0.); 
        weights[0] = 0.44743748480718115;
        weights[1] = 0.308113021726503;
        weights[2] = 0.3545945272231849;
        weights[3] = 0.48144558799093873;
        weights[4] = 0.7115751223978415;
        weights[5] = 1.0972380466887133;
        weights[6] = 1.7346619531188825;
        weights[7] = 2.8650307913290938;
        weights[8] = 5.149138597672265;
        weights[9] = 12.496583384693299;
        return weights; 
    }

    std::vector<double> get_lambdas(){
        std::vector<double> lambdas(n_modes,0.); 
        lambdas[0] = 0.018511135047264583;
        lambdas[1] = 0.3052947599678512;
        lambdas[2] = 1.2859253208373544;
        lambdas[3] = 4.375919634103344;
        lambdas[4] = 14.70318355356868;
        lambdas[5] = 52.16602677281238;
        lambdas[6] = 198.79227495587386;
        lambdas[7] = 827.1563153983546;
        lambdas[8] = 3952.9489514082443;
        lambdas[9] = 26263.791754441896;
        return lambdas; 
    }

    double get_w_infinity(){
        return 9.72339715e-05; 
    }
#endif 

#if defined(alpha_09)
    const std::string tf_degree = "alpha=0.9"; 
    std::vector<double> get_weights(){
        std::vector<double> weights(n_modes,0.); 
        weights[0] = 0.7804307509939867;
        weights[1] = 0.16004349367916604;
        weights[2] = 0.1332151343665852;
        weights[3] = 0.1399438673526587;
        weights[4] = 0.1608443612802441;
        weights[5] = 0.19033855280506712;
        weights[6] = 0.2269452182080697;
        weights[7] = 0.27545767838709895;
        weights[8] = 0.3538508618735971;
        weights[9] = 0.5806252904038118;
        return weights; 
    }

    std::vector<double> get_lambdas(){
        std::vector<double> lambdas(n_modes,0.); 
        lambdas[0] = 0.006179326206527797;
        lambdas[1] = 0.2786110317351769;
        lambdas[2] = 1.2692161602466712;
        lambdas[3] = 4.465273909609429;
        lambdas[4] = 15.273695030857832;
        lambdas[5] = 54.41622330742792;
        lambdas[6] = 204.85835562154304;
        lambdas[7] = 822.863712357904;
        lambdas[8] = 3668.633075979774;
        lambdas[9] = 22184.82014296115;
        return lambdas; 
    }

    double get_w_infinity(){
        return 3.22631261e-06; 
    }
#endif 

#if defined(alpha_1)
    const std::string tf_degree = "alpha=1.0"; 
    std::vector<double> get_weights(){
        std::vector<double> weights(n_modes,0.); 
        weights[0] = 1; 
        return weights; 
    }

    std::vector<double> get_lambdas(){
        std::vector<double> lambdas(n_modes,0.); 
        return lambdas; 
    }

    double get_w_infinity(){
        return 0; 
    }
#endif 

std::vector<double> get_gammas(double dt){
    std::vector<double> gammas(get_lambdas()); 
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::multiplies<double>(), std::placeholders::_1, dt));
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::plus<double>(), std::placeholders::_1, 1));
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::divides<double>(), 1, std::placeholders::_1));
    return gammas; 
} 

double get_beta(double dt){
    double beta = 0; 
    std::vector<double> gammas = get_gammas(dt); 
    std::vector<double> weights = get_weights(); 
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::multiplies<double>(), std::placeholders::_1, dt));
    std::transform(gammas.begin(), gammas.end(), weights.begin(), gammas.begin(), std::multiplies<double>()); 
    std::for_each(gammas.begin(), gammas.end(), [&] (double d) {beta += d;}); 
    return beta + get_w_infinity(); 
}

double get_delta(double dt){
    double delta = 0; 
    std::vector<double> weights = get_weights(); 
    std::for_each(weights.begin(), weights.end(), [&] (double d) {delta += d;}); 
    return delta; 
}