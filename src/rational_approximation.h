using namespace std; 
using namespace mfem;

#if defined(alpha_07)
    const std::string tf_degree = "alpha=0.7"; 
    std::vector<double> get_weights(){
        std::vector<double> weights(n_modes,0.); 
        weights[0] = 0.32173695029727833;
        weights[1] = 0.1929594378231144;
        weights[2] = 0.17601551293668669;
        weights[3] = 0.1821500169447637;
        weights[4] = 0.20149223702803218;
        weights[5] = 0.2317453634281737;
        weights[6] = 0.27293121959508965;
        weights[7] = 0.32627366365937627;
        weights[8] = 0.3940713976457411;
        weights[9] = 0.47996009994965677;
        weights[10] = 0.5891537753001249;
        weights[11] = 0.7288000983214853;
        weights[12] = 0.9090479214619074;
        weights[13] = 1.1454528026357915;
        weights[14] = 1.4639098576581222;
        weights[15] = 1.9121413560617315;
        weights[16] = 2.5927472729833196;
        weights[17] = 3.780871545584721;
        weights[18] = 6.475591275748234;
        weights[19] = 17.0959196037899;
        return weights; 
    }

    std::vector<double> get_lambdas(){
        std::vector<double> lambdas(n_modes,0.); 
        lambdas[0] = 0.006354300821057348;
        lambdas[1] = 0.09322685279869868;
        lambdas[2] = 0.30796780961156717;
        lambdas[3] = 0.7287673777465724;
        lambdas[4] = 1.5108887601104624;
        lambdas[5] = 2.948033881817028;
        lambdas[6] = 5.591971232914399;
        lambdas[7] = 10.483702390261225;
        lambdas[8] = 19.604313309288244;
        lambdas[9] = 36.77351305737705;
        lambdas[10] = 69.47879657914815;
        lambdas[11] = 132.66495088879452;
        lambdas[12] = 256.77669360574714;
        lambdas[13] = 505.4799927627441;
        lambdas[14] = 1016.8771249656868;
        lambdas[15] = 2107.1703954283744;
        lambdas[16] = 4565.325984299194;
        lambdas[17] = 10680.064047162616;
        lambdas[18] = 29316.045231605534;
        lambdas[19] = 125893.40106678547;
        return lambdas; 
    }

    double get_w_infinity(){
        return 3.70179209e-05; 
    }
#endif 

#if defined(alpha_08)
    const std::string tf_degree = "alpha=0.8"; 
    std::vector<double> get_weights(){
        std::vector<double> weights(n_modes,0.); 
        weights[0] = 0.4832835680252494;
        weights[1] = 0.18578868783192137;
        weights[2] = 0.14895465942149894;
        weights[3] = 0.14156589063388803;
        weights[4] = 0.14602090608272802;
        weights[5] = 0.15747463079391716;
        weights[6] = 0.17421087214493855;
        weights[7] = 0.19582427577109576;
        weights[8] = 0.2227027457487346;
        weights[9] = 0.2556728375142658;
        weights[10] = 0.29563175116631474;
        weights[11] = 0.34340435613796194;
        weights[12] = 0.4001189781516623;
        weights[13] = 0.4681259669239108;
        weights[14] = 0.5522079782608846;
        weights[15] = 0.6619345490846826;
        weights[16] = 0.8193973619178275;
        weights[17] = 1.0853288274668809;
        weights[18] = 1.6678009937626073;
        weights[19] = 3.7493485726060065;
        return weights; 
    }

    std::vector<double> get_lambdas(){
        std::vector<double> lambdas(n_modes,0.); 
        lambdas[0] = 0.004395756101756371;
        lambdas[1] = 0.09296602814480928;
        lambdas[2] = 0.3201476593618963;
        lambdas[3] = 0.774225587457129;
        lambdas[4] = 1.6330293101705726;
        lambdas[5] = 3.236146544079618;
        lambdas[6] = 6.226106014285317;
        lambdas[7] = 11.825020563862855;
        lambdas[8] = 22.38932012255858;
        lambdas[9] = 42.543210020010946;
        lambdas[10] = 81.50529916757719;
        lambdas[11] = 157.86013404948963;
        lambdas[12] = 309.4303125978937;
        lambdas[13] = 614.3439837061134;
        lambdas[14] = 1238.4682331736067;
        lambdas[15] = 2550.804010865617;
        lambdas[16] = 5442.9252180879;
        lambdas[17] = 12430.412934774975;
        lambdas[18] = 33046.111031989705;
        lambdas[19] = 135160.30793447944;
        return lambdas; 
    }

    double get_w_infinity(){
        return 6.06798473e-06; 
    }
#endif 

#if defined(alpha_09)
    const std::string tf_degree = "alpha=0.9"; 
    std::vector<double> get_weights(){
        std::vector<double> weights(n_modes,0.); 
        weights[0] = 0.6881765710067753;
        weights[1] = 0.1257323620352659;
        weights[2] = 0.08613219191302444;
        weights[3] = 0.0735296486046057;
        weights[4] = 0.06958616634651485;
        weights[5] = 0.06972738345522947;
        weights[6] = 0.07223790821026307;
        weights[7] = 0.07632048363971183;
        weights[8] = 0.0815669280590884;
        weights[9] = 0.08784676189524011;
        weights[10] = 0.09529177329783373;
        weights[11] = 0.10419480944556943;
        weights[12] = 0.11486450362928662;
        weights[13] = 0.12755810011091004;
        weights[14] = 0.1426400833631656;
        weights[15] = 0.16113104627813848;
        weights[16] = 0.1859138191257834;
        weights[17] = 0.22535144589736877;
        weights[18] = 0.30794411829153406;
        weights[19] = 0.5815781406410866;
        return weights; 
    }

    std::vector<double> get_lambdas(){
        std::vector<double> lambdas(n_modes,0.); 
        lambdas[0] = 0.0018120445213826767;
        lambdas[1] = 0.07274585929183396;
        lambdas[2] = 0.2569090821465522;
        lambdas[3] = 0.6178129731725768;
        lambdas[4] = 1.283414940451642;
        lambdas[5] = 2.4970894234422136;
        lambdas[6] = 4.718643364784542;
        lambdas[7] = 8.821606078850383;
        lambdas[8] = 16.47925808795416;
        lambdas[9] = 30.93502165688226;
        lambdas[10] = 58.601445109478604;
        lambdas[11] = 112.5091685953163;
        lambdas[12] = 219.95050462169357;
        lambdas[13] = 439.8134625738113;
        lambdas[14] = 903.144630983625;
        lambdas[15] = 1913.9470384234892;
        lambdas[16] = 4228.905007757849;
        lambdas[17] = 10000.139079698192;
        lambdas[18] = 27252.870070461497;
        lambdas[19] = 110930.04995714132;
        return lambdas; 
    }

    double get_w_infinity(){
        return 9.06187911e-07; 
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

double get_theta(double dt){
    double theta = 0; 
    std::vector<double> weights = get_weights(); 
    std::for_each(weights.begin(), weights.end(), [&] (double d) {theta += d;}); 
    return theta * dt + get_w_infinity(); 
}

double get_delta(double dt){
    double delta = 0; 
    std::vector<double> weights = get_weights(); 
    std::for_each(weights.begin(), weights.end(), [&] (double d) {delta += d;}); 
    return delta; 
}