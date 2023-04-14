using namespace std; 
using namespace mfem;

std::vector<double> get_weights_00(){
    std::vector<double> weights(n_modes,0.); 
    return weights; 
}

std::vector<double> get_lambdas_00(){
    std::vector<double> lambdas(n_modes,0.); 
    return lambdas; 
}

double get_w_infinity_00(){
    return 1.; 
}

std::vector<double> get_weights_01(){
    std::vector<double> weights(n_modes,0.); 
    weights[0] = 0.010793652947501125;
    weights[1] = 0.025369224813777753;
    weights[2] = 0.046804267803725805;
    weights[3] = 0.08253719896871657;
    weights[4] = 0.14443941915082775;
    weights[5] = 0.2535010270377004;
    weights[6] = 0.44807563315618576;
    weights[7] = 0.799358111759365;
    weights[8] = 1.4407031430213955;
    weights[9] = 2.6236674300416545;
    weights[10] = 4.826318020788839;
    weights[11] = 8.961469712920204;
    weights[12] = 16.782502301218972;
    weights[13] = 31.74405516518808;
    weights[14] = 61.10762511830639;
    weights[15] = 121.93241875091799;
    weights[16] = 261.63773279711194;
    weights[17] = 651.9297793586896;
    weights[18] = 2273.347327014091;
    weights[19] = 20908.113980193666;
    return weights; 
}

std::vector<double> get_lambdas_01(){
    std::vector<double> lambdas(n_modes,0.); 
    lambdas[0] = 0.0274851516309919;
    lambdas[1] = 0.16907955617276357;
    lambdas[2] = 0.488163055307712;
    lambdas[3] = 1.113544835398026;
    lambdas[4] = 2.2970823799945252;
    lambdas[5] = 4.519827790942369;
    lambdas[6] = 8.703841183009555;
    lambdas[7] = 16.63567442368392;
    lambdas[8] = 31.817871790472765;
    lambdas[9] = 61.191090293142864;
    lambdas[10] = 118.64558769876301;
    lambdas[11] = 232.22623043472143;
    lambdas[12] = 458.9766323686381;
    lambdas[13] = 916.203237939184;
    lambdas[14] = 1851.75818341678;
    lambdas[15] = 3821.2567929340444;
    lambdas[16] = 8214.508582878203;
    lambdas[17] = 19252.955762309404;
    lambdas[18] = 55363.81134097007;
    lambdas[19] = 302326.3979879933;
    return lambdas; 
}

double get_w_infinity_01(){
    return 0.23930518; 
}

std::vector<double> get_weights_02(){
    std::vector<double> weights(n_modes,0.); 
    weights[0] = 0.02814396541624207;
    weights[1] = 0.056454661913285464;
    weights[2] = 0.094335554602779;
    weights[3] = 0.15423639268762612;
    weights[4] = 0.2528668419646304;
    weights[5] = 0.41734521344799386;
    weights[6] = 0.6940201725891882;
    weights[7] = 1.1640388464577658;
    weights[8] = 1.9710918889825286;
    weights[9] = 3.369399492884432;
    weights[10] = 5.805485082163826;
    weights[11] = 10.066350350184866;
    weights[12] = 17.576400672029596;
    weights[13] = 31.01831819938072;
    weights[14] = 55.70226342064903;
    weights[15] = 103.19605128043048;
    weights[16] = 203.59716661938432;
    weights[17] = 459.5021221845185;
    weights[18] = 1405.946956232748;
    weights[19] = 10114.832197076286;
    return weights; 
}

std::vector<double> get_lambdas_02(){
    std::vector<double> lambdas(n_modes,0.); 
    lambdas[0] = 0.023774763069275118;
    lambdas[1] = 0.15931634045122414;
    lambdas[2] = 0.47026996674713545;
    lambdas[3] = 1.0856300503462601;
    lambdas[4] = 2.260924288650546;
    lambdas[5] = 4.487138646427246;
    lambdas[6] = 8.706783217102506;
    lambdas[7] = 16.746659535267778;
    lambdas[8] = 32.19428199899476;
    lambdas[9] = 62.17740701929888;
    lambdas[10] = 120.95540033921253;
    lambdas[11] = 237.16405378990157;
    lambdas[12] = 468.76754100429434;
    lambdas[13] = 935.233855474953;
    lambdas[14] = 1890.2313940312956;
    lambdas[15] = 3900.1511259148338;
    lambdas[16] = 8360.41356336485;
    lambdas[17] = 19429.618557788366;
    lambdas[18] = 54816.91322995089;
    lambdas[19] = 282015.16420763085;
    return lambdas; 
}

double get_w_infinity_02(){
    return 0.05645931; 
}

std::vector<double> get_weights_03(){
    std::vector<double> weights(n_modes,0.); 
    weights[0] = 0.050861591150482645;
    weights[1] = 0.08488770660469716;
    weights[2] = 0.12605654091486035;
    weights[3] = 0.18766378961896524;
    weights[4] = 0.28400570561928223;
    weights[5] = 0.43636168086692445;
    weights[6] = 0.6784780341109015;
    weights[7] = 1.0653239819573321;
    weights[8] = 1.6885510406069548;
    weights[9] = 2.703066936072541;
    weights[10] = 4.369998959661785;
    weights[11] = 7.125738582714147;
    weights[12] = 11.704460488727534;
    weights[13] = 19.38743120745257;
    weights[14] = 32.568294797271484;
    weights[15] = 56.235993405090234;
    weights[16] = 102.80087884639626;
    weights[17] = 212.16998704488302;
    weights[18] = 574.7953932201431;
    weights[19] = 3298.1104941418603;
    return weights; 
}

std::vector<double> get_lambdas_03(){
    std::vector<double> lambdas(n_modes,0.); 
    lambdas[0] = 0.01829305869336155;
    lambdas[1] = 0.13468165329106951;
    lambdas[2] = 0.4028677154585384;
    lambdas[3] = 0.9286786500822911;
    lambdas[4] = 1.923391395581524;
    lambdas[5] = 3.7956723930503866;
    lambdas[6] = 7.333129876910291;
    lambdas[7] = 14.061870114849679;
    lambdas[8] = 26.971419405682894;
    lambdas[9] = 52.00424779243572;
    lambdas[10] = 101.13980194797983;
    lambdas[11] = 198.76584487468463;
    lambdas[12] = 394.882939538117;
    lambdas[13] = 793.1690580177715;
    lambdas[14] = 1613.9785390466782;
    lambdas[15] = 3349.453362294987;
    lambdas[16] = 7210.310709329766;
    lambdas[17] = 16775.047399904262;
    lambdas[18] = 46954.68233526852;
    lambdas[19] = 230731.0818570844;
    return lambdas; 
}

double get_w_infinity_03(){
    return 0.01360695; 
}

std::vector<double> get_weights_04(){
    std::vector<double> weights(n_modes,0.); 
    weights[0] = 0.09535246601176381;
    weights[1] = 0.1325255478900162;
    weights[2] = 0.1799812789429687;
    weights[3] = 0.2523212835071904;
    weights[4] = 0.3633505910517081;
    weights[5] = 0.532025499165395;
    weights[6] = 0.7858012710116449;
    weights[7] = 1.1652288320264905;
    weights[8] = 1.7304350467211287;
    weights[9] = 2.571212637902971;
    weights[10] = 3.8241995881975246;
    weights[11] = 5.701743376644512;
    weights[12] = 8.543399560393517;
    weights[13] = 12.915166674653031;
    weights[14] = 19.82704338261684;
    weights[15] = 31.31269921942847;
    weights[16] = 52.2897793784551;
    weights[17] = 97.9131274449514;
    weights[18] = 235.13729481294942;
    weights[19] = 1095.571458477644;
    return weights; 
}

std::vector<double> get_lambdas_04(){
    std::vector<double> lambdas(n_modes,0.); 
    lambdas[0] = 0.017338196496482848;
    lambdas[1] = 0.144981416700032;
    lambdas[2] = 0.45138262597802975;
    lambdas[3] = 1.0739573652236616;
    lambdas[4] = 2.2948936748914477;
    lambdas[5] = 4.672570607354323;
    lambdas[6] = 9.297979928508173;
    lambdas[7] = 18.29142913498248;
    lambdas[8] = 35.763470047662075;
    lambdas[9] = 69.67377165316839;
    lambdas[10] = 135.46167567640413;
    lambdas[11] = 263.252241188873;
    lambdas[12] = 512.4836559935716;
    lambdas[13] = 1002.6413901358599;
    lambdas[14] = 1981.4134623355885;
    lambdas[15] = 3990.832105583683;
    lambdas[16] = 8342.9572165982;
    lambdas[17] = 18859.569638291418;
    lambdas[18] = 51141.37236156393;
    lambdas[19] = 237318.535411309;
    return lambdas; 
}

double get_w_infinity_04(){
    return 0.00300947; 
}

std::vector<double> get_weights_05(){
    std::vector<double> weights(n_modes,0.); 
    weights[0] = 0.14019054895248748;
    weights[1] = 0.15276201114296667;
    weights[2] = 0.178937938673866;
    weights[3] = 0.22107449799272164;
    weights[4] = 0.2834214471005759;
    weights[5] = 0.37304578686883266;
    weights[6] = 0.5012098046581384;
    weights[7] = 0.6852826334131967;
    weights[8] = 0.9512487253579823;
    weights[9] = 1.3369640636385198;
    weights[10] = 1.8970630277390657;
    weights[11] = 2.7118242096655902;
    weights[12] = 3.9035348997678536;
    weights[13] = 5.6669093886760615;
    weights[14] = 8.333244863243799;
    weights[15] = 12.531718572224108;
    weights[16] = 19.694887123844335;
    weights[17] = 34.012225239060676;
    weights[18] = 72.73429057256043;
    weights[19] = 276.96301730092716;
    return weights; 
}

std::vector<double> get_lambdas_05(){
    std::vector<double> lambdas(n_modes,0.); 
    lambdas[0] = 0.011945230064400096;
    lambdas[1] = 0.11403612443275454;
    lambdas[2] = 0.3555867771099227;
    lambdas[3] = 0.8247675021467651;
    lambdas[4] = 1.6935662593237266;
    lambdas[5] = 3.2861428230351035;
    lambdas[6] = 6.217057630958237;
    lambdas[7] = 11.67243074750414;
    lambdas[8] = 21.98631402998542;
    lambdas[9] = 41.826558457825456;
    lambdas[10] = 80.64089662197266;
    lambdas[11] = 157.74802875371245;
    lambdas[12] = 313.13569218168146;
    lambdas[13] = 631.0242491361266;
    lambdas[14] = 1293.617008712904;
    lambdas[15] = 2713.3866886211918;
    lambdas[16] = 5906.005384788986;
    lambdas[17] = 13818.57887248305;
    lambdas[18] = 38218.53073376462;
    lambdas[19] = 173830.92228173843;
    return lambdas; 
}

double get_w_infinity_05(){
    return 0.00074696; 
}

std::vector<double> get_weights_06(){
    std::vector<double> weights(n_modes,0.); 
    weights[0] = 0.2241112206610028;
    weights[1] = 0.18853448470623824;
    weights[2] = 0.19799240240330232;
    weights[3] = 0.22719317108969442;
    weights[4] = 0.27378199196447933;
    weights[5] = 0.3393638549795428;
    weights[6] = 0.42801144997574725;
    weights[7] = 0.5463401839631146;
    weights[8] = 0.7042358780235634;
    weights[9] = 0.9161206117703352;
    weights[10] = 1.2023748501379974;
    weights[11] = 1.591032665236674;
    weights[12] = 2.1207117486372566;
    weights[13] = 2.8480923410583237;
    weights[14] = 3.8685933118411335;
    weights[15] = 5.369457819779584;
    weights[16] = 7.778384448978061;
    weights[17] = 12.29304961142493;
    weights[18] = 23.510627480292435;
    weights[19] = 74.4585290953474;
    return weights; 
}

std::vector<double> get_lambdas_06(){
    std::vector<double> lambdas(n_modes,0.); 
    lambdas[0] = 0.009970803491035069;
    lambdas[1] = 0.11501260220736824;
    lambdas[2] = 0.3726669719547261;
    lambdas[3] = 0.886070596442776;
    lambdas[4] = 1.8601241615040185;
    lambdas[5] = 3.683334763024911;
    lambdas[6] = 7.0872105176489955;
    lambdas[7] = 13.457505373973765;
    lambdas[8] = 25.449904254920668;
    lambdas[9] = 48.2322399023694;
    lambdas[10] = 92.0316056525584;
    lambdas[11] = 177.40787770680677;
    lambdas[12] = 346.2781492835861;
    lambdas[13] = 685.4322192774556;
    lambdas[14] = 1379.2172758903293;
    lambdas[15] = 2838.9645014236667;
    lambdas[16] = 6069.974621516811;
    lambdas[17] = 13968.878298640953;
    lambdas[18] = 37902.365551856055;
    lambdas[19] = 165401.37159406338;
    return lambdas; 
}

double get_w_infinity_06(){
    return 0.00015986; 
}

std::vector<double> get_weights_07(){
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

std::vector<double> get_lambdas_07(){
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

double get_w_infinity_07(){
    return 3.70179209e-05; 
}

std::vector<double> get_weights_08(){
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

std::vector<double> get_lambdas_08(){
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

double get_w_infinity_08(){
    return 6.06798473e-06; 
}

std::vector<double> get_weights_09(){
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

std::vector<double> get_lambdas_09(){
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

double get_w_infinity_09(){
    return 9.06187911e-07; 
}

std::vector<double> get_weights_10(){
    std::vector<double> weights(n_modes,0.); 
    weights[0] = 1; 
    return weights; 
}

std::vector<double> get_lambdas_10(){
    std::vector<double> lambdas(n_modes,0.); 
    return lambdas; 
}

double get_w_infinity_10(){
    return 0; 
}

std::vector<double> get_weights(double alpha){
    // cout << "Are those equal: " << alpha * 100 << " " << (int)(alpha * 100 + 0.5) << endl; 
    switch((int)(alpha * 100 + 0.5)){
        case 0: 
            return get_weights_00(); 
        case 10: 
            return get_weights_01(); 
        case 20: 
            return get_weights_02(); 
        case 30: 
            return get_weights_03(); 
        case 40: 
            return get_weights_04(); 
        case 50: 
            return get_weights_05(); 
        case 60: 
            return get_weights_06(); 
        case 70: 
            return get_weights_07(); 
        case 80: 
            return get_weights_08(); 
        case 90: 
            return get_weights_09(); 
    }
    return get_weights_10(); 
}

std::vector<double> get_lambdas(double alpha){
    // cout << "Are those equal: " << alpha * 100 << " " << (int)(alpha * 100 + 0.5) << endl; 
    switch((int)(alpha * 100 + 0.5)){
        case 0: 
            return get_lambdas_00(); 
        case 10: 
            return get_lambdas_01(); 
        case 20: 
            return get_lambdas_02(); 
        case 30: 
            return get_lambdas_03(); 
        case 40: 
            return get_lambdas_04(); 
        case 50: 
            return get_lambdas_05(); 
        case 60: 
            return get_lambdas_06(); 
        case 70: 
            return get_lambdas_07(); 
        case 80: 
            return get_lambdas_08(); 
        case 90: 
            return get_lambdas_09(); 
    }    
    return get_lambdas_10();     
}

double get_w_infinity(double alpha){
    // cout << "Are those equal: " << alpha * 100 << " " << (int)(alpha * 100 + 0.5) << endl; 
    switch((int)(alpha * 100 + 0.5)){
        case 0: 
            return get_w_infinity_00(); 
        case 10: 
            return get_w_infinity_01(); 
        case 20: 
            return get_w_infinity_02(); 
        case 30: 
            return get_w_infinity_03(); 
        case 40: 
            return get_w_infinity_04(); 
        case 50: 
            return get_w_infinity_05(); 
        case 60: 
            return get_w_infinity_06(); 
        case 70: 
            return get_w_infinity_07(); 
        case 80: 
            return get_w_infinity_08(); 
        case 90: 
            return get_w_infinity_09(); 
    }    
    return get_w_infinity_10(); 
}

std::vector<double> get_gammas(double alpha, double dt){
    std::vector<double> gammas(get_lambdas(alpha)); 
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::multiplies<double>(), std::placeholders::_1, dt));
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::plus<double>(), std::placeholders::_1, 1));
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::divides<double>(), 1, std::placeholders::_1));
    return gammas; 
} 

double get_beta(double alpha, double dt){
    double beta = 0; 
    std::vector<double> gammas = get_gammas(alpha, dt); 
    std::vector<double> weights = get_weights(alpha); 
    std::transform(gammas.begin(), gammas.end(), gammas.begin(),std::bind(std::multiplies<double>(), std::placeholders::_1, dt));
    std::transform(gammas.begin(), gammas.end(), weights.begin(), gammas.begin(), std::multiplies<double>()); 
    std::for_each(gammas.begin(), gammas.end(), [&] (double d) {beta += d;}); 
    return beta + get_w_infinity(alpha); 
}

double get_theta(double alpha, double dt){
    double theta = 0; 
    std::vector<double> weights = get_weights(alpha); 
    std::for_each(weights.begin(), weights.end(), [&] (double d) {theta += d;}); 
    return theta * dt + get_w_infinity(alpha); 
}

double get_delta(double alpha, double dt){
    double delta = 0; 
    std::vector<double> weights = get_weights(alpha); 
    std::for_each(weights.begin(), weights.end(), [&] (double d) {delta += d;}); 
    return delta; 
}