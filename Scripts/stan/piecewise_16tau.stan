functions {
    // function to set lambdas of preTimes to the same value
    // and binned lambdas during Times
    real[] expandParam( real[] params , int len_preTimes){
        int ntimes = 37 + len_preTimes;
        real outParams[ntimes];
        // hard code first 2*8 years of Times being piecewise constant for 2 years
        // followed by being piecewise constant for 3 years after that
        int pointer = 1;
        for (n in 1:ntimes){
            if(n >= (len_preTimes - 1 + 2*pointer + (pointer>9 ? 1 : 0)*(pointer-9))){
                pointer += 1;
            }
            outParams[n] = params[pointer];
        }
        return outParams;
    }
    
    // function to difference, with first element differencing from zero
    vector difference( vector x ){
        vector[num_elements(x)] out;
        out[1] = 0.0;
        out[2: ] = x[1:(num_elements(x)-1)];
        out = x - out ;
        return out;
    }
 
    // function to get Qe_t: expand binnedQe matrix
    matrix getQe_t(matrix binned, int len_expanded, int binwidth){
        matrix[cols(binned), len_expanded] expanded;
        for(i in 1:len_expanded){
            if( i <= ( rows(binned) * binwidth ) ){
                expanded[ , i] = binned[ ( i / binwidth ) + ( ((i % binwidth) > 0) ? 1 : 0 ) , ]' ;
            } else if( rows(binned) ) {
                expanded[ , i] = binned[ rows(binned) , ]' ;
            } else {
                expanded[ , i] = [1.0, 1.0, 1.0, 1.0]';
            }
        }
        return expanded;
    }

    // function to get Qe_a: expand binnedQe matrix
    matrix getQe_a(matrix binned, int len_T, int len_preT, int len_A, int cohort, int binwidth){
        matrix[cols(binned), len_T] expanded;
        for(T in 1:len_T){
            int age = T + len_preT - cohort;                        
            if((age >= 0) && (age <= len_A)){
                expanded[ ,T] = binned[ min(( (age+1) / binwidth )+ ( (((age+1) % binwidth) > 0) ? 1 : 0 ), rows(binned)) , ]' ;
            } else {
                expanded[ ,T] = [0.0,0.0,0.0,0.0]';
            }
        }
        return expanded;
    }

// 
//     // function to expand bins of K into individual ages
//     real[] expandK( real[] binnedK, int len_ages, int binwidth){
//         real K[len_ages];
//         for(a in 1:len_ages){
//             if( (a <= ((size(binnedK)-1)*binwidth + 1) && size(binnedK)  )){
//                 if( a==1 ){
//                     K[a] = binnedK[1];
//                 } else {
//                     // positive integer division yields floored results
//                     K[a] = binnedK[ ((a-2) / binwidth) + 2 ];
//                 }
//             } else {
//                 K[a] = 1;
//             }
//         }
//         return K;
//     }

    // function to expand bins: to be used with PhiA, PhiT
    vector expandPhi( real[] binned, int len_expanded, int binwidth){
        vector[len_expanded] expanded;
        for(i in 1:len_expanded){
            if( i <= ( size(binned) * binwidth ) ){
                expanded[i] = binned[ ( i / binwidth ) + ( ((i % binwidth) > 0) ? 1 : 0 ) ] ;
            } else if( size(binned) ) {
                expanded[i] = binned[ size(binned) ] ;
            } else {
                expanded[i] = 1.0;
            }
        }
        return expanded;
    }

    // function to expand bins of PhiAT
    // : this assumes that PhiAT is used, i.e. use of PhiAT must be checked in prior to calling this function
    matrix expandPhiAT( matrix binned, int[] dim_expanded, int[] binwidth){
        matrix[dim_expanded[1], dim_expanded[2]] expanded;
        for(j in 1:dim_expanded[2]){
            for(i in 1:dim_expanded[1]){
                if( cols(binned) ){
                    expanded[i,j] = binned[
                        (i <= ( rows(binned)*binwidth[1] )) ? 
                        (( i / binwidth[1] ) + ( ((i % binwidth[1]) > 0) ? 1 : 0 )) : 
                         rows(binned)
                    ,(j <= ( cols(binned)*binwidth[2] )) ? 
                        (( j / binwidth[2] ) + ( ((j % binwidth[2]) > 0) ? 1 : 0 )) : 
                        cols(binned)
                    ];
                } else {
                    expanded[i,j] = 1.0;
                }
            }    
        }
        return expanded;
    }

}
data {  

    int<lower=0> len_times;
    int<lower=0> len_cohorts;
    int<lower=0> len_ages;
    int<lower=0> len_Ag;

    int<lower=0> len_binnedQeT; // per infection (i.e. total = 4x)
    int<lower=1> binwidth_QeT;
    int<lower=0> len_binnedQeA; // per infection (i.e. total = 4x)
    int<lower=1> binwidth_QeA;
    
    int<lower=0> len_binnedK;
    int<lower=1> binwidth_K;
    int<lower=0> len_binnedPhiT;
    int<lower=1> binwidth_PhiT;
    int<lower=0> len_binnedPhiA;
    int<lower=1> binwidth_PhiA;
    int<lower=0> dim_binnedPhiAT[2];
    int<lower=0> binwidth_PhiAT[2];

    int naVal; // values for NA counts/popsize

    // observed case counts
    int countObserved[ len_Ag , len_times ];
    // cohort population size
    real countPop[ len_cohorts , len_times ];
    // mapping between Ag and age
    matrix[ len_Ag , len_ages ] mapAge;
}
transformed data {
    int len_preTimes = len_cohorts - len_times ;
}

parameters {

    // hyperparameters
	real<lower=0.0> omega;									// Over-dispersion
	real<lower=0.0> theta;									// Over-dispersion

    // model parameters
    real<lower=0.0,upper=1.0> lambda[ 16 ];
    real<lower=0.0, upper=5.0> binnedK[len_binnedK];
    real<lower=0.0, upper=1.0> binnedPhiT[len_binnedPhiT];
    real<lower=0.0, upper=1.0> binnedPhiA[len_binnedPhiA];
    matrix<lower=0.0, upper=1.0>[dim_binnedPhiAT[1],dim_binnedPhiAT[2]] binnedPhiAT;
    matrix<lower=0.0,upper=1.0>[max({len_binnedQeT,len_binnedQeA,1}), 4] binnedQe;
}
transformed parameters {
    vector[len_ages] K = expandPhi(binnedK, len_ages, binwidth_K);
    vector[len_times] PhiT = expandPhi(binnedPhiT, len_times, binwidth_PhiT);
    vector[len_ages] PhiA = expandPhi(binnedPhiA, len_ages, binwidth_PhiA);
    matrix[len_ages, len_times] PhiAT = expandPhiAT( binnedPhiAT, {len_ages, len_times}, binwidth_PhiAT );
}

model {
    // get expected case counts
#include /Scripts/stan/calcExpectedCaseCounts_piecewise.stan

    // compare with observed counts
    for(T in 1:len_times){
        for(ag in 1:len_Ag){
            if(countObserved[ag,T]==naVal){ continue; }
            if(countPop[ag,T]==naVal){ continue; }
            
            // sampling
            countObserved[ag,T] ~ neg_binomial_2(countExpected[ag,T], countExpected[ag,T]/(  omega + (theta*countExpected[ag,T]) - 1  )) ;
        }
    }

	
	// Weakly informative priors
	omega ~ normal(0.0, 10.0);
	theta ~ normal(0.0, 10.0);

    for(i in 1:size(lambda)){
        lambda[i] ~ normal(0.05, 0.5);
    }
    for(i in 1:len_binnedK){
        binnedK[i] ~ normal(1.0, 3.0);
    }
    binnedQe[ ,1] ~ beta(1.0, 9.0);   // skewed towards 0.0, E(x)=0.1
    binnedQe[ ,2] ~ beta(19.0, 1.0);   // skewed towards 1.0, E(x)=0.95
    binnedQe[ ,3] ~ beta(1.0, 19.0);   // skewed towards 0.0, E(x)=0.05
    binnedQe[ ,4] ~ beta(1.0, 19.0);
    
    if(len_binnedPhiT){
        binnedPhiT ~ normal(0.5, 1.0);
    }
    if(len_binnedPhiA){
        binnedPhiA ~ normal(0.5, 1.0);
    }
    for(j in 1:dim_binnedPhiAT[2]){
        binnedPhiAT[ ,j] ~ normal(0.5, 1.0);    
    }

    // print(target());
}
generated quantities {
    // log-likelihood for the whole dataset
    matrix[ len_Ag, len_times] log_lik;
    real log_lik_sum = 0.0;
    // root mean squared error
    real rmse = 0.0;
    int nEntries = 0;
    // get expected case counts
#include /Scripts/stan/calcExpectedCaseCounts_piecewise.stan

    for(T in 1:len_times){
        for(ag in 1:len_Ag){
            if(
            countObserved[ag,T]==naVal ||
            countPop[ag,T]==naVal
            ){
                log_lik[ag, T] = 1.1;
                continue;
            }
            
            // store log-likelihood
            log_lik[ag, T] = neg_binomial_2_lpmf(
            countObserved[ag, T] | 
            countExpected[ag, T]
            , countExpected[ag, T]/(omega + (theta*countExpected[ag, T]) - 1)
            );
            log_lik_sum += log_lik[ag, T]; 
            
            // store squared error
            rmse += (countExpected[ag, T] - countObserved[ag, T])^2;
            nEntries += 1;
        }
    }
    
    rmse = (rmse/nEntries) ^ .5;

}


