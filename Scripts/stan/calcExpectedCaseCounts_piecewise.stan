
    // expanded lambdas: every time point has an associated lambda
    real lambda_times_expanded[ len_cohorts ];
    // time-dependent Qe;
    matrix[4,len_times] Qe_t;
    // expected counts for each age
    matrix[ len_ages , len_times ] ageCountExpected;
    // expected counts for each age group
    matrix[ len_Ag , len_times ] countExpected;

    lambda_times_expanded = expandParam(lambda, len_preTimes);
    // if fitting Q(i,a)
    if(len_binnedQeA){
        // set Qe_t to matrix of ones
        matrix[0,4] emptymat;
        Qe_t = getQe_t(emptymat, len_times, binwidth_QeT);
    } else {
        Qe_t = getQe_t(binnedQe, len_times, binwidth_QeT);
    }

    // for each cohort
    for(co in 1:len_cohorts){
        real lambda_cohort[ len_cohorts ];
        matrix[ len_cohorts , 5 ] propInfected;
        row_vector[ len_times ] caseExpected;
        matrix[ 4, len_times ] Qe_a;
        if(len_binnedQeA){
            // 4 x len_times matrix specific to ages of this cohort
            Qe_a = getQe_a(binnedQe, len_times, len_preTimes, len_ages, co, binwidth_QeA);
        }
        
        // calculate probability of being infected by serotype i at time T
        for(i in 1:len_cohorts){
            if((i < co) || ((i - co + 1) > len_ages)){
                lambda_cohort[i] = 0.0;
            } else {
                int age = i - co;
                lambda_cohort[i] = lambda_times_expanded[i];
                // multiplicative age modifier
                lambda_cohort[i] *= K[age+1];
            }
        }
        lambda_cohort = cumulative_sum(lambda_cohort);
        for(i in 1:len_cohorts){
            // store the calculated probability
            lambda_cohort[i] = exp(-lambda_cohort[i]);
        }
        
        // calculate proportion with i-th infection
        for(i in 0:4){
            for(T in 1:len_cohorts){
                propInfected[T,i+1] = (lambda_cohort[T]^(4-i)) * 
                    (1-lambda_cohort[T])^(i) *
                    choose(4,i);
            }
        }

        // calculate proportion with > i-th infection occurred
        for(T in 1:size(lambda_cohort)){
            propInfected[T, ] = 1 - cumulative_sum(propInfected[T, ]) ;
        }
        // difference to calculate incidence
        for(i in 0:3){
            propInfected[ ,i+1] = difference(propInfected[ ,i+1]) ;
        }
        
        // calculate expected case counts
        if(len_binnedQeA){
            caseExpected = [1,1,1,1] * (Qe_a .* propInfected[ (len_preTimes+1): , 1:4 ]') ;
        } else {
            caseExpected = [1,1,1,1] * (Qe_t .* propInfected[ (len_preTimes+1): , 1:4 ]') ;        
        }
    
        for(T in 1:len_times){
            // skip ages that are not in the dataset
            if( ((T + len_preTimes) - co +1) >= len_ages){ continue; }
            if( ((T + len_preTimes) - co +1) <= 0){ continue; }
            ageCountExpected[ (T + len_preTimes) - co +1, T ] = caseExpected[T] * countPop[co,T] ;
        }
    }

    // fill in counts of lower triangle with zero
    for(a in 1:len_ages){
        for(T in 1:len_times){
            if(is_nan(ageCountExpected[a,T])){
                ageCountExpected[a,T] = 0;
            }
        }
    }
    
    if( dim_binnedPhiAT[1] ){
        countExpected = mapAge * ((diag_matrix(PhiA) * ageCountExpected * diag_matrix(PhiT)) .* PhiAT);
    } else {
        countExpected = mapAge *  (diag_matrix(PhiA) * ageCountExpected * diag_matrix(PhiT));
    }
