

- `DHF_add`: Provincial age-stratified dengue hemorrhagic fever (DHF) case counts, age strata as rows and years as columns. Counts in annual surveillance reports published by Bureau of Epidemiology, Ministry of Public Health, Thailand were extracted from hard copies of the reports (1981 to 2005) and from their web portal at http://www.boe.moph.go.th/surdata (2003 to 2017). Counts from provinces which splitted off over the time-series were added back to the original provinces to maintain consistency. 


- `pop`: Population size, age as rows and years as columns. Census data between 1993 to 2017 were extracted from data portal of Department of Provincial Administration, Thailand (http://stat.dopa.go.th/stat/statnew/upstat_age.php). Census in 1980, 1990, and 2000 from the same source was available from a former study. Data in non-census years were imputed as described in the Supplementary Materials of our manuscript.


- `config.RDS`: R object which stores (1) mappings between age and age groups/strata, (2) ages included in the dataset, (3) years included in the dataset, (4) birth cohorts involved, and (5) value in the dataset which denotes a missing value.


