### Sampling probability proportional to size
### Scenario: 1,000 schools
### Two sizes, L and S
### 2023.05.15
### 2023.10.16 ## added ICC
### 2023.12.05 
### redefined w2, dependent on size 
### using mix for single and no weights condition 
### added DuMouchel and Duncan test
### 2024.01.10 Testing out minority %s w NR
### 2024.01.28 added classroom sampling
### added raw variance components
### 2024.04.30 added house weight--> same as normalized weight
### Auth: Francis Huang / Umut Atasever
### huangf@missouri.edu

# If I just want to run one condition
# start.c <- end.c <- fh <- 1; reps <- 1
# radical <- FALSE #not used in sim
# class_sampling <- TRUE

simwgt <- function(start.c, end.c, reps = 5, seed = 123){
  ## remove seed option later...
  ## will modify later to save output better
  sta <- Sys.time()
  res <- list()
  ic.c <- c(1:2) #took out 3; had boundary issues 
  sa <- c(100, 150, 200)
  nr.sch <- c("TRUE", "FALSE")
  nr.stu <- c("TRUE", "FALSE")
  radical <- c("FALSE") #just make it default
  class_sampling <- c("TRUE", "FALSE")
  cond <- expand.grid(ic.c = ic.c, nr.sch = nr.sch, radical = radical,
                      class_sampling = class_sampling,
                      nr.stu = nr.stu, rps = reps, sa = sa)
  cond$id <- 1:nrow(cond)
  cond <- cond[start.c:end.c, ]
  
  ic <- cond$ic.c
  nr.sch <- cond$nr.sch
  nr.stu <- cond$nr.stu
  rps <- cond$rps
  sa <- cond$sa
  radical <- cond$radical
  class_sampling <- cond$class_sampling

library(dplyr)
library(WeMix)
library(lme4)
library(survey)


## number of schools
j2 <- 1000 # large schools
j1 <- 9000 # small schools
totj <- j2 + j1 # total schools
# ic <- 2 #change to 1 for ICC = .5, 2 for ICC = .25, 3 for ICC = .05

dfr <- data.frame(id = 1:totj) # this is the school id

# set a random seed for reproducibility for the population
set.seed(1112) #1112

### radical size #2024.01.22
if (radical == TRUE){
  dfr$size[1:j2] <- sample(101:500, size = j2, replace = TRUE)
  dfr$size[(j2 + 1):totj] <- sample(20:100, size = j1, replace = TRUE)
} else {
  dfr$size <- round(rlnorm(totj, 4.5, .4))
  #follows a log normal distribution
}

summary(dfr$size)

dfr$w1 <- rnorm(totj, 0, 1) # l2 continuous normal predictor
dfr$w2 <- rbinom(totj, 1, plogis(scale(dfr$size) - .25))
totsize <- sum(dfr$size) # total number of students, based on size
dfr$private <- rbinom(totj, 1, plogis(dfr$w1 - 2.5)) # make it related to w1
table(dfr$private) # just checking
dfr$e2 <- rnorm(totj, 0, 2) # tau00 is 4

dfr$prob <- dfr$size / sum(dfr$size) # probability of selection
sum(dfr$prob) # should sum to 1

### dfr is a school-only dataset
### expand out the dataset (to students)
### This data dataset is equivalent to total number of students

dat <- data.frame(
  id = rep(dfr$id, dfr$size),
  size = rep(dfr$size, dfr$size),
  private = rep(dfr$private, dfr$size),
  w1 = rep(dfr$w1, dfr$size),
  w2 = rep(dfr$w2, dfr$size),
  e2 = rep(dfr$e2, dfr$size)
)

Ns <- nrow(dat)

## insert classrooms... 2024.01.22

#class_sampling <- TRUE # Say TRUE for class sampling; change later
avg <- 25 # average class size for class sampling

num_main_groups <- round(nrow(dat) / avg)

# Add a column for the main group ID
dat <- dat %>%
  group_by(id) %>%
  mutate(class = rep(1:num_main_groups, each = avg, length.out = n())) %>%
  group_by(id, class) %>%
  mutate(n_std = n()) %>% #how many students in the class
  mutate(class = ifelse(n_std < avg / 2, 1, class)) %>% #if less than avg, then put in class 1
  mutate(n_std = n()) %>%
  ungroup() #count again?

classes <- dat %>% group_by(id, class) %>%
  summarise(n_std = n(), .groups = "drop_last") %>% #don't really need n_std
  mutate(n = n()) # how many classrooms

n_cls <- classes %>% tally() # number of classrooms
ncls <- sum(n_cls$n) #number of classrooms
dat <- left_join(dat, n_cls, by = "id") #merge it with main dataset

dat$cl_re <- rep(rnorm(ncls, 0, 2), classes$n_std) #random effects for classroom

### END CLASSROOM addition

sigm <- switch(ic, 2.81, 7.03, 19.7) #sigma2 for level 1
## this controls the ICC; modifying the sigm (in SD vs sigma2)
dat$x1 <- rnorm(Ns) # independent x1
dat$minority <- rbinom(Ns, 1, .1)
dat$e1 <- rnorm(Ns, 0, sigm) 

# combine to make y
dat$y <- 50 + dat$w1 * 3 + dat$private * 3 +
  scale(dat$size) * 1 +  dat$cl_re + #
  dat$w2 * 3 +
  dat$x1 * 3 + dat$minority * -3 +
  dat$e2 + dat$e1
dat$one <- 1 # constant, just in case needed later
dat$clroom <- paste0(dat$id, "-", dat$class)

### THESE are the population values: no need for weights
m0 <- lmer(y ~ (1 | id), data = dat)
#summary(m0)
performance::icc(m0)
tau00 <- data.frame(VarCorr(m0))[1, 4] #true tau00
sigma2 <- data.frame(VarCorr(m0))[2, 4] #true sigma2

m1 <- lmer(y ~ w1 + w2 + x1 + minority +
             (1 | id), data = dat)
summary(m1, cor = FALSE)
tbeta <- fixef(m1) # true beta
# tbeta has the population values

### needed for sampling

# Step 1: Calculate the minority proportion for each school
school_minority_proportions <- dat %>%
  group_by(id) %>%
  summarize(minority_proportion = mean(minority))

# Step 2: Create a new variable indicating if the school has high minority proportions
threshold <- 0.10  # You can adjust this threshold as needed
school_minority_proportions$high_minority <- ifelse(school_minority_proportions$minority_proportion >= threshold, 1, 0)
dfr <- left_join(dfr, school_minority_proportions, by = "id")

table(dfr$high_minority)

# DuMouchel and Duncan (1983) test
ddtest <- function(){
  h0 <- lmer(y ~ w1 + w2 + x1 + minority +
         (1 | id), data = sampled_data, REML = FALSE)
  h1 <- lmer(y ~ (w1 + w2 + x1 + minority) * nwt +
               (1 | id), data = sampled_data, REML = FALSE)
  pv <- data.frame(anova(h0, h1))[2, 8] #pvalue is h1 is better than h0
}
 
pstest <- function(){
  h0 <- lm(y ~ w1 + w2 + x1 + minority, data = sampled_data)
  h1 <- lm(nwt ~ resid(h0), data = sampled_data)
  pv <- summary(h1)$coef[2,4] #pvalue is h1 is better than h0
} 

# Now, dat is our population dataset 
# We want to draw
# n_sample schools (based on PPS). We can oversample private schools
# which is now 10% of the schools (private = 1)
#
# dfr is a school-level only dataset too.

###########################
# School Sample Selection #
###########################

## School Sample Selection (program can do a PPS sample with public & private schools as explicit stratification and also without stratification. It can generate random non-response on school and std levels)

random <- NULL # Initialize an empty variable 'random' to store random numbers
stratification <- TRUE # Say TRUE if you want to use stratification (public / private)
treat_small_schools <- FALSE #  Say TRUE if you want to set the MOS of small schools same 

n_sample <- sa # Set the number of school samples to be taken
reps <- rps # Set the number of iterations to perform
non_response <- nr.stu # Say TRUE if you want non-response at student level: this will make students with non-minority background to tend to refuse more
non_response_school <- nr.sch # Say TRUE for school non response: depends on school size threshold and private schools
class_sampling <- class_sampling
# Allocation check

frame <- dfr

# summary_data <- frame %>% group_by_at(c("private")) %>% 
#   reframe(N = n(), MOS = sum(size))

summary_data <- frame %>% group_by(private) %>%
  summarise(N = n(), MOS = sum(size))

# Calculate a number proportional to total N (e.g., 200) for each category
# We will use the "N" value in each group to calculate the proportion and then multiply by n_sample

summary_data %>%
  mutate(proportional_N = (MOS / sum(MOS)) * n_sample)

####################################
# weights generated in the program #
####################################

# schwgt: base weight for the reverse selection probability of the sampled schools 
# sch_adj: adjustment factor for school non-response
# schwgt: after the non response adjustment is schwgt = schwgt * sch_adj
# stdwgt: base weight for the reverse selection probability of the sampled students
# std_adj: adjustment factor for student non-response
# stdwgt: after the non response adjustment is stdwgt = stdwgt * std_adj
# totwgt: total weight, which is a multiplication of schwgt and stdwgt (including adjustment factors)
# totwgt = schwgt * stdwgt

###############################################
# Stratification and Allocation of the Sample #
###############################################

# Stratification vars (public vs private)
stratification_vars <- c(0, 1)

#For 100:  80 to 20 (oversampling private with 20 instead of 10)
#For 150: 120 to 30 (oversampling private with 30 instead of 15)
#For 200: 160 to 40 (oversampling private with 40 instead of 20)

## how many schools to sample:: oversample private
## based on n_sample:: 80:20, true value is 90:10; we are oversampling
n_sample_str <- c(n_sample * .8, n_sample * .2)

# n_sample_str <- c(120, 30)

############################################################################
# Non-Response Model (see the line 336 for school and and 658 for students #
############################################################################


# Define thresholds to discretion mos into categories
# they are used for non-response simulation in school and student levels

thresholds <- c(100, 300)  # Adjust thresholds as needed; moved up to 300 from 200

# Small, Medium, Large and Other
school_non_response_prob <- c(0.15, 0.10, 0.05)

# Small, Medium, Large and Other
private_school_fac <- c(1.15, 1.25) 

# Minority=1 and Minority=0
# Minority Non Response Prob on Student level
minority_non_response_prob <- c(0.2, 0.05) #changed

# Create an empty list file to record participation rates

overall_participation_rate <- list()

# Set the range of students to sample per school
min_students_per_school <- 24
max_students_per_school <- 26
#### start of simulation

### initialize containers, cfs = coefficients, ses = standard errors
### for weight type of 1 to 4 (see below later for explanation)
### four types of weights are used
### 1: weights at both levels
### 2: weight at level one only
### 3: no weight
### 4: weight only at level 2, level 1 weight is 1.0


cfs1 <- cfs2 <- cfs3 <- cfs4 <- cfs5 <- 
  ses1 <- ses2 <- ses3 <- ses4 <- ses5 <- matrix(NA, nrow = reps, ncol = 5)
icc1 <- icc2 <- icc3 <- icc4 <- ht <- ht2 <- ps <- min1 <- min2 <- numeric(reps) # a vector to track ICC
tmp1 <- tmp2 <- tmp3 <- tmp4 <- list()
vcomp1 <- vcomp2 <- vcomp3 <- vcomp4 <- matrix(NA, nrow = reps, ncol = 2)
vcomp1se <- vcomp2se <- vcomp3se <- vcomp4se <- matrix(NA, nrow = reps, ncol = 2)

summary_stats <- NULL
if(exists("seed")) set.seed(seed = seed) ### check

## loop for simulation

for (fh in 1:reps) {
  
  # Initialize an empty variable to store sampled schools
  sampled_schools_str <- NULL
    
    # Loop through each stratification variable
    for (strv in 1:length(stratification_vars)){
      
    # don't name as a function (str),  changed to strv
      
      # Copy the original dataset to 'schools'
      schools <- frame
      
      # Check if stratification is enabled
      if (stratification == TRUE) {
        
      # Filter schools based on the current stratification variable
      schools <- schools %>% 
        filter(private == stratification_vars[strv])
      }
      # Get the new sample size for the current stratification variable
      # how many to sample
      n_sample_new <- n_sample_str[strv]
      
      # Copy the sample size to a variable for model estimation
      n_sample_model <- n_sample_new
      
      # Generate a random number between 0 and 1
      random_number <- runif(1)
      
      # Append the random number to the 'random' vector
      random <- rbind(random, random_number)
      
      # Print random numbers for PPS sampling
      # cat("printing random numbers to generate PPS samples:", random_number, '\n')

      
      ### IGNORABLE FOR NOW...
      
      if (treat_small_schools == TRUE){
      
          # treat the small schools: I turned off this. 
          # Schools that are smaller than 84 are affected...
    
          # Calculate the mean of schools$size
          mean_size <- mean(schools$size)
    
          # Update the size of schools below the mean
          schools$size[schools$size < mean_size] <- mean(schools$size[schools$size < mean_size])

      }
      
      ### CONTINUE FROM HERE...
      
      # sort schools by size in descending order
      schools <- schools[order(schools$size, decreasing = TRUE), ]
      
      # Calculate the sampling interval for PPS sampling
      sampling_interval <- sum(schools$size) / n_sample_model
      
      # Mark schools as 'S' (selected) based on PPS sampling
      
      for (chosen in 1:nrow(schools)){
        if((schools$size[chosen] > sampling_interval) == TRUE){
          schools[chosen, "sampling_status"] <- "S"
          certainity_schools <- "Yes"
          print("xx") #just a test, no school has size > sampling interval
        }
      }
      
      # Check if certainty_schools exist
      if(exists("certainity_schools")){
        
        # Filter certainity_schools based on sampling status
        certainity_schools <- schools %>% filter(sampling_status == "S")
        
        # Filter schools based on certainity_schools' ids
        schools <- schools %>% filter("id" %in% certainity_schools$id)
        
        # Update the sample size for the current stratification variable
        n_sample_model <- n_sample_model - nrow(certainity_schools)
        
        # Update the sampling interval based on the remaining schools
        sampling_interval <- sum(schools$size) / n_sample_model
      }
      
      # Calculate the selection number for each school
      selection_number <- random_number * sampling_interval
      
      # Update schools dataframe with additional columns for PPS sampling
      schools <- schools %>% 
        ungroup() %>% 
        mutate(
          mos = size,
          cumsum = cumsum(mos),
          sampling.interval = sum(mos) / n_sample_model
      )
      
      # Initialize an empty vector for selection numbers
      selection_numbers <- NULL

      sel <- rep(sampling_interval, n_sample_model)
      sel[1] <- selection_number
      selection_numbers2 <- cumsum(sel)
      
      # Initialize a variable for rows
      rows <- NA
      
      # Initialize 'sampling_status' column in schools dataframe
      schools$sampling_status <- NA
      
      # Mark selected schools as 'S' based on PPS sampling
      for (i in 1:n_sample_model){
        diff <- round(schools$cumsum - selection_numbers2[i]) #changed to my sel2
        selected_row <- which(diff == min(diff[diff >= 0]))
        if(is.na(schools[selected_row, "sampling_status"])){
          schools[selected_row, "sampling_status"] <- "S" 
        } else {
          selected_row <- which(diff == diff[diff >= 0][2]) #if chosen, get next one
          schools[selected_row, "sampling_status"] <- "S"
        }
      }
      
      # Calculate total_N and other probability-related columns
      schools$total_N <- sum(schools$mos)
      sampled_schools <- schools %>% 
        mutate(
          prob = n_sample_model * mos / total_N,
          summed_prob = sum(prob)
          ) %>% 
          filter(sampling_status == "S") %>% 
          mutate(schwgt = 1 / prob)
      
      
      # Check if certainty_schools exist and bind them to sampled_schools
      if(exists("certainity_schools")){
        certainity_schools <- certainity_schools %>% 
          mutate(
            mos = size,
            prob = 1,
            schwgt = 1,
            total_estimation = NA
        )
        
        sampled_schools <- bind_rows(sampled_schools, certainity_schools)
      }
      
####################################
## Non-response model for schools ##
####################################
      
      # Check if non_response_school is TRUE
      if(non_response_school == TRUE){
        
        # Create a new column 'school_size_category' based on mos
        sampled_schools <- sampled_schools %>%
          mutate(
            school_size_category = case_when(
              mos <= thresholds[1] ~ "Small",
              mos <= thresholds[2] ~ "Medium", #end up having little
              TRUE ~ "Large"
            )
          )
        
        # Simulate no-response probabilities based on minority and school size
        
        sampled_schools <- sampled_schools %>%
          mutate(
            school_size_prob = case_when(
              school_size_category == "Small" ~ school_non_response_prob[1],
              school_size_category == "Medium" ~ school_non_response_prob[2],
              school_size_category == "Large" ~ school_non_response_prob[3],
              TRUE ~ school_non_response_prob[1]
            ),    
            private_prob = ifelse(private == 1, private_school_fac[2], private_school_fac[1]),
            no_response_prob = school_size_prob * private_prob,
            no_response = rbinom(n(), 1, no_response_prob),  # Simulate no-response based on probability
            num_schools = n()
          ) 
        
        # Non response adjustment
        sampled_schools <- sampled_schools %>% 
          mutate(
          no_response_sum = sum(no_response == 1),
          sch_adj = (num_schools) / (num_schools - no_response_sum),
          schwgt = schwgt * sch_adj
        ) %>% 
          filter(no_response == 0) #only schools that responded
        
      }
      
      # Bind the current sampled_schools to the overall sampled_schools_str
      sampled_schools_str <- rbind(sampled_schools_str, sampled_schools)
      
      # Check if certainty_schools exist and remove them
      if(exists("certainity_schools")){
        rm(certainity_schools)
      }
    }
    
  # Combine all sampled_schools_str into a single dataframe
  sampled_schools <- sampled_schools_str
    
  # Calculate the total estimation based on schwgt and mos
  sampled_schools <- sampled_schools %>% 
    mutate(total_estimation = sum(schwgt * mos))
  
  # estimated total population size
  sum(unique(sampled_schools$total_estimation))
  
  # simulate non-response randomly for schools then calculate the adjustment factor
  

  ##########################
  ### Student Sampling #####
  ##########################
  
  # Select relevant columns from sampled_schools dataframe
  sampled_schools <- sampled_schools %>% 
    select(c("id", "high_minority", "mos", "cumsum", "schwgt"))
  
  # Left join sampled_schools with the original data 'dat'
  sampled_data <- left_join(sampled_schools, dat, by = c("id"), multiple = 'all')
  
  ##############################
  ### classroom sampling #######
  ##############################
  
  if(class_sampling == TRUE){
  class2 <- filter(classes, id %in% sampled_data$id)  
  sel <- class2 %>% group_by(id) %>%
    slice_sample(n = 1) %>% #randomly selects
    mutate(clswgt = 1/(1/n)) %>% # or just n
    select(-n_std, -n)
  
  sampled_data <- left_join(sel,
    sampled_data, by = c("id", "class"), 
    multiple = 'all') %>%
    mutate(stdwgt = clswgt,
           totwgt = stdwgt * schwgt)
    
    
  } else { ### the following is sampling the students per school
  
  # Generate a random number of students to sample for each school
  sampled_data <- sampled_data %>% 
    group_by(id) %>% 
    mutate(
      students_to_sample = sample(min_students_per_school:max_students_per_school, 1)
    ) %>% 
    mutate(
      students_to_sample = ifelse(students_to_sample > mos, mos, students_to_sample)
    ) %>% 
    sample_n(size = unique(students_to_sample), replace = TRUE) %>% 
    mutate(num_students = n())
  
  # Calculate the selection probability and student weight
  sampled_data <- sampled_data %>% mutate(
    selection_probability = num_students / mos,
    stdwgt = 1 / selection_probability,
    totwgt = schwgt * stdwgt  # calculate the total weight
  )
}
  
  
  # Count the number of participating students
  students_sampled <- nrow(sampled_data)
  
#####################################
## Non-response model for students ##
#####################################
  
  # Check if non-response modeling is enabled for students
  if (non_response == TRUE) {
    
    # Simulate no-response probabilities based on minority
    sampled_data <- sampled_data %>%
      mutate(
        minority_prob = ifelse(minority == 1, minority_non_response_prob[1], minority_non_response_prob[2]),
        no_response_prob = minority_prob,
        no_response = rbinom(n(), 1, no_response_prob)  # Simulate non-response based on probability
      ) %>%
      select(-minority_prob)
    
    if (class_sampling == TRUE){
      
      # Non response adjustment
      sampled_data <- sampled_data %>% 
        group_by(id, class) %>% #has class
        mutate(
          num_students = n(),
          no_response_sum = sum(no_response == 1),
          std_adj = (num_students) / (num_students - no_response_sum),
          stdwgt = stdwgt * std_adj,
          totwgt = schwgt * stdwgt  # re-calculate the total weight
        ) %>% 
        filter(no_response == 0)
      
    } else { #stu.nr = T and cl_s = F
    
    # Non response adjustment
    sampled_data <- sampled_data %>% 
      group_by_at("id") %>% #does not have class
      mutate(
      no_response_sum = sum(no_response == 1),
      std_adj = (num_students) / (num_students - no_response_sum),
      stdwgt = stdwgt * std_adj,
      totwgt = schwgt * stdwgt  # re-calculate the total weight
    ) %>% 
      filter(no_response == 0)
  }
}
  ###################################################
  ## Calculate participation rates and print stats ##
  ###################################################
  
  cat("Condition:", cond$id, "\n")
  # Count the number of participating schools
  number_of_participating_schools <- nrow(sampled_schools)
  
  # Placeholder for school participation rate
  cat("Number of participating schools:", number_of_participating_schools, "\n")

  school_participation_rate <- nrow(sampled_schools) / n_sample
  cat("Participation rates\n School:", school_participation_rate, "\n")
  
  # Print a message for student sampling
 
  # Count the number of participating students
  number_of_participating_students <- nrow(sampled_data)
  
  # Calculate the student participation rate
  student_participation_rate <- (nrow(sampled_data) / students_sampled)
  
  # Print the student participation rate
  cat(" Student:", student_participation_rate, "\n")
  
  # Calculate the overall participation rate for this replication
  overall_participation_rate[fh] <- school_participation_rate * student_participation_rate * 100
  

  # Print the overall participation rate for this replication
  cat(" Overall:", as.numeric(overall_participation_rate[fh]), "\n")
  
  # Print summary of total weight
  cat("Summary of total weight\n")
  print(summary(sampled_data$totwgt))
  
  # Calculate the total sum of weights
  sum_of_total_weight <- sum(sampled_data$totwgt)
  
  # Print the sum of total weights
  
  cat("Sum of total weight:", sum(sampled_data$totwgt), "\n")
  
  population_mean <- mean(dat$y)
  cat("Population mean:", population_mean,'\n')
  
  unweighted_sample_mean <- mean(sampled_data$y)
  cat("Unweighted sample mean:", unweighted_sample_mean,'\n')

  weighted_sample_mean <- weighted.mean(x = sampled_data$y, w = sampled_data$totwgt)
  cat("Weighted sample mean:", weighted_sample_mean,'\n')

  cat(tau00,":", sigma2, '\n')
  
  # Gather statistics into a dataframe
  sample_stats <- cbind(fh,
                        number_of_participating_schools,
                        school_participation_rate,
                        number_of_participating_students,
                        student_participation_rate,
                        overall_participation_rate,
                        sum_of_total_weight,
                        population_mean,
                        unweighted_sample_mean,
                        weighted_sample_mean
  )
  
  # print(sample_stats)
  
  summary_stats <- rbind(summary_stats, sample_stats)
  
  # print(summary_stats)
  
  #normalized weight for lmer with weights at level one only
  sampled_data$nwt <- sampled_data$totwgt / mean(sampled_data$totwgt)
  sampled_data$one <-1 #just a weight of one to be used when only using L2 weights
  err <- 0 #no error, this is an error counter to track non convergence
  
  # Calculate scaling factors
  scaling_factor_senate <- 500 / sum(sampled_data$totwgt)
  scaling_factor_house <- nrow(sampled_data) / sum(sampled_data$totwgt)
  
  # Scale the variables
  sampled_data$senate_wgt <- sampled_data$totwgt * scaling_factor_senate
  sampled_data$house_wgt <- sampled_data$totwgt * scaling_factor_house
  
  # Check the sums
  sum(sampled_data$senate_wgt)
  sum(sampled_data$house_wgt)
  
  # weights at twolevels; tryCatch is used to track nonconvergence errors
  # only mix has issues, lmer does not in our cases
  # a null model is also run to track unconditional ICCs
  
  ### four types of weights are used
  ### 1: weights at both levels
  ### 2: weight at level one only
  ### 3: no weight
  ### 4: weight only at level 2, level 1 weight is 1.0
  
  ### adding ecluster
  ### see https://www.statmodel.com/download/Scaling3.pdf
  wttmp <- sampled_data %>% group_by(id) %>%
    summarise(nj = sum(stdwgt)^2 / sum(stdwgt^2),
              tots = sum(stdwgt),
              eadj = nj / tots)
  
  sampled_data <- left_join(sampled_data, wttmp, by = 'id')
  sampled_data$eclust <- with(sampled_data, stdwgt * eadj)
    
  # weights at two levels
  wgt2levs <-
    tryCatch(
      mix(
        y ~ w1 + w2 + x1 + minority +
          (1 | id),
        weights = c('stdwgt', 'schwgt'),
        data = sampled_data, cWeights = TRUE
      ),
      error = function(x) {
        err <<- 1 #trigger an error counter if no convergence
      }
    )
  
  
  # null model for ICC
  wgt2levs.icc <- mix(y ~ (1 | id),
    weights = c('stdwgt', 'schwgt'),
    data = sampled_data, cWeights = TRUE)
  
  if (err == 0) {
    wgt1a <- mix(
      y ~ w1 +  w2 + x1 + minority +
        (1 | id),
      weights = c('house_wgt', 'one'),
      data = sampled_data,
      cWeights = TRUE
    )
    
    wgt1a.icc <- mix(
      y ~ 1 +
        (1 | id),
      weights = c('house_wgt', 'one'),
      data = sampled_data,
      cWeights = TRUE
    )
    
    wgt2only <- mix(y ~ w1 + w2 + x1 + minority +
                       (1 | id),
                     weights = c('one', 'schwgt'),
                     data = sampled_data, 
                    cWeights = TRUE
    )
    wgt2only.icc <-
      mix(y ~ (1 | id),
          weights = c('one', 'schwgt'),
          data = sampled_data, 
          cWeights = TRUE
      )
    
    nowgt1a <- mix(
      y ~ w1 +  w2 + x1 + minority +
        (1 | id),
      weights = c('one', 'one'),
      data = sampled_data,
      cWeights = TRUE
    )
    
    nowgt1a.icc <- mix(
      y ~ 1 +
        (1 | id),
      weights = c('one', 'one'),
      data = sampled_data,
      cWeights = TRUE
    )
    
    des <- svydesign(ids = ~id, weights = ~house_wgt, data = sampled_data)
    slev <- svyglm(y ~ w1 + w2 + x1 + minority, des)
    
    ## Hausmann DC test
    # ht[fh] <- haustest(wgt2levs, nowgt1a)
    ht[fh] <- ddtest()
    ps[fh] <- pstest()
  }

### extract data
  
  if (err == 0) {
    #if no error using mix
    tmp1[[fh]] <- summary(wgt2levs)$coef
    tmp4[[fh]] <- summary(wgt2only)$coef
    cfs1[fh,] <- tmp1[[fh]][, 1] #2 levels
    cfs4[fh,] <- tmp4[[fh]][, 1] #level 2 only
    ses1[fh,] <- tmp1[[fh]][, 2]
    ses4[fh,] <- tmp4[[fh]][, 2]
    icc1[fh] <- summary(wgt2levs.icc)$ICC
    icc4[fh] <- summary(wgt2only.icc)$ICC
    
    #using survey
    cfs5[fh, ] <- slev$coef
    ses5[fh, ] <- sqrt(diag(vcov(slev)))
    
    tmp2[[fh]] <- summary(wgt1a)$coef
    tmp3[[fh]] <- summary(nowgt1a)$coef
    
    cfs2[fh,] <- tmp2[[fh]][, 1]
    cfs3[fh,] <- tmp3[[fh]][, 1]
    
    ses2[fh,] <- tmp2[[fh]][, 2]
    ses3[fh,] <- tmp3[[fh]][, 2]
    
    icc2[fh] <- as.numeric(wgt1a.icc$ICC) #l1 only
    icc3[fh] <- as.numeric(nowgt1a.icc$ICC) #no weights
    
    # variance components
    
    vcomp1[fh, ] <- wgt2levs.icc$varDF[,4] #vcomp
    vcomp1se[fh, ] <- wgt2levs.icc$varDF[,7] #se of vcomp
    
    vcomp2[fh, ] <- wgt1a.icc$varDF[,4] #vcomp
    vcomp2se[fh, ] <- wgt1a.icc$varDF[,7] #se of vcomp
    
    vcomp3[fh, ] <- nowgt1a.icc$varDF[,4] #vcomp
    vcomp3se[fh, ] <- nowgt1a.icc$varDF[,7] #se of vcomp
    
    vcomp4[fh, ] <- wgt2only.icc$varDF[,4] #vcomp
    vcomp4se[fh, ] <- wgt2only.icc$varDF[,7] #se of vcomp
    
    
  }
  
  # Print the replication number
  cat("- Replication number", fh, '\n')
  
  ### testing
  min1[fh] <- mean(sampled_data$minority)
  min2[fh] <- weighted.mean(sampled_data$minority, sampled_data$totwgt)
  
}


#### COMPUTE SIMULATION RESULTS-- this is for one condition at a time

biasb1 <- (colMeans(cfs1, na.rm = T) - tbeta) / tbeta #weights at l1 and l2
biasb2 <- (colMeans(cfs2, na.rm = T) - tbeta) / tbeta #weight at l1
biasb3 <- (colMeans(cfs3, na.rm = T) - tbeta) / tbeta #no weight
biasb4 <- (colMeans(cfs4, na.rm = T) - tbeta) / tbeta #l2 weight only
biasb5 <- (colMeans(cfs5, na.rm = T) - tbeta) / tbeta #l2 weight only

tse1 <- apply(cfs1, 2, sd, na.rm = T)
tse2 <- apply(cfs2, 2, sd, na.rm = T)
tse3 <- apply(cfs3, 2, sd, na.rm = T)
tse4 <- apply(cfs4, 2, sd, na.rm = T)
tse5 <- apply(cfs5, 2, sd, na.rm = T)

biasse1 <- (colMeans(ses1, na.rm = T) - tse1) / tse1
biasse2 <- (colMeans(ses2, na.rm = T) - tse2) / tse2
biasse3 <- (colMeans(ses3, na.rm = T) - tse3) / tse3
biasse4 <- (colMeans(ses4, na.rm = T) - tse4) / tse4
biasse5 <- (colMeans(ses5, na.rm = T) - tse5) / tse5

# raw

rawcfs <- list("l2wgts" = tmp1, "l2only" = tmp4, "l1only" = tmp2, 
     'nowgt' = tmp3)
rawvc <- list('l2wgts.vc' = vcomp1, 'l2only.vc' = vcomp4, 'l1only.vc' = vcomp2,
              'nowgt.vc' = vcomp3,
              'l2wgts.vc.se' = vcomp1se, 'l2only.vc.se' = vcomp4se, 
              'l1only.vc.se' = vcomp2se,
              'nowgt.vc.se' = vcomp3se)

#coefs:::

biascfs <-
  round(
    data.frame(
      'two_wgts' = biasb1,
      'l1wgt' = biasb2,
      'no_wgt' = biasb3,
      'l2only' = biasb4,
      'slev' = biasb5
    ) * 100,
    3
  )
nms <- rownames(biascfs) #names of coefficients

biasse <- round(
  data.frame(
    'two_wgts' = biasse1,
    'l1wgt' = biasse2,
    'no_wgt' = biasse3,
    'l2only' = biasse4,
    'slev' = biasse5
  ) * 100,
  3
)
rownames(biasse) <- nms

esticc <-
  data.frame(
    'two_wgts' = mean(icc1, na.rm = T),
    'l1wgt' = mean(icc2),
    'no_wgt' = mean(icc3),
    'l2only' = mean(icc4, na.rm = T)
  )

# function to compute coverage probabilities
coverage <- function(b,
                     se,
                     true,
                     level = .95,
                     df = Inf) {
  qtile <- level + (1 - level) / 2 # Compute the proper quantile
  lower.bound <- b - qt(qtile, df = df) * se # Lower bound
  upper.bound <- b + qt(qtile, df = df) * se # Upper bound
  # Is the true parameter in the confidence interval? (yes = 1)
  true.in.ci <-
    ifelse(true >= lower.bound & true <= upper.bound, 1, 0)
  cp <- round(mean(true.in.ci), 4) # The coverage probability
  mc.lower.bound <-
    cp - 1.96 * sqrt((cp * (1 - cp)) / length(b)) # Monte Carlo error
  mc.upper.bound <- cp + 1.96 * sqrt((cp * (1 - cp)) / length(b))
  return(list(
    coverage.probability = round(cp, 4),
    mc.lb = mc.lower.bound,
    mb.ub = mc.upper.bound
  ))
}

#which have omitted data (nonconvergence):: only an issue with lower ICCs
wmiss <- which(is.na(cfs1[, 1]))

stats <- as.data.frame(summary_stats)
stats <- stats[!duplicated(stats[,1]), ]

l1n <- unlist(stats$number_of_participating_students)
l2n <- unlist(stats$number_of_participating_schools)

l2df <- l2n - 3
l1df <- l1n - 5 - l2n

# totalsch <- n_sample
# totaln <- nrow(sampled_data)
# l2df <- totalsch - 3 # J - no of cluster level predictors - 1
# l1df <- totaln - 5 - totalsch

df <- matrix(l1df, ncol = length(tbeta), nrow = reps)
df[, 1:3] <- l2df # rep(l2df, 4)
#df[1:4] <- l2df #degrees of freedom for coverage

k <- ncol(cfs1)
if (length(wmiss > 0)) {
  #if some reps did not converge
  cov1 <- sapply(1:k, function(x) {
    coverage(cfs1[-wmiss, x], ses1[-wmiss, x], tbeta[x], df = df[,x])
  })
  
  cov4 <- sapply(1:k, function(x) {
    coverage(cfs4[-wmiss, x], ses4[-wmiss, x], tbeta[x], df = df[,x])
  })
} else {
  #if all reps converged
  cov1 <- sapply(1:k, function(x) {
    coverage(cfs1[, x], ses1[, x], tbeta[x], df = df[,x])
  })
  
  cov4 <- sapply(1:k, function(x) {
    coverage(cfs4[, x], ses4[, x], tbeta[x], df = df[,x])
  })
  
}

cov2 <- sapply(1:k, function(x) {
  coverage(cfs2[, x], ses2[, x], tbeta[x], df = df[,x])
})

cov3 <- sapply(1:k, function(x) {
  coverage(cfs3[, x], ses3[, x], tbeta[x], df = df[,x])
})

cov5 <- sapply(1:k, function(x) {
  coverage(cfs5[, x], ses5[, x], tbeta[x], df = df[,x])
})

## coverage probabilities
covp <-
  round(
    data.frame(
      'two_wgts' = unlist(cov1[1, ]),
      'l1wgt' = unlist(cov2[1, ]),
      'no_wgt' = unlist(cov3[1, ]),
      'l2only' = unlist(cov4[1, ]),
      'svy' = unlist(cov5[1, ])
    ) * 100,
    1
  )
rownames(covp) <- nms

## results

biascfs #bias in coefficients
biasse #bias in standard errors
esticc #actual estimated ICCs
covp #coverage rates: we can use between 92.5 to 97.5 are acceptable (Bradley's condition)
#summary_stats # summary stats

cond$time <- Sys.time() - sta

res <- list(cond = cond,
 "biascfs" = biascfs , 
 "biasse" = biasse,
 "esticc" = esticc, 
 "covp" = covp,
 "summary_stats" = stats,
 "dm" = mean(ht < .05),
 "ps" = mean(ps < .05),
 "raw" = rawcfs,
 "df" = df,
 "tbeta" = tbeta,
 "m1" = mean(min1),
 "m2" = mean(min2),
 "min1" = min1, 
 "min2" = min2,
 'T00' = tau00,
 'S2' = sigma2,
 'rawvc' = rawvc #added variance components
 )

return(res)

}

# Example: NOT RUN
# This runs condition 3 for 50 reps with a seed of 321
# xx <- simwgt(3, 3, reps = 50, seed = 321)
# xx$covp
# xx$biascfs
# xx$biasse

