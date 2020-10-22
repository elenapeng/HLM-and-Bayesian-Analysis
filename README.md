# HLM-and-Bayesian-Analysis

Hierarchical Linear Models (HLM) used for the publication Peng et al. 2018 

- Contributor
```
Elena Peng
```
- Goals
```
Linear mixed models with an unstructured variance covariance structure and Tukey post hoc tests were used 
to evaluate the effect of muscle (VM and VMO), sex (male and female), and hip position (no hip rotation 
and 30 hip lateral rotation) on the two outcome variables, initial firing rate and recruitment threshold, 
during the isometric straight leg raise tasks. 
```
- Reason of using HLM
```
Mixed effects models were used because multiple motor units from an individual are statistically correlated 
and mixed models account for this correlation. 
```
- Model structure
```
In the linear mixed models for both outcome variables, the first level was single motor unit. Single motor 
units were nested according to each subject to form the second level, which was defined as the subject 
level. Muscle and hip position were the predictor variables for the motor unit level (level 1), while sex 
was the predictor variable for the subject level (level 2). 
```
- Equations for initial firing rate
```
Level 1: 

Initial Firing Rateij = β0j + β1jMuscleij + β2jHip positionij + β3jMuscle*Hip position + β4jForce + etij 

Level 2: 

β0j = γ00 + γ01Sex + u0j

β1j = γ10 + γ11Sex

β2j = γ20 + γ21Sex

β3j = γ30 + γ31Sex

β4j = γ40 

The absolute force at which a single motor unit was recruited was used as a covariate for initial firing 
rate due to the correlation between recruitment threshold and initial firing rate.  
```
- Equations for recruitment threshold
```
Level 1: 

Recruitment Thresholdij = β0j + β1jMuscleij + β2jHip positionij + β3jMuscle*Hip position + etij 

Level 2:

β0j = γ00 + γ01Sex + u0j 

β1j = γ10 + γ11Sex  

β2j = γ20 + γ21Sex

β3j = γ30 + γ31Sex

We did not use the averaged absolute force as a covariate because of the clear mathematical relation.
```
- Code Summary
```
Loading libraries and data, data pre-processing, assumption tests for multilevel model, HLM for initial 
firing rate and recruitment threshold, post-hoc analysis for interaction and main effects, box plots.
```
- Paper reference
```
Peng YL, Tenan MS, Griffin L. Hip position and sex differences in motor unit firing patterns of the 
vastus medialis and vastus medialis oblique in healthy individuals. J Appl Physiol (1985) 1;124(6):
1438-1446, 2018.
https://journals.physiology.org/doi/pdf/10.1152/japplphysiol.00702.2017
```
Bayesian Multilevel Modeling with brm for the paper in preparation

- Reason of using general HLM models in a Bayesian framework
```
The experimental design and outcome variables, initial firing rate and recruitment threshold, are similar 
to the previous HLM models. The difference is that we are interested in how knee pain affects the single 
motor unit performance. Therefore, the “Sex” factor was replaced by the “Group (healthy females and 
females with patellofemoral pain syndrome)” in the second level.

Bayesian mixed model approach allows parameters to vary through marginalization. We are trying to locate 
the relevant predictors for the outcome variables. Given our observed data, to generate the models with 
parameters that have 95% probability to lie within the credible region. 
```
