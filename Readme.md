# Virtual in vitro distribution model (VIVD) 
Date: 20/11/2022

Authors: Dunja Dimitrijevic, Dr. Varun Giri

The objective of this study was the prediction of the distribution of test
substances with varying physicochemical properties in the compartments 
contained in a typical in vitro test system, i.e. medium and its ingredients,
head space, cell layer, culture vessel. In the following R Scripts, the 
distribution of test substances within an in vitro system is presented in 
line with the equations and model published mass balance model by Fisher et 
al., 2019. In the current version, the model is valid for neutral substances.
The model was further validated with published data from Armitage et al., 2014.

The described model is divided in the below presented main sections:

Universal Constants:
This R Script contains generic, universal values, e.g. the universal gas 
constant R.
 
System Constants:
In this R Script, the in vitro system is characterized with specific 
information about the test system. The relevant parameters are defined by the
user.

Compound Constants:
In this file, relevant parameters of the test compound are provided. The 
relevant parameters are defined by the user. If more than one test substance 
need to be considered, a separate Excel file with relevant information of the
test substances can be included and read by the "execute model" file. 
 
Initialization:
Several pre-calculation steps, e.g. partitioning to lipids, are performed in 
this file (e.g. inclusion of ionization of substances). In this section, the
ionization term of substances and its influence on the uptake into different
cellular compartments based on the pH value and ionization state can be
considered. As mentioned above, the model is at the current state only 
validated for neutral test substances.

Partitioning:
The partitioning to the main constituents within the in vitro system, i.e. 
proteins and lipids deriving from the culture medium, plastic vessels, 
head space and cells.

Execute model:
The R Script file "Execute model" performs the actual calculation including 
all the equations from the Script files presented above. Furthermore, an 
output file with the results is generated.

Furthermore, several parameters have been included from different 
publications. In the R script files, following abbreviations refer to the 
corresponding reference:

Ref 1: Fisher et al., 2019, DOI: 10.1016/j.tiv.2018.12.017.
Ref 2: Armitage et al., 2014, DOI: 10.1021/es501955g.
Ref 3: Kupke et al., 1972, DOI: 10.1073/pnas.69.8.2258.
Ref 4: Deckelbaum et al., 1984, DOI: 10.1161/01.atv.4.3.225.
Ref 5: Rodgers et al., 2005, DOI: 10.1002/jps.20322.
Ref 6: Rodgers and Rowland, 2006; DOI: DOI:10.1002/jps.20502.      
Ref 7: Trapp et al., 2008, DOI: 10.1007/s00249-008-0338-4.
Ref 8: Glunde et al., 2003; DOI: 10.1016/s1476-5586(03)80037-4.
