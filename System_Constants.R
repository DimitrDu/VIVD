# System specific constants ----------------------------------------------------
# AP, f, LipidFBS, XFBS, mAlbumin, pH, PSV, Temp, fSerum, SAMediumPlastic, V
# ------------------------------------------------------------------------------
SysConstants <- list(
  # Volume fractions - values from software, Certara)
  f = c('IW' = 0.568,               # fraction of intracellular water in  L
        'Lyso' = 0.01,              # fraction of lysosome in cells in L
        'Mito' = 0.1,               # fraction of mitochondria in cells in L
        'nl' = 0.0348,              # fraction of neutral lipids
        'np' = 0.0252),             # fraction of neutral phospholipid
  
  # Concentration of lipid in FBS in g/L (Ref 2); if LipidFBS is not available 
  # see Eq. 7 in Ref 1.
  LipidFBS = 0.46,  
  
  # Mass of albumin in g/L (From Ref 2)
  mAlbumin = 24, 
  
  # Partial specific volume
  PSV = c('Albumin' = 0.00073,      # of albumin in L/g (Ref 3)
          'TAG' = 0.001093),        # of lipids in L/g (Ref 4)
  
  # Temperature
  Temp = c('Ref' = 298.15,          # reference temperature in K
           'Sys' = 310.15),         # system temperature in K
  
  
  # Input variables ------------------------------------------------------------
  fSerum = 0.1,                     # Dimensionless fraction of serum in culture 
                                    # medium
  SAMediumPlastic = 0.000354,       # Area in m^2; see Eq (3) in Ref 1 if 
                                    # AMediumPlastic not available
  
  # Volume
  V = c('Medium' = 0.005,           # culture medium in L (= 1000 cm^3)
        'TotalCell' = 3e-6,         # total volume of cells in L; see Eq 1 if
                                    # VTotalCell is not available
        'Well' = 0.025)             # volume of a culture vessel in L

)


SysConstants$V[['Air']] <- with (
  SysConstants, 
  V[['Well']] - (V[['Medium']] + V[['TotalCell']])
)
