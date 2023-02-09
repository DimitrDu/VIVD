# 1. PARTITIONING TO SERUM -----------------------------------------------------

binding_serum <- function (SysConstants, CpdConstants) {
  # * Eq. 6: Albumin-water partition coefficient: KProtein ---------------------
  if (CpdConstants$log_Pow < 4.5) {
    log_KProtein <- 1.08 * CpdConstants$log_PowTSys - 0.7
  } else {
    log_KProtein <- 0.37 * CpdConstants$log_PowTSys + 2.56
  }
  
  KProtein <- 10^(log_KProtein)
  
  # * Eq. 7: Neutral lipids in FBS: fnlFBS -------------------------------------
  fnlFBS <- SysConstants$LipidFBS * SysConstants$PSV[['TAG']]
  
  # * Eq. 5: Protein fraction: fProtein in L -----------------------------------
  fProtein <- SysConstants$mAlbumin * SysConstants$PSV[['Albumin']]
  
  # * Eq. 4: Fraction unbound in serum: fuFBS ----------------------------------
  fuFBS <- 1 / (1 + 
                  (KProtein * fProtein) + 
                  (CpdConstants$Pnl * fnlFBS)
                )

  # * Eq. 8, 9: Correction of fuFBS by dilution factor: D ----------------------
  fuFBSdilu <- fuFBS / (SysConstants$fSerum * (1 - fuFBS) + fuFBS)
  
  return (fuFBSdilu)
}


# 2. PARTITIONING TO PLASTICS --------------------------------------------------
# * Eq. 10 ---------------------------------------------------------------------

partitioning_plastics <- function (CpdConstants) {
  # KPlasticu in m
  K_Plasticu <- 10^(0.97 * CpdConstants$log_PowTSys - 6.94)
  
  return (K_Plasticu)
}


# 3. PARTIONING TO AIR ---------------------------------------------------------
# * Eq. 12: Partition coefficient to air: KAir ---------------------------------

partitioning_air <- function (UConstants, SysConstants, CpdConstants) {
  KAir <- CpdConstants$KH /(UConstants$R * SysConstants$Temp[['Ref']])
  
  # * Eq 13: Correction of KAir by the system temperature: KAiru ---------------
  log_KAiru <- log10(KAir) - (UConstants$deltaUaw / (log(10) * UConstants$R)) *
                  ((1 / SysConstants$Temp[['Sys']]) - 
                     (1 / SysConstants$Temp[['Ref']]))
  K_Airu <- 10^(log_KAiru)
  
  return (K_Airu)
}


# 4: PARTITIONING INTO CELLS ---------------------------------------------------
# Inputs: KaAP, KCellEWIW, KCellIWMito, KCellIWLyso, KCellu1, KCellu2, KCellu3
# Outputs: K['Cellu'], K['IW'], K['Lyso'], K['Mito']
partitioning_cells <- function (SysConstants, CpdConstants) {
  K <- c('CellU' = 0, 'CellIW' = 0, 'Mito' = 0, 'Lyso' = 0)
  K['CellU'] <- SysConstants$f[['IW']] + CpdConstants$Pnl * 
    SysConstants$f[['nl']] + CpdConstants$Pnp * SysConstants$f[['np']]

  # * Eq 1, 2, 3, Suppl.: Partitioning between extracellular water and cell 
  # compartments 
  K['CellIW'] <- K['Mito'] <- K['Lyso'] <- K[['CellU']]

  return (K)
}


MassBalance <- function(SysConstants, CpdConstants, fuFBSdilu, K) {
  C <- c('CellU' = 0, 'CellIW' = 0, 'Mito' = 0, 'Lyso' = 0, 'Air' = 0, 
         'Plastic' = 0)
  
  # 5. MASS BALANCE ------------------------------------------------------------
  # * Eq 20: Freely dissolved concentration in culture medium ------------------
  CMediumfree <- (CpdConstants$CNom * SysConstants$V[['Medium']]) / 
    (
      SysConstants$V[['Medium']] / fuFBSdilu + 
      K[['Air']] * SysConstants$V[['Air']] * CpdConstants$fui +
      K[['CellU']] * SysConstants$V[['TotalCell']] + 
      K[['Plastic']] * SysConstants$SAMediumPlastic * 1000
    )
  
  # * Eq 21 - 23, 4- 6 in Suppl.: Concentration in headspace, plastic, cells, 
  # intracellular water, lysosomes and mitochondria ----------------------------
  C <- K * CMediumfree
  C['MediumTotal'] <- CMediumfree / fuFBSdilu
  C['MediumFree'] <- CMediumfree
  C['MediumFBS'] <- CMediumfree * (1 - fuFBSdilu) / fuFBSdilu
  C['Air'] <- C['Air'] * CpdConstants$fui
  C['Plastic'] <- C['Plastic'] * 1000
  
  if(SysConstants$V[['Air']] == 0){
    C['Air'] <- 0
  }
  
  return (C)
}


perctMass <- function(C, SystemConstants) {
  M <- c()
  M['Nom'] <- CpdConstants$CNom * SysConstants$V[['Medium']]
  M['MediumTotal'] <- C[['MediumTotal']] * SysConstants$V[['Medium']]  / M[['Nom']] * 100
  M['MediumFree'] <- C[['MediumFree']] * SysConstants$V[['Medium']] / 
    M[['Nom']] * 100
  M['MediumFBS'] <- C[['MediumFBS']] * SysConstants$V[['Medium']] /
    M[['Nom']] * 100
  M['Air'] <- C[['Air']] * SysConstants$V[['Air']] / M[['Nom']] * 100
  M['Plastic'] <- C[['Plastic']] * SysConstants$SAMediumPlastic / M['Nom'] * 100
  
  
  M['CellT'] <- C[['CellU']] * SysConstants$V[['TotalCell']] / M['Nom'] * 100
  M['Lyso'] <- C[['Lyso']] * SysConstants$V[['TotalCell']] * 
    SysConstants$f[['Lyso']] / M[['Nom']] * 100
  M['Mito'] <- C[['Mito']] * SysConstants$V[['TotalCell']] * 
    SysConstants$f[['Mito']] / M[['Nom']] * 100
  M['CellIW'] <- M['CellT'] - M['Lyso'] - M['Mito']
  
  return (M)
}