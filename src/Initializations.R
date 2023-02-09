# Initialization steps

init_Compound_Paras <- function (UConstants, SysConstants, CpdConstants) {

  # * Eq. 19: Correction of logPow by system temperature TSys ------------------
  # log(10) - to convert all log values base 10 (log_PowTSys and log_Pow) to 
  # natural log. In case the value is already in natural log, then we do not 
  # need this factor.
  CpdConstants$log_PowTSys <- CpdConstants$log_Pow -
    (UConstants$deltaUow / (log(10) * UConstants$R)) * 
    ((1 / SysConstants$Temp[['Sys']]) - 
       (1 / SysConstants$Temp[['Ref']]))
  
    
  # Partitioning between IW and neutral phospholipid defined by octanol to water
  # partition coefficent, dimensionless (Ref 5)
  # CpdConstants$Pnp <- 10^(CpdConstants$log_Pow)
  CpdConstants$Pnp <- 10^(CpdConstants$log_PowTSys)
  
  # Eq 13, 14 (Ref 9): Calculation of partitioning in neutral lipids -----------
  if (CpdConstants$OxygenMol == T){
    # log_Pvow <- 1.099 * CpdConstants$log_Pow - 1.31
    log_Pvow <- 1.099 * CpdConstants$log_PowTSys - 1.31
  } else if (CpdConstants$OxygenMol == F){
    # log_Pvow <- 1.0654 * CpdConstants$log_Pow - 0.232
    log_Pvow <- 1.0654 * CpdConstants$log_PowTSys - 0.232
  } else {
    stop ('Incorrect Moltype.')
  }
  
  # log_Pvow <- 1.115 * CpdConstants$log_PowTSys - 1.35
  
  # Partitioning between IW and neutral lipids, defined by olive oil to water 
  CpdConstants$Pnl <- 10^(log_Pvow)
  # partition coefficient, dimensionless (Ref 5)

  return(CpdConstants)
}