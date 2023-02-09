library(data.table)
library(ggplot2)

source('src/Universal_Constants.R')
source('src/System_Constants.R')
# source('src/Compound_Constants.R')
source('src/Initializations.R')
source('src/Partitioning.R')

cpds <- readxl::read_xlsx('./20221120_Compounds.xlsx')

results <- data.table()

for (i in 1:nrow(cpds)) {
  CpdConstants <- as.list(cpds[i,])
  
  CpdConstants <- init_Compound_Paras(UConstants, SysConstants, CpdConstants)
  fuFBSdilu <- binding_serum(SysConstants, CpdConstants)
  KPlastic <- partitioning_plastics(CpdConstants)
  KAir <- partitioning_air(UConstants, SysConstants, CpdConstants)
  K <- partitioning_cells(SysConstants, CpdConstants)
  
  K[['Air']] <- KAir
  K[['Plastic']] <- KPlastic
  
  C <- MassBalance(SysConstants, CpdConstants, fuFBSdilu, K)
  M <- perctMass(C, SysConstants)
  
  names(C) <- paste0('C_', names(C))
  names(M) <- paste0('M_', names(M))

  results <- rbind(results, cbind(as.data.table(cpds[i,]), 
                                  # as.data.table(t(K)), 
                                  as.data.table(t(C)), 
                                  as.data.table(t(M))))
}

writexl::write_xlsx(results, './results.xlsx')

# To plot Cell Concentration with logPow

results <- data.table()
for (i in seq(-200, 1000, 2)) {
  CpdConstants$log_Pow <- i/100.0
  
  CpdConstants <- init_Compound_Paras(UConstants, SysConstants, CpdConstants)
  fuFBSdilu <- binding_serum(SysConstants, CpdConstants)
  KPlastic <- partitioning_plastics(CpdConstants)
  KAir <- partitioning_air(UConstants, SysConstants, CpdConstants)
  K <- partitioning_cells(SysConstants, CpdConstants)
  K[['Air']] <- KAir
  K[['Plastic']] <- KPlastic

  C <- MassBalance(SysConstants, CpdConstants, fuFBSdilu, K)

  results <- rbind(results, cbind(as.data.table(CpdConstants), as.data.table(t(C))))
}

show(
  ggplot(results, aes(log_Pow, CellU)) +
    geom_point(color = 'red') + geom_hline(yintercept = 1e-9) +
   #  geom_point(aes(y = Plastic), color = 'blue', size = 1) +
    scale_y_log10(limits = c(1e-15, 1e-3)) + ggtitle('Cell Concentration vs logPow')
)
