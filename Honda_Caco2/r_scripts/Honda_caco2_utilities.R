
#### Model Summary Utils ####
get_descriptors <- function(option.desc = c("mordred", "padel")){ 
 # loc.wd <- getwd()
  #load(paste0(loc.wd,"/r_data/processed/chembl_14FEB2019.RData"))
  
  ##### Initialize ####
  load(paste0(loc.wd, "/r_data/chemprop/gut_chems_all.RData")) # gut.chems.all; dtxsid, qsar_ready_smiles, matched to reference chem names
  load(paste0(loc.wd, "/r_data/chemprop/combined_chemprop.RData")) # chemprops.dt
  
  # Determine ionization and logD
  logD.dt <- unique(chemprop.dt[,
                                .(dtxsid, 
                                  logP = LogP_use, 
                                  pKa_a = pKa_a_use, 
                                  pKa_b = pKa_b_use)]) %>% 
    .[,c("fraction_neutral",
         "fraction_charged",
         "fraction_negative",
         "fraction_positive",
         "fraction_zwitter") := calc_ionization(pKa_Donor = pKa_a,
                                                   pKa_Accept = pKa_b,
                                                   pH = 7.4), .(pKa_a,pKa_b)] %>% 
    .[, logD_opera := log10(calc_dow(Pow = 10^logP, fraction_charged = fraction_charged)), .(logP, fraction_charged)] %>% 
    .[, logP_opera := logP] %>% 
    .[, fneutral_logit := logit(fraction_neutral, delta = 1e-6)] %>% 
    .[, fcharged_logit := logit(fraction_charged, delta = 1e-6)] %>% 
    .[, fpositive_logit := logit(fraction_positive, delta = 1e-6)] %>% 
    .[, fnegative_logit := logit(fraction_negative, delta = 1e-6)] %>% 
    .[, fzwitter_logit := logit(fraction_zwitter, delta = 1e-6)] %>% 
    .[, c("logP","pKa_a", "pKa_b") := NULL]
  
  
  # Load mordred descriptors
  load(paste0(loc.wd, "/r_data/descriptors/mordred_desc_26MAR2019.RData")) # mordred.desc
  mordred.names <- names(mordred.desc)[!(names(mordred.desc) %in% c("qsar_ready_smiles", "casrn", "dtxsid"))]
  setnames(mordred.desc,c("qsar_ready_smiles", mordred.names),c("smiles",fixnames(mordred.names)))
  mordred.names <- names(mordred.desc)[!(names(mordred.desc) %in% c("smiles", "casrn", "dtxsid"))]
  mordred.desc2 <- cbind(mordred.desc[,lapply(mordred.names, function(x) as.numeric(eval(parse(text = x))))],
                         mordred.desc[,.(smiles,casrn,dtxsid)])
  setnames(mordred.desc2, names(mordred.desc2)[1:length(mordred.names)],mordred.names)
  mordred.desc2 <- logD.dt[mordred.desc2, on = "dtxsid"]
  mordred.desc2 <- mordred.desc2[, c("dtxsid","casrn") := NULL] 
#  mordred.desc2 <- mordred.desc2[,dtxsid := NULL] %>% 
#    .[,casrn:= NULL] %>% 
#    unique()
  
  # Load PaDEL descriptors
  load(paste0(loc.wd, "/r_data/descriptors/padel_desc_26MAR2019.RData")) # mordred.desc
  padel.names <- names(padel.desc)[!(names(padel.desc) %in% c("qsar_ready_smiles","casrn","dtxsid","chem_id"))]
  setnames(padel.desc, c("qsar_ready_smiles", padel.names), c("smiles", fixnames(padel.names)))
  padel.names <- names(padel.desc)[!(names(padel.desc) %in% c("smiles","casrn","dtxsid","chem_id"))]
  padel.desc2 <- cbind(padel.desc[,lapply(padel.names, function(x) as.numeric(eval(parse(text = x))))],
                       padel.desc[,.(smiles,casrn,dtxsid,chem_id)])
  setnames(padel.desc2, names(padel.desc2)[1:length(padel.names)],padel.names)
  padel.desc2 <- padel.desc2[, c("dtxsid","casrn","chem_id") := NULL] 
  #%>% 
   # unique()
  # Drop repeats
  padel.desc2[, X := 1:length(smiles)] %>% 
    .[,X2 := 1:length(X),.(smiles)]
  padel.desc2 <- padel.desc2[X2 == 1,]
  padel.desc2[,c("X","X2") := NULL]
  
  # Combine mordred & padel, prioritize mordred descriptors  
  unique.padel.nm <- paste0(padel.names[!padel.names %in% mordred.names],"_padel")
  setnames(padel.desc2,padel.names,paste0(padel.names,"_padel"))
  unique.padel.nm.smi <- c("smiles",unique.padel.nm)
  
  # Choose descriptor set
  if(all(option.desc %in% c("mordred","padel"))){
    desc.dt <- merge(padel.desc2,subset(mordred.desc2,!(is.na(TopoPSA))),by="smiles")
  }
  if(all(option.desc %in% c("mordred"))){
    desc.dt <- copy(mordred.desc2)
  }
  
  return(desc.dt)  
}

calc_catstat <- function(ymeas, ycalc){
  temp <- as.data.table(melt(table(ycalc,ymeas), value.name = "ncount"))
  highacc <- temp[ycalc == "high" & ymeas == "high", sum(ncount)]/temp[ymeas == "high", sum(ncount)]
  lowacc <- temp[ycalc == "low" & ymeas == "low", sum(ncount)]/temp[ymeas == "low", sum(ncount)]
  balacc <- (highacc + lowacc)/2
  nright <- temp[ymeas == ycalc, sum(ncount)]
  ntot <- temp[, sum(ncount)]
  nacc <- nright/ntot
  randprob <- (temp[ymeas == "high", sum(ncount)]*temp[ycalc == "high", sum(ncount)] +
                 temp[ymeas == "low", sum(ncount)]*temp[ycalc == "low", sum(ncount)])/ntot^2
  cohen.kappa <- (nacc-randprob)/(1-randprob)
  return(list("nacc" = nacc, "ntot" = ntot,"cohen.kappa" = cohen.kappa,"balacc" = balacc,  "lowacc" = lowacc, "highacc" = highacc))
}

get_rfres <- function(this.model){
  this.model$results[as.numeric(rownames(this.model$bestTune)),]
}

get_cattrain <- function(this.dt){
  temp <- this.dt[,as.list(c(quantile(mpab, c(0.05,0.95)),
                            "med.pab" = median(mpab, na.rm = TRUE),
                            "mean.pab" = mean(mpab, na.rm = TRUE),
                            "sd.pab" = sd(mpab, na.rm = TRUE))),
                  .(mpab.cat)] %>% 
    melt(., id.vars = "mpab.cat") %>% 
    .[, variable := paste0("cat",mpab.cat, "_",gsub("%","perc", variable))] %>% 
    .[, mpab.cat := NULL] %>% 
    dcast(... ~ variable) %>% 
    .[, "." := NULL]
  return(temp)
  
}

build_qsar_desc <- function(alldesc.dt, this.caco2){
  padel.nm <- names(alldesc.dt)[regexpr("padel",names(alldesc.dt))!=-1]
  mordred.nm <- names(alldesc.dt)[!names(alldesc.dt) %in% padel.nm &
                                    !names(alldesc.dt) %in% c("smiles",
                                                              "fraction_neutral", "fraction_charged", "fraction_negative", "fraction_positive", "fraction_zwitter",
                                                              "fneutral_logit", "fcharged_logit", "fnegative_logit", "fpositive_logit", "fzwitter_logit",
                                                              "logD_opera","logP_opera")]
  qsar.dt <- desc.dt[this.caco2[,.(setx, smiles, mpab, mpab.cat,
                                   fraction_neutral, fraction_charged, fraction_negative, fraction_positive, fraction_zwitter,
                                   fneutral_logit, fcharged_logit, fnegative_logit, fpositive_logit, fzwitter_logit,
                                   logD_opera, logP_opera)],on="smiles"]
  
  
}


#### random forest model ####
create_qsardt <- function(alldesc.dt, 
                          this.caco2, 
                          option.desc = "mordred"){
  padel.nm <- names(alldesc.dt)[regexpr("padel",names(alldesc.dt))!=-1]
  mordred.nm <- names(alldesc.dt)[!names(alldesc.dt) %in% padel.nm &
                                    !names(alldesc.dt) %in% c("smiles",
                                                              "fraction_neutral", "fraction_charged", "fraction_negative", "fraction_positive", "fraction_zwitter",
                                                              "fneutral_logit", "fcharged_logit", "fnegative_logit", "fpositive_logit", "fzwitter_logit",
                                                              "logD_opera","logP_opera")]
  
  # Choose descriptor set
  if(all(c("mordred","padel") %in% option.desc)){
    xnames <- c(padel.nm, mordred.nm)
    use.names <- c("smiles", xnames)
    desc.dt <- unique(alldesc.dt[,..use.names])
    desc_out <- "mord_padl"
  }
  if(all(option.desc %in% "mordred")){
    xnames <- c(mordred.nm)
    use.names <- c("smiles", xnames)
    desc.dt <- unique(alldesc.dt[,..use.names])
    desc_out <- "mord"
  }
  if(all(option.desc %in% "padel")){
    xnames <- c(padel.nm)
    use.names <- c("smiles", xnames)
    desc.dt <- copy(alldesc.dt)
    desc.dt <- unique(desc.dt)
    desc_out <- "padl"
  }
  
  # Combine with descriptors
  qsar.dt <- desc.dt[this.caco2, on = "smiles"]
  xnames <- c(xnames, c("fraction_neutral", "fraction_charged", "fraction_negative", "fraction_positive", "fraction_zwitter",
                        "fneutral_logit", "fcharged_logit", "fnegative_logit", "fpositive_logit", "fzwitter_logit",
                        "logD_opera","logP_opera"))
  return(list("qsar.dt" = qsar.dt, "xnames" = xnames))
  
}

prep_caco2 <- function(caco2_qc_file = "caco2_QC_rec0d4to2",
                       prep_options = list(
                         "pab.breaks" = c(0.9),
                         "pab.lims" = c(-2,2.5),
                         "efflux.break" = 0.5,
                         "efflux.train.size" = 10,
                         "option.gap" = FALSE,
                         "option.round" = "ten",
                         "option.flatten" = NA,
                         "option.drop.high_efflux" = FALSE,
                         "train.size" = NA,
                         "nval.smiles" = 50)){
  
  pab.breaks <- prep_options[["pab.breaks"]]
  pab.lims <- prep_options[["pab.lims"]]
  efflux.break <- prep_options[["efflux.break"]]
  efflux.train.size <- prep_options[["efflux.train.size"]]
  option.gap <- prep_options[["option.gap"]]
  option.round <- prep_options[["option.round"]]
  option.flatten <- prep_options[["option.flatten"]]
  option.drop.high_efflux <- prep_options[["option.drop.high_efflux"]]
  train.size <- prep_options[["train.size"]]
  nval.smiles <- prep_options[["nval.smiles"]]
  #loc.wd <- getwd()
  
  
  #### Initialize ####
  {
    load(paste0(loc.wd, "/r_data/processed/all_gut_data.RData")) # gut.data; compiled data table of absorption and permeability data 
    load(paste0(loc.wd, "/r_data/chemprop/gut_chems_all.RData")) # gut.chems.all; dtxsid, qsar_ready_smiles, matched to reference chem names
    load(paste0(loc.wd, "/r_data/chemprop/combined_chemprop.RData")) # chemprops.dt
    
    # Load our Caco-2 data 
    load(paste0(loc.wd,"/r_data/caco2_qc/",caco2_qc_file,".RData")) # caco2.dt
    caco2.dt[,dtxsid := NULL]
    caco2.dt <- unique(gut.chems.all[repref == "honda_caco2", 
                                     .(dtxsid, casrn = ref.casrn, qsar_ready_smiles)])[caco2.dt, on = "casrn"]
    # Check relationship with recovery (mass balance)
    # ggplot(caco2.dt)+
    #   geom_point(aes(x = rec_ab, y = log10(Pab)))
    caco2.dt[, Pab := Pab/rec_ab] %>%  # Normalize with respect to recovery
      .[, Pba := Pba/rec_ba] %>% 
      .[, efflux_ratio := Pba/Pab]
    
    # Determine ionization and logD
    caco2.dt <- unique(chemprop.dt[,
                                   .(dtxsid, 
                                     logP = LogP_use, 
                                     pKa_a = pKa_a_use, 
                                     pKa_b = pKa_b_use)])[caco2.dt, on = "dtxsid"]
    
    caco2.dt[,c("fraction_neutral",
                "fraction_charged",
                "fraction_negative",
                "fraction_positive",
                "fraction_zwitter") := calc_ionization(pKa_Donor = pKa_a,
                                                          pKa_Accept = pKa_b,
                                                          pH = 7.4), .(pKa_a,pKa_b)] %>% 
      .[,logD := log10(gh_calc_dow(Pow = 10^logP, fraction_charged = fraction_charged)), .(logP, fraction_charged)]
    
    # Load Lanevskij and Obringer data
    load(paste0(loc.wd, "/r_data/processed/lit_caco2_26MAR2019.RData")) # lit.caco2.dt; Lanevskij and Obringer
    # Retain measured values with pH near our pH
    lit.caco2.dt[,near74 := (abs(lit_exp_pH - 7.4) == min(abs(lit_exp_pH - 7.4))), .(dtxsid)]
    # Determine ionization and logD
    lit.caco2.dt <- unique(chemprop.dt[,
                                       .(dtxsid, 
                                         logP = LogP_use, 
                                         pKa_a = pKa_a_use, 
                                         pKa_b = pKa_b_use)])[lit.caco2.dt, on = "dtxsid"]
    
    lit.caco2.dt[,c("fraction_neutral",
                    "fraction_charged",
                    "fraction_negative",
                    "fraction_positive",
                    "fraction_zwitter") := gh_calc_ionization(pKa_Donor = pKa_a,
                                                              pKa_Accept = pKa_b,
                                                              pH = lit_exp_pH), .(pKa_a, pKa_b, lit_exp_pH)] %>% 
      .[,logD := log10(gh_calc_dow(Pow = 10^logP, fraction_charged = fraction_charged)),.(logP, fraction_charged)]
    
    # Load pkstats from Wambaugh et al. 2018
    load(paste0(loc.wd, "/r_data/tk_data/PKstats-2018-01-16.RData")) #
    pkstats3 <- unique(gut.chems.all[,.("CAS" = casrn, dtxsid, qsar_ready_smiles)])[PKstats, on = "CAS"]
    reserve.smiles <- unique(pkstats3[, qsar_ready_smiles])
    
    
    
    
  }
  
  #### pharma-pab ####
  pharma.pab <- lit.caco2.dt[near74 == TRUE,
                             .(mpab = mean(log10(lit_pab),na.rm = T),
                               sdpab = sd(log10(lit_pab),na.rm = T),
                               logD_opera = mean(logD,na.rm = T),
                               logP_opera = logP[1],
                               fraction_neutral = fraction_neutral[1],
                               fraction_charged = fraction_charged[1],
                               fraction_negative = fraction_negative[1],
                               fraction_positive = fraction_positive[1],
                               fraction_zwitter = fraction_zwitter[1]),
                             .(smiles = qsar_ready_smiles)] %>% 
    .[!is.na(mpab),]
  pharma.pab[,fneutral_logit := logit(fraction_neutral, delta = 1e-6)] %>% 
    .[, fcharged_logit := logit(fraction_charged, delta = 1e-6)] %>% 
    .[, fpositive_logit := logit(fraction_positive, delta = 1e-6)] %>% 
    .[, fnegative_logit := logit(fraction_negative, delta = 1e-6)] %>% 
    .[, fzwitter_logit := logit(fraction_zwitter, delta = 1e-6)]
  
  set.seed(5302019)
  validation.smiles <- pharma.pab[!(smiles %in% unique(caco2.dt$smiles)) & 
               !(smiles %in% reserve.smiles) & !is.na(mpab), sample(smiles, nval.smiles, replace = FALSE)]
  reserve.smiles <- c(reserve.smiles, validation.smiles)
  # Define Training set
  pharma.pab2 <- def.trainx(pab.dt = pharma.pab, 
                            reserve.smiles = reserve.smiles,
                            pab.breaks = pab.breaks,
                            pab.lims = pab.lims,
                            efflux.break = NA,
                            byref = FALSE,
                            option.gap = option.gap,
                            option.round = option.round,
                            option.flatten = option.flatten,
                            train.size = NA)
  
  ##### pab-combined #####
  rcaco2.med <- caco2.dt[!is.na(Pab) & !is.na(qsar_ready_smiles),
                         .(mpab=mean(log10(Pab),na.rm=T),
                           sdpab = sd(log10(Pab),na.rm=T),
                           mpba = mean(log10(Pba),na.rm = T),
                           sdpba = sd(log10(Pba),na.rm = T),
                           mrecab=mean((rec_ab),na.rm=T),
                           sdrecab = sd((rec_ab),na.rm=T),
                           mrecba = mean((rec_ba),na.rm = T),
                           sdrecba = sd((rec_ba),na.rm = T),
                           mefflux = mean(log10(efflux_ratio),na.rm = T),
                           sdefflux = sd(log10(efflux_ratio),na.rm = T),
                           logD_opera = mean(logD,na.rm = T),
                           logP_opera = logP[1],
                           fraction_neutral = fraction_neutral[1],
                           fraction_charged = fraction_charged[1],
                           fraction_negative = fraction_negative[1],
                           fraction_positive = fraction_positive[1],
                           fraction_zwitter = fraction_zwitter[1]),
                         .(smiles = qsar_ready_smiles)]
  rcaco2.med[,fneutral_logit := logit(fraction_neutral, delta = 1e-6)] %>% 
    .[, fcharged_logit := logit(fraction_charged, delta = 1e-6)] %>% 
    .[, fpositive_logit := logit(fraction_positive, delta = 1e-6)] %>% 
    .[, fnegative_logit := logit(fraction_negative, delta = 1e-6)] %>% 
    .[, fzwitter_logit := logit(fraction_zwitter, delta = 1e-6)]
  
  this.ourcaco2 <- def.trainx(pab.dt = rcaco2.med, 
                              reserve.smiles = reserve.smiles,
                              efflux.break = efflux.break,
                              efflux.train.size = efflux.train.size,
                              pab.breaks = pab.breaks,
                              pab.lims = pab.lims,
                              byref = FALSE,
                              option.gap = option.gap,
                              option.round = option.round,
                              option.flatten = option.flatten,
                              option.drop.high_efflux = option.drop.high_efflux,
                              train.size = NA)
  
  caco2.comb <- rbind(copy(rcaco2.med)[, c("setx", "mpab.cat") := NULL][,repref := "honda"],
                      copy(pharma.pab[!(smiles %in% rcaco2.med$smiles),])[, c("setx", "mpab.cat") := NULL][, repref := "lit"], fill = TRUE)
  caco2.comb[repref == "lit", mefflux := 0.00001]
  caco2.comb2 <- def.trainx(pab.dt = caco2.comb,
                            reserve.smiles = reserve.smiles,
                            pab.breaks = pab.breaks,
                            pab.lims = pab.lims,
                            efflux.break = efflux.break,
                            efflux.train.size = efflux.train.size,
                            byref = TRUE,
                            option.gap = option.gap,
                            option.round = option.round,
                            option.flatten = option.flatten,
                            option.drop.high_efflux = option.drop.high_efflux,
                            train.size = train.size)
  
  return(list("this.litcaco2" = pharma.pab2,
              "this.ourcaco2" = this.ourcaco2,
              "this.combcaco2" = caco2.comb2,
              "options.caco2" = list("pab.breaks" = pab.breaks,
                                     "pab.lims" = pab.lims,
                                     "")))
  
}



# Build randomforest classification and regression models
rfmodel_build <- function(qsar.dt,
                          desc.names,
                          option.filter.chems = TRUE,
                          option.model = "pab.cat",
                          pval.cutoff = 0.1) # Optional, removes chemical with the most missing (NA or Inf) descriptors
{
  # Define seedval for rf functions
  seedval <- 22102018
  
  # Define randomforest function for caret
  # rfFuncs2 <- rfFuncs
  # rfFuncs2$fit <- function (x, y, first, last, ...) {
  #   loadNamespace("randomForest")
  #   randomForest::randomForest(x, y, ...)
  # }
  # 
  # Filter chemicals with missing values (all sets)
  if(option.filter.chems == TRUE){
    print("Filtering chemicals with missing values")
    temp <-unique(which(is.na(as.matrix(qsar.dt[, ..desc.names])) |
                          is.inf(as.matrix(qsar.dt[, ..desc.names])), arr.ind = TRUE)[, 1])
    tempa <- qsar.dt[, ..desc.names]
    temp2 <- data.table(x = temp)[, nalength := sum(is.na(unlist(tempa[x,])) |
                                                      is.inf(unlist(tempa[x,]))), .(x)]
    if(length(temp2$nalength) > 0 && 
       (max(temp2$nalength) > min(temp2$nalength))){
      skiprow <- temp2[nalength == max(nalength), x]
      qsar.dt2 <- qsar.dt[-skiprow,]
    }else{
      skiprow <- NA
      qsar.dt2 <- copy(qsar.dt)
    }
  }else{
    skiprow <- NA
    qsar.dt2 <- copy(qsar.dt)
  }
  
  # Filter descriptors with missing values  
  skipcol <- unique(which(is.na(qsar.dt2[, ..desc.names]) | 
                            is.inf(qsar.dt2[, ..desc.names]), arr.ind = TRUE)[, 2])
  if(length(skipcol) > 0){
    desc.usenames <- desc.names[-skipcol]
  }else{
    skipcol <- NA
    desc.usenames <- copy(desc.names)
  }
  
  # Random Forest Classification Models
  if(any(c("pab.cat","pab.cat.glm") %in% option.model)){
    
    # Define qsar recipes
    pab.cat.dt <- as.data.frame(data.table(qsar.dt2[setx == "trainx",.(mpab.cat)],qsar.dt2[setx == "trainx",..desc.usenames]))

    qsar.recipe.cat <- recipe(mpab.cat ~ ., pab.cat.dt) %>% 
      step_nzv(all_predictors()) %>% 
      step_corr(all_predictors(), threshold = 0.9) 
    qsar.recipe.cat.glm <- recipe(mpab.cat ~ ., pab.cat.dt) %>% 
      step_nzv(all_predictors()) %>% 
      step_corr(all_predictors(), threshold = 0.9) %>% 
      step_glm(all_predictors(), all_outcomes(), yfamily = "binomial", pval_cutoff = pval.cutoff)
    
    # Train models    
    if("pab.cat" %in% option.model){
      print("Training classification model, nzv & highcor filter")
      set.seed(seed = seedval)
      rf.pab.cat <- train(qsar.recipe.cat, pab.cat.dt,
                          model = "rf",
                          tuneLength = 5,
                          trControl = trainControl(method = "cv",
                                                   number = 5,
                                                   classProbs = TRUE,
                                                   savePredictions = "final"),
                          verbose = TRUE,
                          importance = TRUE)
    }else{
      rf.pab.cat <- list()
    }
    print(rf.pab.cat)
    
    if("pab.cat.glm" %in% option.model){
    print("Training classification model, nzv, highcor, & glm filter")
    set.seed(seed = seedval)
    rf.pab.cat.glm <- train(qsar.recipe.cat.glm, pab.cat.dt,
                            model = "rf",
                            tuneLength = 5,
                            trControl = trainControl(method = "cv",
                                                     number = 5,
                                                     classProbs = TRUE,
                                                     savePredictions = "final"),
                            verbose = TRUE,
                            importance = TRUE)
    }else{
      rf.pab.cat.glm <- list()
    }
    print(rf.pab.cat.glm)
  
  }else{
    rf.pab.cat <- list()
    rf.pab.cat.glm <- list()
  }
  
  # Random Forest Regression Models
  if(any(c("pab.reg","pab.reg.glm") %in% option.model)){
    
    # Define recipe inputs
    pab.reg.dt <- as.data.frame(data.table(qsar.dt2[setx == "trainx",.(mpab)],qsar.dt2[setx == "trainx",..desc.usenames]))
    qsar.recipe.reg <- recipe(mpab ~ ., pab.reg.dt) %>% 
      step_nzv(all_predictors()) %>% 
      step_corr(all_predictors(), threshold = 0.9) 
    qsar.recipe.reg.glm <- recipe(mpab ~ ., pab.reg.dt) %>% 
      step_nzv(all_predictors()) %>% 
      step_corr(all_predictors(), threshold = 0.9) %>% 
      step_glm(all_predictors(), all_outcomes(), yfamily = "gaussian", pval_cutoff = pval.cutoff)
    
    # Train regression models
    if("pab.reg" %in% option.model){
      print("Training regression model, nzv & highcor filter")
      set.seed(seed = seedval)
      rf.pab.reg <- train(qsar.recipe.reg, pab.reg.dt,
                          model = "qrf",
                          tuneLength = 5,
                          trControl = trainControl(method = "cv",
                                                   number = 5,
                                                   savePredictions = "final"),
                          verbose = TRUE,
                          importance = TRUE,
                          keep.inbag = TRUE)
    }else{
      rf.pab.reg <- list()
    }
    print(rf.pab.reg)
    
    if("pab.reg.glm" %in% option.model){
      print("Training regression model, nzv, highcor & glm filter")
      set.seed(seed = seedval)
      rf.pab.reg.glm <- train(qsar.recipe.reg.glm, pab.reg.dt,
                              model = "qrf",
                              tuneLength = 5,
                              trControl = trainControl(method = "cv",
                                                       number = 5,
                                                       savePredictions = "final"),
                              verbose = TRUE,
                              importance = TRUE,
                              keep.inbag = TRUE)
    }else{
      rf.pab.reg.glm <- list()
    }
    print(rf.pab.reg.glm)
  }else{
    rf.pab.reg <- list()
    rf.pab.reg.glm <- list()
  }
  
  pab.res <- list(rf.prefilter = list(input.dt = qsar.dt,
                                      use.dt = qsar.dt2,
                                      desc.usenames = desc.usenames,
                                      skipcol = skipcol,
                                      skiprow = skiprow),
                  rf.pab.cat = rf.pab.cat,
                  rf.pab.cat.glm = rf.pab.cat.glm,
                  rf.pab.reg = rf.pab.reg,
                  rf.pab.reg.glm = rf.pab.reg.glm)
  return(pab.res)
  
}


# Define training set
def.trainx <- function(pab.dt,
                       reserve.smiles,
                       efflux.break =  NA,
                       efflux.train.size = 10,
                       pab.breaks = c(0.25,1), # Length 1 or 2 (2 values or 3 values)
                       pab.lims = c(-3,3), # Lower and Upper limits for mpab
                       byref = TRUE,
                       option.gap = TRUE,
                       option.flatten = NA, # Use sub-bins (quartiles) to flatten distributions
                       option.round = "one",
                       option.drop.high_efflux = FALSE,
                       train.size = NA) 
{
  
  # prep_options <- caco2_prep_options[[1]]
  # caco2_qc_file <- "caco2_QC_rec0d4to2"
  # 
  # pab.dt <- copy(caco2.comb)
  # reserve.smiles <- reserve.smiles
  # pab.breaks <- 0.6
  # pab.lims <- pab.lims
  # byref <- TRUE
  # option.gap <- option.gap
  # option.round <- option.round
  # train.size <- train.
  # 
  
  seedval <- 10192018
  seedval2 <- 5212019
  if(("repref" %in% names(pab.dt)) & byref == TRUE){
    pab.dt[,this.ref := repref]
  }else{
    pab.dt[,this.ref := "thisref"]
  }
  pab.dt[,setx := as.character(NA)] %>% 
    .[smiles%in%reserve.smiles,setx:="reservex"] # Reserve smiles from pkstats

  if(!is.na(efflux.break)){
    set.seed(seedval2)
    pab.dt[!is.na(mefflux) & mefflux <= efflux.break, mefflux.cat := "low"] %>% 
      .[!is.na(mefflux) & mefflux > efflux.break, mefflux.cat := "high"]
    efflux.trainx <- pab.dt[!is.na(mefflux),.("smiles" = sample(smiles, efflux.train.size, replace = FALSE)) ,.(mefflux.cat)]
    pab.dt[smiles %in% efflux.trainx[, smiles], setx := "trainx.efflux"]
    if(option.drop.high_efflux == TRUE){
      
      pab.dt[!is.na(mefflux.cat) & mefflux.cat == "high", setx := "high.efflux"]
      
    }
  }
  
  if(length(pab.breaks) == 2){
    
    pab.dt[mpab < pab.breaks[1] & mpab >= pab.lims[1], mpab.cat := "low"] %>% 
      .[mpab >= pab.breaks[1] & mpab < pab.breaks[2], mpab.cat := "med"] %>% 
      .[mpab >= pab.breaks[2] & mpab <= pab.lims[2], mpab.cat := "high"]
    
    if(!is.na(option.flatten)){
      warning("'option.flatten = TRUE' only works with a single break value")
    #   temp <- pab.dt[!is.na(mpab.cat),.("qbreaks" = list(quantile(mpab, c(0,0.25,0.5,0.75,1)))),.("this.cat" = mpab.cat)]
    #   pab.dt[!is.na(mpab.cat) & mpab.cat == "high","bin.sub" := makebins(mpab,temp[this.cat == "high", qbreaks])] %>% 
    #     .[!is.na(mpab.cat) & mpab.cat == "low", "bin.sub" := makebins(mpab, tempb[this.cat == "low", qbreaks])]
    #   pab.dt[!is.na(mpab.cat) & !is.na(bin.sub), 
    #          mpab.cat := paste0(mpab.cat, "_", bin.sub), .(mpab.cat, bin.sub)]
    #   pab.dt[, bin.sub := NULL]
    }
    
    pab.breakdown <- pab.dt[!is.na(mpab.cat) & is.na(setx),.N,.(mpab.cat,this.ref)]    
    

  }else{
    pab.dt[mpab < pab.breaks[1] & mpab >= pab.lims[1], mpab.cat := "low"] %>%
      .[mpab >= pab.breaks[1] & mpab <= pab.lims[2], mpab.cat := "high"]
    
    if(!all(is.na(option.flatten))){
      
      if(length(option.flatten) == 1){
        temp <- data.table("qbreaks" = list(seq(pab.lims[1], pab.breaks, length.out = option.flatten + 1),
                                            seq(pab.breaks, pab.lims[2], length.out = option.flatten + 1)),
                           "this.cat" = c("low", "high"))
        pab.dt[, bin.sub := as.integer(NA)] %>% 
          .[!is.na(mpab.cat) & mpab.cat == "high", "bin.sub" := makebins(mpab, temp[this.cat == "high", unlist(qbreaks)])] %>% 
          .[!is.na(mpab.cat) & mpab.cat == "low", "bin.sub" := makebins(mpab, temp[this.cat == "low", unlist(qbreaks)])]
        pab.dt[!is.na(mpab.cat) & !is.na(bin.sub), 
               mpab.cat := paste0(mpab.cat, "_", bin.sub)] %>% 
          .[is.na(bin.sub), mpab.cat := as.character(NA)] %>% 
          .[, bin.sub := NULL]
      }else{
        pab.dt[, bin.sub := as.integer(NA)] %>% 
          .[!is.na(mpab.cat), "bin.sub" := makebins(mpab, option.flatten)]
        pab.dt[!is.na(mpab.cat) & !is.na(bin.sub), 
               mpab.cat := paste0(mpab.cat, "_", bin.sub)] %>% 
          .[is.na(bin.sub), mpab.cat := as.character(NA)] %>% 
          .[, bin.sub := NULL]
      }
    }
    
    pab.breakdown <- pab.dt[!is.na(mpab.cat) & is.na(setx),.N,.(mpab.cat,this.ref)]    

    
    
  }
  print(pab.breakdown)
  g1 <- ggplot(pab.dt[is.na(setx),]) + 
    geom_histogram(aes(x = mpab, fill = mpab.cat)) +
    labs(title = "all data")
  print(g1)
  if(option.gap != TRUE){
    if(option.round == "ten"){
      max.samplesize <- round(min(as.numeric(unlist(pab.breakdown$N)),na.rm = T)-5,-1)
    }else{
      max.samplesize <- round(min(as.numeric(unlist(pab.breakdown$N)),na.rm = T),1)
    }
  }else{
    pab.breakdown.gap <- copy(pab.breakdown)[mpab.cat != "med",]
    if(option.round == "ten"){
      max.samplesize <- round(min(as.numeric(unlist(pab.breakdown.gap$N)),na.rm = T)-5,-1)
    }else{
      max.samplesize <- round(min(as.numeric(unlist(pab.breakdown.gap$N)),na.rm = T),1)
    }
  }
  if(!is.na(train.size)){
    max.samplesize  <- train.size
  }
  
  
  print(paste0("Sample size per cat: ",max.samplesize))
  if(option.gap != TRUE){
    set.seed(seed = seedval)
    train.smiles <- unlist(pab.dt[is.na(setx) & !is.na(mpab.cat),
                                  .("xsmiles" = sample(smiles,size = max.samplesize, replace = FALSE)),
                                  .(mpab.cat,this.ref)][,xsmiles])# Use max.samplesize in training
    pab.dt[!is.na(mpab.cat), mpab.cat := unlist(str_split(mpab.cat, "_"))[1],.(mpab.cat)] %>% 
      .[,mpab.cat := factor(mpab.cat, levels = unique(mpab.cat))]
  }else{
    set.seed(seed = seedval)
    train.smiles <- unlist(pab.dt[is.na(setx) & !is.na(mpab.cat) & mpab.cat %in% c("low","high"),
                                  .("smiles" = sample(smiles,size = max.samplesize, replace = FALSE)),
                                  .(mpab.cat,this.ref)][,smiles])# Use max.samplesize in training
    pab.breaks <- mean(pab.breaks)
    pab.breakdown <- pab.dt[,.("low" = sum(mpab < pab.breaks[1]), 
                               "high" = sum(mpab >= pab.breaks[1])),.(this.ref)]
    pab.dt[mpab < pab.breaks[1] & mpab >= pab.lims[1], mpab.cat := "low"] %>%
      .[mpab >= pab.breaks[1] & mpab <= pab.lims[2], mpab.cat := "high"] %>% 
      .[,mpab.cat := factor(mpab.cat, levels = c("low", "high"))]
    
  }
  pab.dt[smiles%in%train.smiles,setx:="trainx"] %>%
    .[is.na(setx),setx:="testx"]
  g2 <- ggplot(pab.dt[setx == "trainx",]) +
    geom_histogram(aes(x = mpab, fill = mpab.cat)) +
    labs(title = "training set")
  print(g2)
  pab.dt[, this.ref := NULL]
  return(pab.dt)
}


make.predict <- function(desc.dt,this.model,smiles.eval=NA,smiles.col="smiles", type = "response"){
  setnames(desc.dt,smiles.col,"smiles")
  desc.names <- rownames(this.model$importance)
  if(all(is.na(smiles.eval))){smiles.eval <- desc.dt[,smiles]}
  desc.dt <- desc.dt[smiles %in% smiles.eval,]
  skiprow <- unique(which(is.na(desc.dt[,..desc.names])|is.inf(desc.dt[,..desc.names]),arr.ind=T)[,1])
  if(length(skiprow) > 0){
    pred.res <- predict(this.model,desc.dt[-skiprow,..desc.names], type = type)[match(smiles.eval,desc.dt[-skiprow,smiles])]
  }else{
    pred.res <- predict(this.model,desc.dt[,..desc.names], type = type)[match(smiles.eval,desc.dt[,smiles])]
  }
  return(pred.res)
}

make.predict.qrf <- function(desc.dt,this.model,smiles.eval=NA,smiles.col="smiles", what = c(0.05,0.5,0.95)){
  setnames(desc.dt,smiles.col,"smiles")
  desc.names <- rownames(this.model$importance)
  if(all(is.na(smiles.eval))){smiles.eval <- desc.dt[,smiles]}
  desc.dt <- desc.dt[smiles %in% smiles.eval,]
  skiprow <- unique(which(is.na(desc.dt[,..desc.names])|is.inf(desc.dt[,..desc.names]),arr.ind=T)[,1])
  if(length(skiprow) > 0){
    pred.res <- predict(this.model,desc.dt[-skiprow,..desc.names], what = what)[match(smiles.eval,desc.dt[-skiprow,smiles]),]
  }else{
    pred.res <- predict(this.model,desc.dt[,..desc.names], what = what)[match(smiles.eval,desc.dt[,smiles]),]
  }
  return(pred.res)
}

##### other utilities #####
makebins <- function(xvector,xlimbins){
  lx <- length(xlimbins)
  if(typeof(xvector) != "character" & typeof(xvector) != "list"){
    xvuse <- which(!is.nan(xvector) & 
                     !is.na(xvector) &
                     !is.inf(xvector) &
                     xvector <= last(xlimbins) &
                     xvector > xlimbins[1])
    resbin <- array(dim = c(length(xvector), 1))
    resbin[xvuse] <- unlist(lapply(xvector[xvuse], 
                                   function(x) which(xlimbins[1:(lx-1)] < x & xlimbins[2:lx] >= x)))
  }
  if(typeof(xvector) == "list"){
    xvuse <- which(!is.na(xvector) &
                     xvector %in% unlist(xlimbins))
    resbin <- array(dim = c(length(xvector), 1))
    resbin[xvuse] <- unlist(lapply(xvector[xvuse],
                                   function(x) which(unlist(lapply(xlimbins, 
                                                                   function(y) x %in% y)))))
  }
  if(typeof(xvector) == "character"){
    xvuse <- which(!is.na(xvector) &
                     xvector %in% xlimbins)
    resbin <- array(dim = c(length(xvector), 1))
    resbin[xvuse] <- unlist(lapply(xvector[xvuse],
                                   function(x) which(xlimbins %in% x)))
    
  }
  return(resbin)
}


is.inf <- function(x){
  
  if(is.data.table(x)|is.data.frame(x)){
    x <- as.matrix(x)
  }
  if(typeof(x)=="character"){
    print("Character type not allowed")
  }else{
    return(is.infinite(x))
  }
  
}

pred_int <- function(x0,y0,alpha = 0.95){
  x <- as.numeric(unlist(x0))[!is.na(x0) & !is.na(y0)]
  y <- as.numeric(unlist(y0))[!is.na(x0) & !is.na(y0)]
  if(length(x)!=length(y)){
    print("x and y must be equal length")
    break
  }
  
  xmin <- min(x)
  xmax <- max(x)
  xn <- seq(xmin,xmax,length.out = 1000)
  pi.mse <- sum((x-y)^2,na.rm=T)/(length(x)-2)
  tqt <- qt((1-alpha)/2,df=length(x)-2,lower.tail=FALSE)
  
  res.param <- list("MSE" = pi.mse,
                     "t12" = tqt,
                     "n" = length(x),
                     "xbar" = mean(x),
                     "denom" = sum((x-mean(x))^2)
                     )
  
  res.eqn <- paste0("yint = ",signif(res.param$t12,3), 
                    "*sqrt(", signif(pi.mse,3),
                    "(1 + 1/", length(x), 
                    " + (xn - ", signif(mean(x),3),
                    ")^2/", signif(res.param$denom,3),
                    "))"
                    )
  print(res.eqn)
  yn <-tqt*sqrt(pi.mse*(1+1/length(x)+(xn-mean(x))^2/sum((x-mean(x))^2)))
  
  res.dt <- data.table(xn,yn1=xn+yn,yn2=xn-yn, resint = yn)
  res.list <- list(res.dt = res.dt,
                   res.param = res.param,
                   res.eqn = res.eqn)
  return(res.list)
  
  
  
}
get_param <- function(param.name,Params,calling.func,default=NULL){
  if (is.null(Params[[param.name]]))
  {
    if (is.null(default)) stop(paste("Parameter",param.name,"needed in",calling.func))
    else return(default)
  } else return(Params[[param.name]])
}
tlsq <- function(x,y){
  xa <- x[!is.na(x) & !is.na(y)]
  ya <- y[!is.na(x) & !is.na(y)]
  v <- prcomp(data.table(x=xa,y=ya))$rotation
  beta <- v[2,1]/v[1,1]
  b0 <- mean(ya)-beta*mean(xa)
  xtrans<-beta*xa+b0
  yxresidual <- ya-(beta*xa+b0)
  coefficients <- c("intercept"=b0,"slope"=beta)
  res<- list("coefficients" = coefficients,"xtrans"=xtrans,"residuals"=yxresidual)
  return(res)
}

calc_Ffp <- function(species="Human",
                     chem.cas=NULL,
                     Kp.model = "schmitt",
                     default.to.human = TRUE,
                     restrictive.clearance=TRUE,
                     clu = NA){
  
  parameters <- parameterize_steadystate(
    chem.cas=chem.cas,
    species = species,
    default.to.human = default.to.human,
    suppress.messages = TRUE)
  
  
  if(is.na(clu)){
    cl <- calc_hepatic_clearance(
      chem.cas=chem.cas,
      species=species,
      default.to.human=default.to.human,
      hepatic.model = "unscaled",
      suppress.messages=TRUE)
  }else{
    cl <- copy(clu) # L/hr/kg
  }
  fup <- get_param("Funbound.plasma",parameters,"calc_Ffp") # unitless fraction
  if(!restrictive.clearance) cl <- cl/fup
  
  Qtotal.liverc <- get_param("Qtotal.liverc",parameters,"calc_Ffp") # L/h/kg^3/4BW
  Qtotal.liverc <- Qtotal.liverc / parameters[['BW']]^0.25 # L/h/kgBW
  
  rb2p <- available_rblood2plasma(
    chem.cas=chem.cas,
    species=species,
    suppress.messages=TRUE)
  #cl <- cl*parameters[['BW']] # convert to L/hr
  Ffp <- Qtotal.liverc/(Qtotal.liverc+cl*fup/rb2p)
  
  # if(Kp.model == "schmitt"){
  #   Kp.schmitt <- predict_partitioning_schmitt(chem.cas=chem.cas,species=species,default.to.human = default.to.human)
  # }
  
  return(Ffp)
}




xysummary <- function(ycalc=NA,ymeas=NA){
  ycalc1 <- ycalc[!is.na(ycalc) & !is.na(ymeas)]
  ymeas1 <- ymeas[!is.na(ycalc) & !is.na(ymeas)]
  ymeas.mean <- mean(ymeas1,na.rm=T)
  sst <- sum((ymeas1-ymeas.mean)^2,na.rm=T)
  sse <- sum((ymeas1-ycalc1)^2,na.rm=T)
  r2 <- 1-sse/sst
  xyrmse <- (mean((ymeas1-ycalc1)^2,na.rm=T))^0.5
  xycor <- cor(ycalc1,ymeas1)
  xyN <- length(ycalc1)
  fresult <- list("xysst"=sst,"xysse"=sse,"xyr2"=r2,"xyrmse"=xyrmse,"xycor" = xycor,"xyN" = xyN)
  return(fresult)
}

#fixnames
fixnames <- function(vnames){
  unlist(lapply(vnames,
                function(x) stri_replace_all_fixed(x,c("(",")","/","*"," ","%","?","-",".",":","="),c("","","per","x","_","percent","q","","","coln","eq"),vectorize_all=F)))
  
}

logit <- function(x=NA,delta=0.00001){
  y <- log((x+delta)/(1-x+delta))
  return(y)
}
sigmoid <- function(x=NA){
  y <- 1/(1+exp(-x))
  return(y)
}

na_rm <- function(x){
  if(!all(is.na(x))){
    y <- x[!is.na(x)]
  }else{
    y <- as.numeric(NA)
  }
  return(y)
}

