#### Collect literature data ####
oral_lit_collect <- function(){
  
  #### httk data ####
  chem.dt <- as.data.table(copy(chem.physical_and_invitro.data))
  setnames(chem.dt,c("CAS"),"casrn")
  chem.dt <- chem.dt[regexpr("Pearce",casrn)==-1 & regexpr("casrn",casrn)==-1,]
  chem.dt[SMILES.desalt=="FAIL",SMILES.desalt:=NA]
  chem.dt2 <- chem.dt[,.(preferred_name=Compound,Name=Compound,casrn=casrn,tk_fabh=Human.Fgutabs,tk_fabh_ref=Human.Fgutabs.Reference,fuph=Human.Funbound.plasma,clinth=Human.Clint,
                         tk_fabr=Rat.Fgutabs,tk_fabr_ref=Rat.Fgutabs.Reference,fupr=Rat.Funbound.plasma,clintr=Rat.Clint,logP=logP,MW=MW,
                         SMILES.desalt,SMILES.desalt.Reference)]
  
  

  
  
  
  #### Lanevskij - Pab ####
  lanv.pab <- fread(paste0(loc.wd, "/r_data/lit_data/lanevskij_2019_pab.csv"))
  lanv.ref <- fread(paste0(loc.wd, "/r_data/lit_data/lanevskij_2019_ref.csv"))
  lanv.pab <- lanv.pab[,c(1:14)]
  setnames(lanv.pab,c(5,7,8,10,14),
           c("log10_pab","stir_rpm","ref_num","r_angstrom","log10_pab_pred"))
  setnames(lanv.pab,names(lanv.pab),tolower(fixnames(names(lanv.pab))))
  setnames(lanv.ref, c(1:2),c("ref_num","ref_full"))
  lanv.pab2 <- lanv.pab[stir_rpm == "0",c(1:13)] %>% 
    .[ph >= 6.5 & ph <= 7.5,]

  #### Obringer - Pab and Fbior ####
  obringer.pab <- fread(paste0(loc.wd, "/r_data/lit_data/obringer2016_pab.csv"))
  obringer.fbior <- fread(paste0(loc.wd, "/r_data/lit_data/obringer2016_fbior.csv"))
  obringer.pab[,pab:=as.numeric(pab)] # NAs are non-detects for obringer pab
  obringer.fbior[, fbior_plasma := as.numeric(fbior_plasma)/100] %>%  # NAs are not determinable
    .[, fbior_urine := as.numeric(fbior_urine)/100]
  obringer.pab[casrn == "526-66-6", casrn := "525-66-6"] %>%  # correct casrn for propranolol
    .[casrn == "143-67-9", casrn := "865-21-4"] # vinblastine sulfate casrn corrected to vinblastine casrn
  
  ##### Chembl Pab; 16OCT2018 ####
  chembl1 <- fread(paste0(loc.wd,"/r_data/lit_data/chembl_caco2.csv"))
  setnames(chembl1,names(chembl1),tolower(fixnames(names(chembl1))))
  
  
  chembl2 <- copy(chembl1)
  chembl2[,gunits:=standard_units] %>% 
    .[,gval:=standard_value] %>%
    .[gunits=="",gval:=published_value] %>% 
    .[gunits=="",gunits:=published_units] %>% 
    .[,gunits:=gsub("'","^",gunits),.(gunits)] %>% 
    .[,gunits:=gsub("s-1","/s",gunits),.(gunits)] %>% 
    .[,gunits:=gsub(" ","",gunits),.(gunits)]
  
  chembl2[,gtype:=standard_type] %>% 
    .[,gtype:=gsub("Log","log",gtype),.(gtype)]
  
  chembl2[,use_val:=FALSE] %>% 
    .[regexpr("app",tolower(gtype))!=-1,use_val:=TRUE] %>% 
    .[regexpr("perm",tolower(gtype))!=-1,use_val:=TRUE] %>% 
    .[regexpr("pc",tolower(gtype))!=-1,use_val:=TRUE] %>% 
    .[regexpr("caco",tolower(gtype))!=-1,use_val:=TRUE] %>% 
    .[regexpr("pm",tolower(gtype))!=-1,use_val:=TRUE] 
  chembl3 <- chembl2[use_val==TRUE & 
                       !(gunits%in%c("","-","%")) & 
                       regexpr("%",gunits)==-1 &
                       (regexpr("min",gunits)!=-1 | regexpr("s",gunits)!=-1 | regexpr("hr",gunits)!=-1) &
                       !is.na(pubmed_id) &
                       gunits != "nM/s" &
                       relation %in% c("","=") &
                       pubmed_id %in% chembl2[,.N,.(pubmed_id)][N>10,pubmed_id],]
  chembl3[,gunits:=gsub("u","10^-6",gunits)] %>% 
    .[,gunits:=gsub("nm","10^-7cm",gunits)] %>% 
    .[gunits=="cm/s*10E6",gunits:="10^-6cm/s"] %>% 
    .[gunits=="10^6cm/s",gunits:="10^-6cm/s"] %>% 
    .[gunits=="10^7cm/s",gunits:="10^-7cm/s"] %>% 
    .[gunits=="10e-5*cm/min",gunits:="(10^-5)/60cm/s"] %>% 
    .[gunits=="10^-5cm/min",gunits:="(10^-5)/60cm/s"] %>% 
    .[gunits=="10e-6cm/s",gunits:="10^-6cm/s"]
  
  chembl3[gunits=="cm/s",conv_val:=1] %>% 
    .[gunits!="cm/s",conv_val:=eval(parse(text=gsub("cm/s","",gunits))),.(gunits)] %>% 
    .[,use_val:=NULL] %>% 
    .[,use_val:=conv_val*gval*10^6] %>% 
    .[,log10_val:=log10(use_val)]
  
  chembl3[,pubmed_id:=factor(pubmed_id,levels=unique(pubmed_id))]
  chembl4 <- chembl3[!(pubmed_id %in% unique(chembl3[log10_val > 4,pubmed_id])),]
  
  chembl_smiles <- get_babel(chembl4[,.(canonical_smiles)],smiles.col = "canonical_smiles",babel.args="-r")
  chembl4[,use_smiles:=unlist(chembl_smiles)]
  chembl4[1:4,.(canonical_smiles,use_smiles)]
  med.chembl <- chembl4[,.(medval=median(log10_val,na.rm=T),sdval=sd(log10_val,na.rm=T),.N),.(use_smiles)]
  chembl5 <- med.chembl[chembl4,on="use_smiles"]
  
  
  #### Ohagan et al. 2015 - Pab ####
  ohagan.pab <- fread(paste0(loc.wd,"/r_data/lit_data/ohagan_2015_drugcaco2.csv"))
  setnames(ohagan.pab,fixnames(names(ohagan.pab)))
  setnames(ohagan.pab,"Caco2_Papp_x_10^6_cmpers","pab_lit")
  ohagan.pab <- ohagan.pab[,pab_lit:=as.numeric(pab_lit)] %>% 
    .[!is.na(pab_lit)]
  ohagan.pab[,Name:=tolower(Name)]
  
  
  #### Dahlgren Lennernas et al. 2015 - Peff ####
  peff.dt <- fread(paste0(loc.wd,"/r_data/lit_data/dahlgren2015_peff.csv"))
  peff_loc_key <- peff.dt[1,1]
  setnames(peff.dt,as.character(unlist(peff.dt[2,])))
  peff.dt <- peff.dt[3:dim(peff.dt)[1],]
  setnames(peff.dt,c("Substance","Peff (x10-4 cm/s)"),c("Name","peff_lit"))
  setnames(peff.dt,fixnames(names(peff.dt)))
  peff.dt[,peff_lit:=as.numeric(peff_lit)]
  peff.dt <- peff.dt[!is.na(peff_lit),]
  setnames(peff.dt,c("Organ","Conc_mM","percent_coefficient_of_variation_CV","Reference","chemtype"),c("peff_organ","peff_conc","peff_CV","peff_ref","peff_chemtype"))
  peff.dt[,Name:=tolower(Name)]
  peff.dt2 <- peff.dt[peff_organ %in% c("P-SI","D-SI","Duodenum")]
  
  #### Varma, Obach et al. 2010 - Fah, Fgh, Ffph, Fbioh ####
  fglit.dt <- fread(paste0(loc.wd,"/r_data/lit_data/varma_obach_2010.csv"))
  setnames(fglit.dt,fixnames(names(fglit.dt)))
  setnames(fglit.dt,names(fglit.dt)[2:13],paste0("vo_",names(fglit.dt)[2:13]))
  fglit.dt[,Name:=tolower(Name)]
  fglit.dt <- fglit.dt[Name!=""]
  
  
  #### Kim et al. 2014 - Fbio ####
  kim.fbio <- fread(paste0(loc.wd, "/r_data/lit_data/Kim2014_Fbio.csv"))
  setnames(kim.fbio, names(kim.fbio),fixnames(names(kim.fbio)))
  kim.fbio[, fbioh := percentF]
  
  #### Wambaugh et al. 2018 rat-tk data ####
  load(paste0(loc.wd,"/r_data/tk_data/FitData-Dec2017.RData"))
  setnames(FitData,c("CAS","SelectedFgutabs","Fgutabs.pred"),c("casrn","fbior_meas","fbior_gplus"))
  FitData <- as.data.table(FitData)
  setnames(FitData, "casrn","fd_casrn")
  
  load(paste0(loc.wd,"/r_data/tk_data/PKstats-2018-01-16.RData"))
  pkstats2 <- as.data.table(PKstats)
  
  # dashboard.dt <- fread(paste0(loc.wd,"/r_data/chemprop/caco2lit_comptoxdashboard_06SEP2018v2.csv"))
  # setnames(dashboard.dt,tolower(fixnames(names(dashboard.dt))))
  # setnames(dashboard.dt,"input","Name")
  # dashboard.dt[Name=="prasugrel",
  #              c("casrn","smiles","preferred_name"):=list("150322-43-3",
  #                                                         "CC(=O)OC1=CC2=C(CCN(C2)C(C(=O)C2CC2)C2=CC=CC=C2F)S1",
  #                                                         "prasugel")]
  
  #### all our caco2 data ####
  load(paste0(loc.wd, "/r_data/caco2_qc/caco2dt_processed.RData"))

  
  
  #### collect chems ####
  gut.chems <- unique(rbind(chem.dt2[,.(ref.chnm = Name, ref.casrn = casrn, ref.smiles = SMILES.desalt, repref = "httk")],
                            unique(obringer.pab[,.(ref.chnm = chnm, ref.casrn = casrn, ref.smiles = NA, repref = "Obringer")]),
                            lanv.pab2[,.(ref.chnm = name, ref.casrn = NA, ref.smiles = smiles, repref = "Lanevskij")],
                            ohagan.pab[,.(ref.chnm = Name, ref.casrn = NA, ref.smiles = SMILES, repref = "Ohagan")],
                            peff.dt2[,.(ref.chnm = Name, ref.casrn = NA, ref.smiles = NA, repref = "Dahlgren")],
                            fglit.dt[,.(ref.chnm = Name, ref.casrn = NA, ref.smiles = NA, repref = "Varma&Obach")],
                            kim.fbio[,.(ref.chnm = Name, ref.casrn = NA, ref.smiles = Updated_SMILES, repref = "Kim")],
                            pkstats2[,.(ref.chnm = Compound, ref.casrn = CAS, ref.smiles = NA, repref = "Wambaugh")]
                            ))
  
  temp <- unique(caco2.dt[!(casrn %in% gut.chems$ref.casrn),.(casrn)])
  temp[,casrn:=paste0(" ",casrn),.(casrn)]
  #fwrite(temp, file = paste0(loc.wd,"/r_data/lit_data/dashboard_input_caco2_18MAR19.csv"))
  
  temp <- unique(gut.chems[repref == "Obringer",]) %>% 
    .[,casrn:=paste0(" ",ref.casrn),.(ref.casrn)]
  #fwrite(temp, file = paste0(loc.wd,"/r_data/lit_data/dashboard_input_caco2_25MAR19.csv"))

  gut.chems <- unique(rbind(gut.chems,
                      caco2.dt[,.(ref.chnm = compound, ref.casrn = casrn, ref.smiles = NA, repref = "honda_caco2")]))
  gut.chems[,X:= 1:length(ref.chnm)]
  dashboard.input.chnm <- unique(gut.chems[,.(ref.chnm)]) %>% 
    .[regexpr("Compound",ref.chnm)==-1,]
  dashboard.input.chnm[,adj.chnm := gsub("_", " ",ref.chnm)] %>% 
    .[, adj.chnm := tolower(adj.chnm)] %>% 
    .[, adj.chnm := gsub("anulline","aniline",adj.chnm)] %>% 
    .[, adj.chnm := gsub("anilline","aniline",adj.chnm)]
  dashboard.input.chnm <- unique(dashboard.input.chnm)
  # fwrite(unique(dashboard.input.chnm[,.(adj.chnm)]), file = paste0(loc.wd, "/r_data/lit_data/dashboard_input_chnm.csv"))
  
  dashboard.input.casrn <- copy(unique(gut.chems[!is.na(ref.casrn),.(ref.casrn)]))[,ref.casrn:=paste0(" ",ref.casrn),.(ref.casrn)]
  # fwrite(dashboard.input.casrn, file = paste0(loc.wd, "/r_data/lit_data/dashboard_input_casrn.csv"), quote = TRUE)
    
  dashboard.output.chnm <- fread(file = paste0(loc.wd, "/r_data/lit_data/dashboard_output_chnm_07MAR2019.csv"))
  setnames(dashboard.output.chnm,tolower(fixnames(names(dashboard.output.chnm))))

  
  dashboard.output.casrn <- fread(file = paste0(loc.wd, "/r_data/lit_data/dashboard_output_casrn_07MAR2019.csv"))
  setnames(dashboard.output.casrn,tolower(fixnames(names(dashboard.output.casrn))))
  
  dashboard.output.casrn2 <- fread(file = paste0(loc.wd, "/r_data/lit_data/dashboard_output_caco2_18MAR2019.csv"))
  setnames(dashboard.output.casrn2,tolower(fixnames(names(dashboard.output.casrn2))))
  
  dashboard.output.casrn3 <- fread(file = paste0(loc.wd, "/r_data/lit_data/dashboard_output_obringer_25MAR2019.csv"))
  setnames(dashboard.output.casrn3,tolower(fixnames(names(dashboard.output.casrn3))))
  
  dashboard.output.casrn <- rbind(dashboard.output.casrn, dashboard.output.casrn2, dashboard.output.casrn3, fill = TRUE)
  
  temp1 <- unique(dashboard.output.chnm[found_by != "NO_MATCH",.(adj.chnm = input,chnm_found_by = found_by)]) %>% 
    .[dashboard.input.chnm,on = "adj.chnm"] %>% 
    .[gut.chems, on = "ref.chnm"]
  gut.check <- unique(dashboard.output.casrn[found_by != "NO_MATCH",.(ref.casrn = input,casrn_found_by = found_by)]) %>% 
    .[temp1, on = "ref.casrn"]
  gut.check[!is.na(casrn_found_by), ctox_db := "use_casrn"] %>% 
    .[is.na(casrn_found_by) & !is.na(chnm_found_by), ctox_db := "use_chnm"]
  gut.check[,c("casrn_found_by","chnm_found_by") := NULL]
  temp1 <- unique(dashboard.output.casrn[found_by != "NO_MATCH",])[,ref.casrn := input][gut.check[ctox_db == "use_casrn",],on = "ref.casrn"]
  temp2 <- unique(dashboard.output.chnm[found_by != "NO_MATCH",])[,adj.chnm := input][gut.check[ctox_db == "use_chnm",],on = "adj.chnm"]

  gut.chems2 <- unique(rbind(temp1,temp2,gut.check[is.na(ctox_db),], fill = TRUE))
  fix.temp <- gut.chems2[X %in% c(3388,3594,3596),.(input,casrn,preferred_name,found_by,adj.chnm, ref.chnm, ref.casrn, ref.smiles, smiles, qsar_ready_smiles)]
  
  gut.chems2[is.na(casrn), casrn := paste0("NO_CAS_",adj.chnm),.(adj.chnm)] %>% 
    .[is.na(adj.chnm), adj.chnm := ref.chnm]

  gut.chems3 <- gut.chems2[!(adj.chnm == "ibutilide" & casrn == "122647-32-9") &
                       !(adj.chnm == "norfenfluramine" & casrn == "1886-26-6"),]
  
  gut.chems3[,use.smiles := smiles] %>% 
    .[!(is.na(use.smiles)), use.smiles.ref := "ctox_db"] %>% 
    .[is.na(use.smiles), use.smiles.ref := repref] %>% 
    .[is.na(use.smiles), use.smiles := ref.smiles]
  
  gut.chems3[use.smiles == "-", use.smiles.ref := repref] %>% 
    .[use.smiles == "-", use.smiles := ref.smiles] %>% 
    .[regexpr("[*]",use.smiles)!=-1,use.smiles:=NA] %>% 
    .[use.smiles %in% c("[Li]", "[K]"),use.smiles:=NA] %>% 
    .[use.smiles == "-",use.smiles := NA] %>% 
    .[qsar_ready_smiles == "-", qsar_ready_smiles := NA] %>% 
    .[ref.smiles == "-", ref.smiles := NA]
  
  gut.chems3[repref == "honda_caco2",.N,.(good_qsar_smi = !is.na(qsar_ready_smiles), good_smi = !is.na(use.smiles),repref)]
  
  gut.chems.dt <- gut.chems3[!is.na(qsar_ready_smiles),]
  gut.chems.all <- copy(gut.chems3)
  save(gut.chems.dt, file = paste0(loc.wd, "/r_data/chemprop/gut_chems.RData"))  
  save(gut.chems.all, file = paste0(loc.wd, "/r_data/chemprop/gut_chems_all.RData")) 

  #### save some data ####
  # save caco2 data
  
  lanv.pab3 <- copy(lanv.pab2)
  setnames(lanv.pab3,c("name","smiles"),c("ref.chnm","ref.smiles"))
  lanv.pab3 <- gut.chems.dt[repref == "Lanevskij", .(ref.chnm,ref.smiles,qsar_ready_smiles,casrn,dtxsid)][lanv.pab3, on = c("ref.chnm","ref.smiles")]
  
  obringer.pab2 <- copy(obringer.pab)
  setnames(obringer.pab2,c("casrn","chnm"),c("ref.casrn","ref.chnm"))
  obringer.pab2 <- unique(gut.chems.dt[repref == "Obringer", .(ref.casrn,ref.chnm,qsar_ready_smiles,casrn,dtxsid)])[obringer.pab2, on = c("ref.casrn","ref.chnm")]
  
  names(lanv.pab3)
  lit.caco2.dt <- rbind(lanv.pab3[!is.na(qsar_ready_smiles),
                                  .(dtxsid,casrn,qsar_ready_smiles,
                                    lit_pab = 10^log10_pab,
                                    lit_exp_pH = ph,
                                    lanv_refnum = ref_num,
                                    lanv_refset = set,
                                    repref = "Lanevskij"
                                  )],
                        obringer.pab2[!is.na(qsar_ready_smiles),
                                      .(dtxsid,casrn,qsar_ready_smiles,
                                        lit_pab = pab,
                                        lit_exp_pH = pH,
                                        obri_batch_id = batch_id,
                                        repref = "Obringer"
                                      )],
                        fill = TRUE)[,lit_pab_units := "10^-6 cm/s"]
  
  save(lit.caco2.dt, file = paste0(loc.wd, "/r_data/processed/lit_caco2_26MAR2019.RData"))
  
  caco2.dt2 <- copy(caco2.dt)
  setnames(caco2.dt2, c("casrn","dtxsid","compound"),c("ref.casrn","ref.dtxsid","ref.chnm"))
  caco2.dt2 <- gut.chems.all[repref == "honda_caco2", .(ref.chnm,ref.casrn,qsar_ready_smiles,casrn,dtxsid)][caco2.dt2, on = c("ref.chnm","ref.casrn")]
  
  save(lanv.pab3, file = paste0(loc.wd, "/r_data/processed/lanevskij_proc_25MAR2019.RData"))
  save(caco2.dt2, file = paste0(loc.wd, "/r_data/processed/our_caco2_25MAR2019.RData"))
  
  chem.dt3 <- copy(chem.dt2)
  setnames(chem.dt3,c("Name","SMILES.desalt","casrn"),c("ref.chnm","ref.smiles","ref.casrn"))
  chem.dt3 <- gut.chems.dt[repref == "httk", .(ref.chnm,ref.casrn,ref.smiles,qsar_ready_smiles,casrn,dtxsid)][chem.dt3, on = c("ref.chnm","ref.smiles","ref.casrn")]
  
  peff.dt3 <- copy(peff.dt2)
  setnames(peff.dt3,c("Name"),c("ref.chnm"))
  peff.dt3 <- gut.chems.dt[repref == "Dahlgren", .(ref.chnm,qsar_ready_smiles,casrn,dtxsid)][peff.dt3, on = c("ref.chnm")]
  meanpeff.dt <- peff.dt3[,.(mpeff_lit = mean(peff_lit,na.rm = T), peff_lit_unit = "10^-4 cm/s"),.(casrn,dtxsid,qsar_ready_smiles)]
  
  save(peff.dt3, file= paste0(loc.wd,"/r_data/processed/dahlgren_peff_26MAR2019.RData"))
  
  fglit.dt2 <- copy(fglit.dt)
  setnames(fglit.dt2,c("Name"),c("ref.chnm"))
  fglit.dt2 <- gut.chems.dt[repref == "Varma&Obach", .(ref.chnm,qsar_ready_smiles,casrn,dtxsid)][fglit.dt2, on = c("ref.chnm")]
  
  kim.fbio2 <- copy(kim.fbio)
  setnames(kim.fbio2,c("Name","Updated_SMILES"),c("ref.chnm","ref.smiles"))
  kim.fbio2 <- gut.chems.dt[repref == "Kim", .(ref.chnm,qsar_ready_smiles,casrn,dtxsid)][kim.fbio2, on = c("ref.chnm")]
  
  pkstats3 <- copy(pkstats2)
  setnames(pkstats3,c("Compound","CAS"),c("ref.chnm","ref.casrn"))
  pkstats3 <- gut.chems.dt[repref == "Wambaugh", .(ref.chnm,ref.casrn,qsar_ready_smiles,casrn,dtxsid)][pkstats3, on = c("ref.chnm","ref.casrn")]
  
  gut.data <- meanpeff.dt[!is.na(qsar_ready_smiles),
                          ][fglit.dt2[!is.na(qsar_ready_smiles),
                                      .(casrn, dtxsid, qsar_ready_smiles,
                                        vo_Vdss_lperkg, vo_fu,vo_CLt,vo_CLr,vo_F,
                                        vo_Fa,vo_Fg,vo_Fh, vo_Reference)
                                      ][kim.fbio2[!is.na(qsar_ready_smiles),
                                                  .(casrn, dtxsid,qsar_ready_smiles,
                                                    kim_fbioh = fbioh/100,
                                                    kim_ref = Source,
                                                    kim_smiles = ref.smiles)
                                                  ][unique(pkstats3[!is.na(qsar_ready_smiles),
                                                                    .(casrn,dtxsid, qsar_ready_smiles,
                                                                      pk_ref.chnm = ref.chnm,
                                                                      pk_ref.casrn = ref.casrn,
                                                                      pk_fbior = Fbio)
                                                                    ]
                                                  )[chem.dt3[!is.na(qsar_ready_smiles) &
                                                               (!is.na(clinth) | !is.na(clintr) | !is.na(fuph) | !is.na(fupr)),
                                                             .(casrn,dtxsid,qsar_ready_smiles,
                                                               tk_ref.chnm = ref.chnm,
                                                               tk_ref.casrn = ref.casrn,
                                                               tk_ref.smiles = ref.smiles,
                                                               tk_fabh,tk_fabh_ref,tk_fabr,tk_fabr_ref,
                                                               clinth,fuph,clintr,fupr)
                                                             ][unique(gut.chems.dt[,.(qsar_ready_smiles,dtxsid,casrn)]),
                                                               on = c("qsar_ready_smiles","dtxsid","casrn")],
                                                    on = c("qsar_ready_smiles","dtxsid","casrn")],
                                                  on = c("qsar_ready_smiles","dtxsid","casrn")],
                                        on = c("qsar_ready_smiles","dtxsid","casrn")],
                            on = c("qsar_ready_smiles","dtxsid","casrn")]
  
  save(gut.data, file = paste0(loc.wd,"/r_data/processed/all_gut_data_25MAR2019.RData"))
  
}

oral_desc_chemprop <- function(){

  #### run opera and padel ####
  #qsar_smiles.padel <- get_padel(smiles.table = qsar_smiles.dt, smiles.col = "smiles", name.col = "compound")  
  qsar_smiles.dt <- unique(gut.chems3[!is.na(qsar_ready_smiles),.(compound = paste0("casrn_",casrn), smiles = qsar_ready_smiles, average_mass,dtxsid),.(casrn)])
  qsar_smiles_obringer.dt <- unique(gut.chems3[!is.na(qsar_ready_smiles) & repref == "Obringer",.(compound = paste0("casrn_",casrn), smiles = qsar_ready_smiles, average_mass),.(casrn)])
  
  propv <-   c('BCF','logBCF','BP','logP','MP',
               'VP','logVP','WS', 'AOH', 'BioDeg', 'RB','ReadyBiodeg','HL','logHL','KM','logKM',
               'KOA','Koc','logKoc', 'RT', 'pKa', 'logD')
 
  setwd(paste0(loc.wd, "/r_data/opera_pred_25MAR2019"))
  #qsar_smiles_obringer.opera <- get_opera(smiles.table = qsar_smiles_obringer.dt, smiles.col = "smiles", name.col = "compound", clear.tempfiles = FALSE, opera.args = "-v 1", properties = propv)
  #save(qsar_smiles.opera, file = paste0(loc.wd, "/r_data/opera_pred/qsar_smiles_opera.RData"))
  load(file = paste0(loc.wd, "/r_data/opera_pred/qsar_smiles_opera.RData"))
  names(qsar_smiles.opera) %in% names(qsar_smiles_obringer.opera)
  qsar_smiles.opera <- unique(rbind(qsar_smiles.opera[MoleculeID != "casrn_865-21-4"],qsar_smiles_obringer.opera))
  qsar_smiles.opera <- qsar_smiles.dt[,.(dtxsid, qsar_ready_smiles = smiles, MoleculeID = compound, casrn)][qsar_smiles.opera, on = "MoleculeID"]
  
  mydtxsid <- unique(gut.chems3[!is.na(dtxsid),dtxsid])
  
  #### get measured values ####
  meas.prop <- query_dsstox.measchemprop(dtxsid.vector = mydtxsid,
                                         chemprop.vector = "theusual",
                                         password = "pass",
                                         user = "_dataminer")

  meas.prop <- as.data.table(meas.prop)  
  meas.prop.desc <- unique(meas.prop[,.(chemprop,chemprop_desc,chemprop_unit)])
  meas.prop.agg <- dcast(meas.prop[,.(dtxsid,chemprop,chemprop_val)], dtxsid ~ chemprop, fun.aggregate = list(mean,sd), value.var = "chemprop_val", drop = FALSE, fill = NA)
  innames <- names(meas.prop.agg)[!names(meas.prop.agg)%in%"dtxsid"]
  meannames <- innames[regexpr("mean",innames) != -1]
  sdnames <- innames[regexpr("sd", innames) != -1]
  setnames(meas.prop.agg,meannames,paste0(gsub("chemprop_val_mean_","",meannames),"_meas"))
  setnames(meas.prop.agg,sdnames, paste0(gsub("chemprop_val_sd_","",sdnames), "_meas_SD"))
  
  meas.prop.list <- list("meas.prop" = meas.prop,
                         "meas.prop.agg" = meas.prop.agg,
                         "meas.prop.desc" = meas.prop.desc)
  save(meas.prop.list, file = paste0(loc.wd, "/r_data/chemprop/meas_chemprop.RData"))
  
  names(qsar_smile.opera)
  which(!(meas.prop.agg$dtxsid %in% qsar_smiles.opera$dtxsid))
  
  chemprop.dt <- merge(meas.prop.agg, qsar_smiles.opera, on = "dtxsid", all = TRUE)
  measprop.nm <- gsub("chemprop_val_mean_","",meannames)
  predprop.nm <- gsub("_pred","",names(qsar_smiles.opera)[regexpr("_pred",names(qsar_smiles.opera))!=-1])
  ggplot(chemprop.dt)+geom_point(aes(x=LogP_meas,y=LogP_pred,color = as.logical(AD_LogP)),alpha = 0.3)
  ggplot(chemprop.dt)+geom_point(aes(x=pkaAcidicApparent_meas,y=pKa_a_pred,color = as.logical(AD_pKa)),alpha = 0.3)
  ggplot(chemprop.dt)+geom_point(aes(x=pkaBasicApparent_meas,y=pKa_b_pred,color = as.logical(AD_pKa)),alpha = 0.3)
  ggplot(chemprop.dt)+geom_point(aes(x=log10(HLC_meas),y=LogHL_pred,color = as.logical(AD_HL)),alpha = 0.3)

  sort(names(chemprop.dt))  

  chemprop.dt[,paste0(predprop.nm,"_use"):=lapply(paste0(predprop.nm,"_pred"),function(x) eval(parse(text = x)))] %>% 
    .[,paste0(predprop.nm,"_use_ref"):="OPERA"] %>% 
    .[is.na(qsar_ready_smiles),paste0(predprop.nm,"_use_ref"):="none"]
  
  warnings()
  sort(names(chemprop.dt)[regexpr("use",names(chemprop.dt))!=-1])
  unique(gsub("meas",names(meas.prop.agg)))
  measprop.nm
  list(1,2)
  
  chemprop.dt[(is.na(qsar_ready_smiles) | AD_BP == 0) & !is.na(BP_meas), c("BP_use","BP_use_ref"):=list(BP_meas,"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_MP == 0) & !is.na(MP_meas), c("MP_use","MP_use_ref"):=list(MP_meas,"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_KOA == 0) & !is.na(LogKOA_meas), c("LogKOA_use","LogKOA_use_ref"):=list(LogKOA_meas,"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_LogKoc == 0) & !is.na(LogKoc_meas), c("LogKoc_use","LogKoc_use_ref"):=list(LogKoc_meas,"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_HL == 0) & !is.na(HLC_meas), c("LogHL_use","LogHL_use_ref"):=list(log10(HLC_meas),"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_LogP == 0) & !is.na(LogP_meas), c("LogP_use","LogP_use_ref"):=list(LogP_meas,"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_VP == 0) & !is.na(VP_meas), c("LogVP_use","LogVP_use_ref"):=list(log10(VP_meas),"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_WS == 0) & !is.na(WSol_meas), c("LogWS_use","LogWS_use_ref"):=list(log10(WSol_meas),"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_pKa == 0) & !is.na(pkaAcidicApparent_meas), c("pKa_a_use","pKa_a_use_ref"):=list(pkaAcidicApparent_meas,"dsstox.measchemprop")] %>% 
    .[(is.na(qsar_ready_smiles) | AD_pKa == 0) & !is.na(pkaBasicApparent_meas), c("pKa_b_use","pKa_b_use_ref"):=list(pkaBasicApparent_meas,"dsstox.measchemprop")]

  chemprop.dt[,.N,.(AD_LogP,LogP_use_ref)]

  #save(chemprop.dt, file = paste0(loc.wd, "/r_data/chemprop/combined_chemprop.RData"))
  
  #save(qsar_smiles.opera, file = paste0(loc.wd,"/r_data/opera_pred/qsar_smiles_opera_25MAR19.RData"))
  #fwrite(qsar_smiles.dt[MW < 1000,.(smiles,compound)], file = paste0(loc.wd,"/r_data/opera_pred/qsar_smiles_smallmw.smi"), col.names = FALSE, sep = " ")
  #fwrite(qsar_smiles.dt[MW >= 1000,.(smiles,compound)], file = paste0(loc.wd,"/r_data/opera_pred/qsar_smiles_largemw.smi"), col.names = FALSE, sep = " ")
  
  #### get cdk and mordred descriptors ####
  cdk.desc <- rbind(fread(paste0(loc.wd,"/r_data/opera_pred/qsar_smiles_smallmw_output.csv")),
                    fread(paste0(loc.wd,"/r_data/opera_pred/qsar_smiles_largemw_output.csv"))
                    )
  padel.desc <- fread(paste0(loc.wd,"/r_data/opera_pred/PadelDesc.csv"))
  padel.fp <- fread(paste0(loc.wd,"/r_data/opera_pred/PadelFP.csv"))
  load(paste0(loc.wd,"/r_data/opera_pred/qsar_smiles_opera.RData"))
  opera.props <- copy(qsar_smiles.opera)
  
  names(cdk.desc)[1]
  names(padel.fp)[1]
  names(padel.desc)[1]
  names(opera.props)[1]
  names(cdk.desc)
  
  mordred.desc <- fread(paste0(loc.wd, "/python/mordred_result_no3d_26MAR2019.csv"))
  
  mordred.desc[750,TopoPSA]  
  padel.desc[750,TopoPSA]
  cdk.desc[750,TopoPSA]
  cdk.desc[751,1]
  mordred.desc[750,TopoPSA]
  
  cdk.desc[749,MW]
  mordred.desc[750,chem_id]
  names(cdk.desc)[1]
  cdk.desc[Title == "casrn_4342-03-4",TopoPSA]
  padel.desc[750,1]
  
  mordred.desc[,V1:= NULL] %>% 
    .[,casrn := gsub("casrn_","",chem_id), .(chem_id)] %>% 
    .[,chem_id := NULL]
  setnames(mordred.desc, "smiles","qsar_ready_smiles")
  mordred.desc <- qsar_smiles.dt[,.(casrn,dtxsid)][mordred.desc, on= c("casrn")]
  names(qsar_smiles.dt)
  dim(gut.chems.red)
  rm(mordred.desc)
  mordred.desc <- copy(mordred.desc2)
  rm(mordred.desc2)
  
  save(mordred.desc, file = paste0(loc.wd, "/r_data/descriptors/mordred_desc_26MAR2019.RData"))
  gut.mordred <- mordred.desc[gut.chems.red, on = c("qsar_ready_smiles", "casrn")]
  names(gut.chems.red)
?fread
  test.smiles1 <- fread(paste0(loc.wd, "/r_data/opera_pred/test.smi"), header = FALSE, col.names = c("qsar_ready_smiles","chem_id"))
  padel.desc1 <- fread(paste0(loc.wd, "/r_data/opera_pred/PadelDesc.csv"))
  padel.desc1 <- test.smi[,.("Name" = chem_id,smiles)][padel.desc1, on = "smiles"]
  
  test.smiles2 <- fread(paste0(loc.wd, "/r_data/opera_pred_25MAR2019/test.smi"), header = FALSE, col.names = c("qsar_ready_smiles","chem_id"))
  padel.desc2 <- fread(paste0(loc.wd, "/r_data/opera_pred_25MAR2019/PadelDesc.csv"))
  
  test.smiles <- unique(rbind(test.smiles1, test.smiles2))
  padel.desc <- unique(rbind(padel.desc1,padel.desc2))

  setnames(padel.desc, "Name","chem_id")
  
  padel.desc <- test.smiles[padel.desc, on = "chem_id"]
  padel.desc <- qsar_smiles.opera[,.(dtxsid,chem_id = MoleculeID, casrn, qsar_ready_smiles)][padel.desc, on=c("qsar_ready_smiles","chem_id")]
  
  save(padel.desc, file = paste0(loc.wd, "/r_data/descriptors/padel_desc_26MAR2019.RData"))
"TopoPSA" %in% names(mordred.desc)
  comp.dt <- padel.desc[,.(dtxsid, tpsa_padel = TopoPSA)][mordred.desc[,.(dtxsid,tpsa_mordred = TopoPSA)], on = "dtxsid"]
ggplot(comp.dt)+
  geom_point(aes(x=tpsa_padel, y=tpsa_mordred))

dt2 <- data.table(padel.nm = sort(names(padel.desc)),
           mordred.nm = sort(names(mordred.desc)))
      
}



{
  # load our processed caco2 measurements
  load(paste0(loc.wd,"/r_data/caco2_qc/",caco2_qc_file,".RData"))
  caco2.dt[,Pab:=Pab/rec_ab]
  
  # additional chem identifiers
  extra_chemid <- fread(paste0(loc.wd,"/r_data/chemprop/peff_missing_chemid.csv")) # pubchem search
  
  # processing
  rxcaco2.dt <- dashboard.dt[,.(Name,preferred_name,casrn)][drugcaco2.dt[,.(Name,pab_lit,rxpab_refn = Reference_number,rxpabref = Endnote)],on="Name"]
  fglit.dt2 <- dashboard.dt[,.(Name,preferred_name,casrn)][fglit.dt,on="Name"][casrn=="-",casrn:=vo_CAS]
  fglit.dt2 <- smiles.dt[,.(smiles,casrn)][fglit.dt2,on="casrn"]
  peff.dt3 <- dashboard.dt[,.(Name,preferred_name,casrn)][peff.dt2,on="Name"]
  peff.dt3 <- extra_chemid[,.(Name,preferred_name,casrn)][peff.dt3,on="Name"] %>% 
    .[is.na(casrn)|casrn=="-",casrn := i.casrn] %>% 
    .[,i.casrn:=NULL] %>% 
    .[is.na(preferred_name)|preferred_name=="-",preferred_name:=i.preferred_name] %>% 
    .[,i.preferred_name:=NULL]
  peff.dt4 <- smiles.dt[,.(casrn,smiles,babel_smiles)][peff.dt3,on="casrn"]
  peff.dt4[is.na(babel_smiles),]
  save(peff.dt4,file=paste0(loc.wd,"/r_data/processed/processed_peff_14FEB2019.RData"))
  peffl.babel <- peff.dt4[,.(peffl= median(log10(peff_lit),na.rm=T)),.(smiles=babel_smiles)]
  peff.dt3 <- smiles.dt[,.(casrn,smiles)][peff.dt3,on="casrn"]
  peff.med.smiles <- peff.dt3[,.(peffl=median(log10(peff_lit),na.rm=T),sdpeffl=sd(log10(peff_lit),na.rm=T)),.(smiles)]
  peff.med.casrn <- peff.dt3[,.(peffl=median(log10(peff_lit),na.rm=T),sdpeffl=sd(log10(peff_lit),na.rm=T)),.(casrn)]
  
  rxcaco2.dt <- smiles.dt[,.(casrn,smiles)][rxcaco2.dt,on="casrn"]
  caco2.dt2 <- smiles.dt[,.(casrn,smiles)][caco2.dt,on="casrn"]
  
  # combine our data with pharma data
  caco2.medref.smiles <- rbind(rxcaco2.dt[rxpab_refn!=10,.(mpab=median(log10(pab_lit),na.rm=T),msdpab=sd(log10(pab_lit),na.rm = T)),.(smiles,pabrefnum=rxpab_refn,pabref=rxpabref)],
                               caco2.dt2[,.(mpab=median(log10(Pab),na.rm=T),msdpab=sd(log10(Pab),na.rm=T),pabrefnum=100,pabref="Honda"),.(smiles)]) %>% 
    .[!is.na(mpab),]
  
  caco2.medref.casrn <- rbind(rxcaco2.dt[rxpab_refn!=10,.(mpab=median(log10(pab_lit),na.rm=T),msdpab=sd(log10(pab_lit),na.rm = T)),.(casrn,pabrefnum=rxpab_refn,pabref=rxpabref)],
                              caco2.dt2[,.(mpab=median(log10(Pab),na.rm=T),msdpab=sd(log10(Pab),na.rm=T),pabrefnum=100,pabref="Honda"),.(casrn)]) %>% 
    .[!is.na(mpab),]
  
  
  ##### calibrating measured pab to peff by lab - casrn #####
  pabpeffref <- peff.med.casrn[caco2.medref.casrn,on="casrn"] %>% 
    .[!is.na(peffl) & !is.na(mpab),]
  caco2.medref.casrn[,.(lcas=length(unique(casrn))),.(pabrefnum,pabref)]
  pabpeffref[,.(lcas=length(unique(casrn))),.(pabrefnum,pabref)]
  pabpeffref[,lcas:=length(unique(casrn)),.(pabrefnum)]
  pabpeffref <- pabpeffref[lcas>=5,]
  
  pabpeffcal <- pabpeffref[,as.list(summary(lm(peffl~mpab))$coefficients[1:2,1]),.(pabref)]
  setnames(pabpeffcal,names(pabpeffcal),c("pabref","intercept","slope"))
  pabpefftls <- pabpeffref[,as.list(tlsq(x=mpab,y=peffl)$coefficients[1:2]),.(pabref)]
  
  caco2.mpefcal <- caco2.medref.casrn[pabpeffcal,on="pabref"]
  caco2.mpefcal[,pefcal:=mpab*slope+intercept]
  caco2.mpefcal <- pabpefftls[,.(tls_int= intercept,tls_slope=slope,pabref)][caco2.mpefcal,on="pabref"]
  caco2.mpefcal[,peftls:=mpab*tls_slope+tls_int]
  
  caco2.pefclean <- caco2.mpefcal[,.(mpefcal=median(pefcal,na.rm=T),mpeftls=median(peftls,na.rm=T),mmpab=median(mpab,na.rm=T)),.(casrn)]
  caco2.pefclean.casrn <- smiles.dt[,.(casrn,smiles)][caco2.pefclean,on="casrn"]
  caco2.pefclean.ours <- smiles.dt[,.(casrn,smiles)][caco2.mpefcal[pabref == "Honda", .(mpefcal=median(pefcal,na.rm=T),mpeftls=median(peftls,na.rm=T),mmpab=median(mpab,na.rm=T)),.(casrn)]]
  caco2.pefclean.smiles <- caco2.pefclean.casrn[,.(mpefcal=median(mpefcal,na.rm=T),mpeftls=median(mpeftls,na.rm=T),mmpab=median(mmpab,na.rm=T)),.(smiles)]
  caco2.pefclean2 <- peff.med.smiles[caco2.pefclean.smiles,on="smiles"]
  save(caco2.pefclean.smiles,file=paste0(loc.wd,"/r_data/processed/caco2_pefclean_smiles_14FEB2019.RData"))
  save(caco2.pefclean.casrn,file=paste0(loc.wd,"/r_data/caco2_pefclean_casrn_14FEB2019.RData"))
  
  setnames(chem.dt3,c("smiles","tk_casrn"),c("tk_smiles","casrn"))
  
  # this is our measured pab only 14FEB2019 GSH
  allgut.casrn <- merge(merge(merge(merge(
    caco2.pefclean.ours,fglit.dt2,by="casrn",all=TRUE),
    peff.med.casrn,by="casrn",all=TRUE),
    FitData[,.(casrn = fd_casrn,fd_smiles=smiles,fbior_meas)],by="casrn",all=TRUE),
    chem.dt3,by="casrn",all=TRUE)
  
  allgut.casrn <- allgut.casrn[!is.na(vo_F)|!is.na(vo_Fa)|!is.na(vo_Fg)|!is.na(vo_Fh)|!is.na(tk_fabh)|!is.na(tk_fabr)|!is.na(fbior_meas)|!is.na(peffl)|!is.na(mmpab),]
  
  allgut.casrn <- smiles.dt[,.(new_smiles = smiles,casrn)][allgut.casrn, on="casrn"] %>% 
    .[,smiles:=new_smiles] %>% 
    .[,new_smiles:=NULL]
  
  all(caco2.pefclean$smiles %in% allgut.casrn$smiles)
  allgut.casrn[,c("smiles.x","smiles.y","preferred_name.x","preferred_name.y","smiles","Name.x","Name.y"):=NULL]
  allgut.casrn <- unique(smiles.dt[,.(preferred_name=preferred_name[1]),.(casrn,smiles)])[allgut.casrn,on="casrn"]
  save(allgut.casrn,file=paste0(loc.wd,"/r_data/processed/allgut_casrn_14FEB2019.RData"))
  save(rxcaco2.dt, file = paste0(loc.wd,"/r_data/processed/rxcaco2_14FEB2019.RData"))
}
