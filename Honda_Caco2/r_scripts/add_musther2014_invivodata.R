# Delete all objects from memory:
rm(list=ls())

packages <- c("readxl")
sapply(packages, require,character.only=TRUE) #Note, the "character.only" argument is necessary here

loc.wd <- "C:/Users/jwambaug/git/Comptox-Caco2/Honda_Caco2"

load(paste0(loc.wd, "/r_data/processed/all_gut_data_25MAR2019.RData")) 
gut.data <- subset(gut.data,!duplicated(gut.data))
dim(gut.data)

musther <- as.data.frame(read_excel(
                      paste0(loc.wd, "/r_data/lit_data/Musther-2014-DTXSID.xlsx"),
                      skip=1)) 
musther <- subset(musther, !is.na(DTXSID))
                      
dim(musther)
sum(musther$DTXSID %in% gut.data$dtxsid)
for (this.chem in unique(musther$DTXSID))
{
  this.row.musther <- musther$DTXSID==this.chem
  # Check if the chemical is already in gut.data:
  if (this.chem %in% gut.data$dtxsid)
  {
    this.row.gut.data <- which(gut.data$dtxsid==this.chem)
  # Otherwise add the new chemical to the table:
  } else {
    this.row.gut.data <- dim(gut.data)[1]+1
    gut.data[this.row.gut.data,"casrn"] <- musther[this.row.musther,"CASRN"]
    gut.data[this.row.gut.data,"dtxsid"] <- musther[this.row.musther,"DTXSID"]
  }
  # Check to see if we already have a value for human for this chemical:
    if (is.na(gut.data[this.row.gut.data, "vo_F"]) & # Varma in vivo Fbio
      is.na(gut.data[this.row.gut.data, "kim_fbioh"])) # Kim in vivo Fbio
    {
      gut.data[this.row.gut.data,"musther_Fbio_human"] <- 
        as.numeric(musther[this.row.musther,"F% Mean...45"])/100
      gut.data[this.row.gut.data,"musther_Fbio_human_ref"] <- 
        musther[this.row.musther,"References(human)"]
    }
    # We average across strains:
    gut.data[this.row.gut.data,"musther_Fbio_mouse"] <- 
      mean(as.numeric(musther[this.row.musther,"F% Mean...23"])/100, na.rm=TRUE)
    gut.data[this.row.gut.data,"musther_Fbio_mouse_ref"] <- 
      paste(musther[this.row.musther,"References(mouse)"], collapse=", ")
    gut.data[this.row.gut.data,"musther_Fbio_rat"] <- 
      mean(as.numeric(musther[this.row.musther,"F% Mean...28"])/100, na.rm=TRUE)
    gut.data[this.row.gut.data,"musther_Fbio_rat_ref"] <- 
      paste(musther[this.row.musther,"References(rat)"], collapse=", ")
    gut.data[this.row.gut.data,"musther_Fbio_dog"] <- 
      mean(as.numeric(musther[this.row.musther,"F% Mean...34"])/100, na.rm=TRUE)
    gut.data[this.row.gut.data,"musther_Fbio_dog_ref"] <- 
      paste(musther[this.row.musther,"References(dog)"], collapse=", ")
    gut.data[this.row.gut.data,"musther_Fbio_monkey"] <- 
      mean(as.numeric(musther[this.row.musther,"F% Mean...40"])/100, na.rm=TRUE)
    gut.data[this.row.gut.data,"musther_Fbio_monkey_ref"] <- 
      paste(musther[this.row.musther,"References(NHP)"], collapse=", ")    
}

dim(gut.data)

save(gut.data, file=paste0(loc.wd, "/r_data/processed/all_gut_data.RData")) 
