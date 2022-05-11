# Package ID: edi.577.2 Cataloging System:https://pasta.edirepository.org.
# Data set title: The Sierra Lakes Inventory Project: Non-Native fish and community composition of lakes and ponds in the Sierra Nevada, California.
# Data set creator:  Roland Knapp - University of California - Santa Barbara, Sierra Nevada Aquatic Research Laboratory 
# Data set creator:  Claire Pavelka - Environmental Data Initiative 
# Data set creator:  Ericka Hegeman - University of California - Santa Barbara, Sierra Nevada Aquatic Research Laboratory 
# Data set creator:  Thomas Smith - University of California - Santa Barbara, Sierra Nevada Aquatic Research Laboratory 
# Contact:  Roland Knapp - Research biologist University of California - Santa Barbara, Sierra Nevada Aquatic Research Laboratory  - roland.knapp@ucsb.edu
# Contact:  Thomas Smith - Research biologist University of California - Santa Barbara, Sierra Nevada Aquatic Research Laboratory  - tcsmith@ucsb.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/436ae7becbd2daaea8575ba0399ef669" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "survey_date",     
                 "amphibian_species",     
                 "amphibian_life_stage",     
                 "amphibian_state",     
                 "amphibian_number",     
                 "amphibian_location",     
                 "amphibian_voucher"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$lake_id)=="factor") dt1$lake_id <-as.numeric(levels(dt1$lake_id))[as.integer(dt1$lake_id) ]               
if (class(dt1$lake_id)=="character") dt1$lake_id <-as.numeric(dt1$lake_id)                                   
# attempting to convert dt1$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-" %Y-%m-%d"
tmp1survey_date<-as.Date(dt1$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1survey_date) == length(tmp1survey_date[!is.na(tmp1survey_date)])){dt1$survey_date <- tmp1survey_date } else {print("Date conversion failed for dt1$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1survey_date) 
if (class(dt1$amphibian_species)!="factor") dt1$amphibian_species<- as.factor(dt1$amphibian_species)
if (class(dt1$amphibian_life_stage)!="factor") dt1$amphibian_life_stage<- as.factor(dt1$amphibian_life_stage)
if (class(dt1$amphibian_state)!="factor") dt1$amphibian_state<- as.factor(dt1$amphibian_state)
if (class(dt1$amphibian_number)=="factor") dt1$amphibian_number <-as.numeric(levels(dt1$amphibian_number))[as.integer(dt1$amphibian_number) ]               
if (class(dt1$amphibian_number)=="character") dt1$amphibian_number <-as.numeric(dt1$amphibian_number)
if (class(dt1$amphibian_location)!="factor") dt1$amphibian_location<- as.factor(dt1$amphibian_location)
if (class(dt1$amphibian_voucher)=="factor") dt1$amphibian_voucher <-as.numeric(levels(dt1$amphibian_voucher))[as.integer(dt1$amphibian_voucher) ]               
if (class(dt1$amphibian_voucher)=="character") dt1$amphibian_voucher <-as.numeric(dt1$amphibian_voucher)

# Convert Missing Values to NA for non-dates

dt1$amphibian_number <- ifelse((trimws(as.character(dt1$amphibian_number))==trimws("NA")),NA,dt1$amphibian_number)               
suppressWarnings(dt1$amphibian_number <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$amphibian_number))==as.character(as.numeric("NA"))),NA,dt1$amphibian_number))
dt1$amphibian_voucher <- ifelse((trimws(as.character(dt1$amphibian_voucher))==trimws("NA")),NA,dt1$amphibian_voucher)               
suppressWarnings(dt1$amphibian_voucher <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$amphibian_voucher))==as.character(as.numeric("NA"))),NA,dt1$amphibian_voucher))


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(amphibian_species)
summary(amphibian_life_stage)
summary(amphibian_state)
summary(amphibian_number)
summary(amphibian_location)
summary(amphibian_voucher) 
# Get more details on character variables

summary(as.factor(dt1$amphibian_species)) 
summary(as.factor(dt1$amphibian_life_stage)) 
summary(as.factor(dt1$amphibian_state)) 
summary(as.factor(dt1$amphibian_location))
detach(dt1)               


inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/8a735c05ed737d8a32833c56fe37356a" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "collect_date",     
                 "id_date",     
                 "identifier",     
                 "shrimp_species",     
                 "genetics"    ), check.names=TRUE)

unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt2$lake_id)=="factor") dt2$lake_id <-as.numeric(levels(dt2$lake_id))[as.integer(dt2$lake_id) ]               
if (class(dt2$lake_id)=="character") dt2$lake_id <-as.numeric(dt2$lake_id)                                   
# attempting to convert dt2$collect_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp2collect_date<-as.Date(dt2$collect_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp2collect_date) == length(tmp2collect_date[!is.na(tmp2collect_date)])){dt2$collect_date <- tmp2collect_date } else {print("Date conversion failed for dt2$collect_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp2collect_date)                                    
# attempting to convert dt2$id_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp2id_date<-as.Date(dt2$id_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp2id_date) == length(tmp2id_date[!is.na(tmp2id_date)])){dt2$id_date <- tmp2id_date } else {print("Date conversion failed for dt2$id_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp2id_date) 
if (class(dt2$identifier)=="factor") dt2$identifier <-as.numeric(levels(dt2$identifier))[as.integer(dt2$identifier) ]               
if (class(dt2$identifier)=="character") dt2$identifier <-as.numeric(dt2$identifier)
if (class(dt2$shrimp_species)!="factor") dt2$shrimp_species<- as.factor(dt2$shrimp_species)
if (class(dt2$genetics)!="factor") dt2$genetics<- as.factor(dt2$genetics)

# Convert Missing Values to NA for non-dates

dt2$genetics <- as.factor(ifelse((trimws(as.character(dt2$genetics))==trimws("NA")),NA,as.character(dt2$genetics)))


# Here is the structure of the input data frame:
str(dt2)                            
attach(dt2)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(collect_date)
summary(id_date)
summary(identifier)
summary(shrimp_species)
summary(genetics) 
# Get more details on character variables

summary(as.factor(dt2$shrimp_species)) 
summary(as.factor(dt2$genetics))
detach(dt2)               


inUrl3  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/071642fa72ba780ee90ed36350d82745" 
infile3 <- tempfile()
try(download.file(inUrl3,infile3,method="curl"))
if (is.na(file.size(infile3))) download.file(inUrl3,infile3,method="auto")


dt3 <-read.csv(infile3,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "survey_date",     
                 "fish_id",     
                 "fish_species",     
                 "fish_length",     
                 "fish_weight",     
                 "fish_sex",     
                 "fish_egg_stage",     
                 "fish_otolith",     
                 "fish_age"    ), check.names=TRUE)

unlink(infile3)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt3$lake_id)=="factor") dt3$lake_id <-as.numeric(levels(dt3$lake_id))[as.integer(dt3$lake_id) ]               
if (class(dt3$lake_id)=="character") dt3$lake_id <-as.numeric(dt3$lake_id)                                   
# attempting to convert dt3$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp3survey_date<-as.Date(dt3$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp3survey_date) == length(tmp3survey_date[!is.na(tmp3survey_date)])){dt3$survey_date <- tmp3survey_date } else {print("Date conversion failed for dt3$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp3survey_date) 
if (class(dt3$fish_id)=="factor") dt3$fish_id <-as.numeric(levels(dt3$fish_id))[as.integer(dt3$fish_id) ]               
if (class(dt3$fish_id)=="character") dt3$fish_id <-as.numeric(dt3$fish_id)
if (class(dt3$fish_species)!="factor") dt3$fish_species<- as.factor(dt3$fish_species)
if (class(dt3$fish_length)=="factor") dt3$fish_length <-as.numeric(levels(dt3$fish_length))[as.integer(dt3$fish_length) ]               
if (class(dt3$fish_length)=="character") dt3$fish_length <-as.numeric(dt3$fish_length)
if (class(dt3$fish_weight)=="factor") dt3$fish_weight <-as.numeric(levels(dt3$fish_weight))[as.integer(dt3$fish_weight) ]               
if (class(dt3$fish_weight)=="character") dt3$fish_weight <-as.numeric(dt3$fish_weight)
if (class(dt3$fish_sex)!="factor") dt3$fish_sex<- as.factor(dt3$fish_sex)
if (class(dt3$fish_egg_stage)!="factor") dt3$fish_egg_stage<- as.factor(dt3$fish_egg_stage)
if (class(dt3$fish_otolith)!="factor") dt3$fish_otolith<- as.factor(dt3$fish_otolith)
if (class(dt3$fish_age)=="factor") dt3$fish_age <-as.numeric(levels(dt3$fish_age))[as.integer(dt3$fish_age) ]               
if (class(dt3$fish_age)=="character") dt3$fish_age <-as.numeric(dt3$fish_age)

# Convert Missing Values to NA for non-dates

dt3$fish_length <- ifelse((trimws(as.character(dt3$fish_length))==trimws("NA")),NA,dt3$fish_length)               
suppressWarnings(dt3$fish_length <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt3$fish_length))==as.character(as.numeric("NA"))),NA,dt3$fish_length))
dt3$fish_weight <- ifelse((trimws(as.character(dt3$fish_weight))==trimws("NA")),NA,dt3$fish_weight)               
suppressWarnings(dt3$fish_weight <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt3$fish_weight))==as.character(as.numeric("NA"))),NA,dt3$fish_weight))
dt3$fish_sex <- as.factor(ifelse((trimws(as.character(dt3$fish_sex))==trimws("NA")),NA,as.character(dt3$fish_sex)))
dt3$fish_egg_stage <- as.factor(ifelse((trimws(as.character(dt3$fish_egg_stage))==trimws("NA")),NA,as.character(dt3$fish_egg_stage)))
dt3$fish_age <- ifelse((trimws(as.character(dt3$fish_age))==trimws("NA")),NA,dt3$fish_age)               
suppressWarnings(dt3$fish_age <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt3$fish_age))==as.character(as.numeric("NA"))),NA,dt3$fish_age))


# Here is the structure of the input data frame:
str(dt3)                            
attach(dt3)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(fish_id)
summary(fish_species)
summary(fish_length)
summary(fish_weight)
summary(fish_sex)
summary(fish_egg_stage)
summary(fish_otolith)
summary(fish_age) 
# Get more details on character variables

summary(as.factor(dt3$fish_species)) 
summary(as.factor(dt3$fish_sex)) 
summary(as.factor(dt3$fish_egg_stage)) 
summary(as.factor(dt3$fish_otolith))
detach(dt3)               


inUrl4  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/bc03219e367d3a4b7b1758c8af598174" 
infile4 <- tempfile()
try(download.file(inUrl4,infile4,method="curl"))
if (is.na(file.size(infile4))) download.file(inUrl4,infile4,method="auto")


dt4 <-read.csv(infile4,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "lake_area_nbr",     
                 "lake_perimeter_nbr",     
                 "lake_type_code",     
                 "lake_drainage_name",     
                 "lake_elevation_nbr",     
                 "lake_quad_name",     
                 "lake_county_name",     
                 "lake_juris_name",     
                 "lake_wilderness_name"    ), check.names=TRUE)

unlink(infile4)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt4$lake_id)=="factor") dt4$lake_id <-as.numeric(levels(dt4$lake_id))[as.integer(dt4$lake_id) ]               
if (class(dt4$lake_id)=="character") dt4$lake_id <-as.numeric(dt4$lake_id)
if (class(dt4$lake_area_nbr)=="factor") dt4$lake_area_nbr <-as.numeric(levels(dt4$lake_area_nbr))[as.integer(dt4$lake_area_nbr) ]               
if (class(dt4$lake_area_nbr)=="character") dt4$lake_area_nbr <-as.numeric(dt4$lake_area_nbr)
if (class(dt4$lake_perimeter_nbr)=="factor") dt4$lake_perimeter_nbr <-as.numeric(levels(dt4$lake_perimeter_nbr))[as.integer(dt4$lake_perimeter_nbr) ]               
if (class(dt4$lake_perimeter_nbr)=="character") dt4$lake_perimeter_nbr <-as.numeric(dt4$lake_perimeter_nbr)
if (class(dt4$lake_type_code)!="factor") dt4$lake_type_code<- as.factor(dt4$lake_type_code)
if (class(dt4$lake_drainage_name)!="factor") dt4$lake_drainage_name<- as.factor(dt4$lake_drainage_name)
if (class(dt4$lake_elevation_nbr)=="factor") dt4$lake_elevation_nbr <-as.numeric(levels(dt4$lake_elevation_nbr))[as.integer(dt4$lake_elevation_nbr) ]               
if (class(dt4$lake_elevation_nbr)=="character") dt4$lake_elevation_nbr <-as.numeric(dt4$lake_elevation_nbr)
if (class(dt4$lake_quad_name)!="factor") dt4$lake_quad_name<- as.factor(dt4$lake_quad_name)
if (class(dt4$lake_county_name)!="factor") dt4$lake_county_name<- as.factor(dt4$lake_county_name)
if (class(dt4$lake_juris_name)!="factor") dt4$lake_juris_name<- as.factor(dt4$lake_juris_name)
if (class(dt4$lake_wilderness_name)!="factor") dt4$lake_wilderness_name<- as.factor(dt4$lake_wilderness_name)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt4)                            
attach(dt4)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(lake_area_nbr)
summary(lake_perimeter_nbr)
summary(lake_type_code)
summary(lake_drainage_name)
summary(lake_elevation_nbr)
summary(lake_quad_name)
summary(lake_county_name)
summary(lake_juris_name)
summary(lake_wilderness_name) 
# Get more details on character variables

summary(as.factor(dt4$lake_type_code)) 
summary(as.factor(dt4$lake_drainage_name)) 
summary(as.factor(dt4$lake_quad_name)) 
summary(as.factor(dt4$lake_county_name)) 
summary(as.factor(dt4$lake_juris_name)) 
summary(as.factor(dt4$lake_wilderness_name))
detach(dt4)               


inUrl5  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/69e6f7ecef2027c21c9261bd8a866bc9" 
infile5 <- tempfile()
try(download.file(inUrl5,infile5,method="curl"))
if (is.na(file.size(infile5))) download.file(inUrl5,infile5,method="auto")


dt5 <-read.csv(infile5,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "survey_date",     
                 "littoral_type",     
                 "littoral_amount"    ), check.names=TRUE)

unlink(infile5)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt5$lake_id)=="factor") dt5$lake_id <-as.numeric(levels(dt5$lake_id))[as.integer(dt5$lake_id) ]               
if (class(dt5$lake_id)=="character") dt5$lake_id <-as.numeric(dt5$lake_id)                                   
# attempting to convert dt5$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp5survey_date<-as.Date(dt5$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp5survey_date) == length(tmp5survey_date[!is.na(tmp5survey_date)])){dt5$survey_date <- tmp5survey_date } else {print("Date conversion failed for dt5$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp5survey_date) 
if (class(dt5$littoral_type)!="factor") dt5$littoral_type<- as.factor(dt5$littoral_type)
if (class(dt5$littoral_amount)=="factor") dt5$littoral_amount <-as.numeric(levels(dt5$littoral_amount))[as.integer(dt5$littoral_amount) ]               
if (class(dt5$littoral_amount)=="character") dt5$littoral_amount <-as.numeric(dt5$littoral_amount)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt5)                            
attach(dt5)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(littoral_type)
summary(littoral_amount) 
# Get more details on character variables

summary(as.factor(dt5$littoral_type))
detach(dt5)               


inUrl6  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/42e5d4dd5b9900828180db3c70549650" 
infile6 <- tempfile()
try(download.file(inUrl6,infile6,method="curl"))
if (is.na(file.size(infile6))) download.file(inUrl6,infile6,method="auto")


dt6 <-read.csv(infile6,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "survey_date",     
                 "shoreline_type",     
                 "shoreline_amount"    ), check.names=TRUE)

unlink(infile6)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt6$lake_id)=="factor") dt6$lake_id <-as.numeric(levels(dt6$lake_id))[as.integer(dt6$lake_id) ]               
if (class(dt6$lake_id)=="character") dt6$lake_id <-as.numeric(dt6$lake_id)                                   
# attempting to convert dt6$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp6survey_date<-as.Date(dt6$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp6survey_date) == length(tmp6survey_date[!is.na(tmp6survey_date)])){dt6$survey_date <- tmp6survey_date } else {print("Date conversion failed for dt6$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp6survey_date) 
if (class(dt6$shoreline_type)!="factor") dt6$shoreline_type<- as.factor(dt6$shoreline_type)
if (class(dt6$shoreline_amount)=="factor") dt6$shoreline_amount <-as.numeric(levels(dt6$shoreline_amount))[as.integer(dt6$shoreline_amount) ]               
if (class(dt6$shoreline_amount)=="character") dt6$shoreline_amount <-as.numeric(dt6$shoreline_amount)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt6)                            
attach(dt6)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(shoreline_type)
summary(shoreline_amount) 
# Get more details on character variables

summary(as.factor(dt6$shoreline_type))
detach(dt6)               


inUrl7  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/eae835e83c0494a376229f254f7d3392" 
infile7 <- tempfile()
try(download.file(inUrl7,infile7,method="curl"))
if (is.na(file.size(infile7))) download.file(inUrl7,infile7,method="auto")


dt7 <-read.csv(infile7,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "survey_date",     
                 "stream_nos",     
                 "stream_type",     
                 "stream_width",     
                 "stream_depth",     
                 "stream_fish_presence",     
                 "stream_barrier_distance",     
                 "stream_spawning_habitat",     
                 "stream_spawning_fish",     
                 "stream_redds",     
                 "stream_fry",     
                 "assoc_lake_spawning_habitat"    ), check.names=TRUE)

unlink(infile7)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt7$lake_id)=="factor") dt7$lake_id <-as.numeric(levels(dt7$lake_id))[as.integer(dt7$lake_id) ]               
if (class(dt7$lake_id)=="character") dt7$lake_id <-as.numeric(dt7$lake_id)                                   
# attempting to convert dt7$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp7survey_date<-as.Date(dt7$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp7survey_date) == length(tmp7survey_date[!is.na(tmp7survey_date)])){dt7$survey_date <- tmp7survey_date } else {print("Date conversion failed for dt7$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp7survey_date) 
if (class(dt7$stream_nos)=="factor") dt7$stream_nos <-as.numeric(levels(dt7$stream_nos))[as.integer(dt7$stream_nos) ]               
if (class(dt7$stream_nos)=="character") dt7$stream_nos <-as.numeric(dt7$stream_nos)
if (class(dt7$stream_type)!="factor") dt7$stream_type<- as.factor(dt7$stream_type)
if (class(dt7$stream_width)=="factor") dt7$stream_width <-as.numeric(levels(dt7$stream_width))[as.integer(dt7$stream_width) ]               
if (class(dt7$stream_width)=="character") dt7$stream_width <-as.numeric(dt7$stream_width)
if (class(dt7$stream_depth)=="factor") dt7$stream_depth <-as.numeric(levels(dt7$stream_depth))[as.integer(dt7$stream_depth) ]               
if (class(dt7$stream_depth)=="character") dt7$stream_depth <-as.numeric(dt7$stream_depth)
if (class(dt7$stream_fish_presence)!="factor") dt7$stream_fish_presence<- as.factor(dt7$stream_fish_presence)
if (class(dt7$stream_barrier_distance)=="factor") dt7$stream_barrier_distance <-as.numeric(levels(dt7$stream_barrier_distance))[as.integer(dt7$stream_barrier_distance) ]               
if (class(dt7$stream_barrier_distance)=="character") dt7$stream_barrier_distance <-as.numeric(dt7$stream_barrier_distance)
if (class(dt7$stream_spawning_habitat)=="factor") dt7$stream_spawning_habitat <-as.numeric(levels(dt7$stream_spawning_habitat))[as.integer(dt7$stream_spawning_habitat) ]               
if (class(dt7$stream_spawning_habitat)=="character") dt7$stream_spawning_habitat <-as.numeric(dt7$stream_spawning_habitat)
if (class(dt7$stream_spawning_fish)!="factor") dt7$stream_spawning_fish<- as.factor(dt7$stream_spawning_fish)
if (class(dt7$stream_redds)!="factor") dt7$stream_redds<- as.factor(dt7$stream_redds)
if (class(dt7$stream_fry)!="factor") dt7$stream_fry<- as.factor(dt7$stream_fry)
if (class(dt7$assoc_lake_spawning_habitat)=="factor") dt7$assoc_lake_spawning_habitat <-as.numeric(levels(dt7$assoc_lake_spawning_habitat))[as.integer(dt7$assoc_lake_spawning_habitat) ]               
if (class(dt7$assoc_lake_spawning_habitat)=="character") dt7$assoc_lake_spawning_habitat <-as.numeric(dt7$assoc_lake_spawning_habitat)

# Convert Missing Values to NA for non-dates

dt7$stream_width <- ifelse((trimws(as.character(dt7$stream_width))==trimws("NA")),NA,dt7$stream_width)               
suppressWarnings(dt7$stream_width <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt7$stream_width))==as.character(as.numeric("NA"))),NA,dt7$stream_width))
dt7$stream_depth <- ifelse((trimws(as.character(dt7$stream_depth))==trimws("NA")),NA,dt7$stream_depth)               
suppressWarnings(dt7$stream_depth <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt7$stream_depth))==as.character(as.numeric("NA"))),NA,dt7$stream_depth))
dt7$stream_barrier_distance <- ifelse((trimws(as.character(dt7$stream_barrier_distance))==trimws("NA")),NA,dt7$stream_barrier_distance)               
suppressWarnings(dt7$stream_barrier_distance <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt7$stream_barrier_distance))==as.character(as.numeric("NA"))),NA,dt7$stream_barrier_distance))
dt7$stream_spawning_habitat <- ifelse((trimws(as.character(dt7$stream_spawning_habitat))==trimws("NA")),NA,dt7$stream_spawning_habitat)               
suppressWarnings(dt7$stream_spawning_habitat <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt7$stream_spawning_habitat))==as.character(as.numeric("NA"))),NA,dt7$stream_spawning_habitat))
dt7$assoc_lake_spawning_habitat <- ifelse((trimws(as.character(dt7$assoc_lake_spawning_habitat))==trimws("NA")),NA,dt7$assoc_lake_spawning_habitat)               
suppressWarnings(dt7$assoc_lake_spawning_habitat <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt7$assoc_lake_spawning_habitat))==as.character(as.numeric("NA"))),NA,dt7$assoc_lake_spawning_habitat))


# Here is the structure of the input data frame:
str(dt7)                            
attach(dt7)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(stream_nos)
summary(stream_type)
summary(stream_width)
summary(stream_depth)
summary(stream_fish_presence)
summary(stream_barrier_distance)
summary(stream_spawning_habitat)
summary(stream_spawning_fish)
summary(stream_redds)
summary(stream_fry)
summary(assoc_lake_spawning_habitat) 
# Get more details on character variables

summary(as.factor(dt7$stream_type)) 
summary(as.factor(dt7$stream_fish_presence)) 
summary(as.factor(dt7$stream_spawning_fish)) 
summary(as.factor(dt7$stream_redds)) 
summary(as.factor(dt7$stream_fry))
detach(dt7)               


inUrl8  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/9fd9f9ccd630cd4b6894051c35710572" 
infile8 <- tempfile()
try(download.file(inUrl8,infile8,method="curl"))
if (is.na(file.size(infile8))) download.file(inUrl8,infile8,method="auto")


dt8 <-read.csv(infile8,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "survey_date",     
                 "site_code_nbr",     
                 "air_temp",     
                 "air_temp_time",     
                 "water_temp",     
                 "water_temp_time",     
                 "wind",     
                 "sun",     
                 "lake_max_depth",     
                 "zoo_sample_ind",     
                 "zoo_sample_time",     
                 "zoo_tow_number",     
                 "zoo_tow_type",     
                 "zoo_tow_depth",     
                 "benthic_sample_ind",     
                 "benthic_sample_percent",     
                 "nbr_benthic_sweeps",     
                 "lake_fairy_shrimp_ind",     
                 "lake_shrimp_collection",     
                 "pool_fairy_shrimp_ind",     
                 "pool_shrimp_collection",     
                 "amphib_survey_starttime",     
                 "amphib_survey_endtime",     
                 "amphib_survey_duration",     
                 "amphib_survey_desc",     
                 "amphib_survey_fish_presence",     
                 "actual_fish_presence",     
                 "fish_survey_type",     
                 "fish_net_location_type",     
                 "fish_net_set_datetime",     
                 "fish_net_pull_datetime"    ), check.names=TRUE)

unlink(infile8)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt8$lake_id)=="factor") dt8$lake_id <-as.numeric(levels(dt8$lake_id))[as.integer(dt8$lake_id) ]               
if (class(dt8$lake_id)=="character") dt8$lake_id <-as.numeric(dt8$lake_id)                                   
# attempting to convert dt8$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp8survey_date<-as.Date(dt8$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp8survey_date) == length(tmp8survey_date[!is.na(tmp8survey_date)])){dt8$survey_date <- tmp8survey_date } else {print("Date conversion failed for dt8$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp8survey_date) 
if (class(dt8$site_code_nbr)!="factor") dt8$site_code_nbr<- as.factor(dt8$site_code_nbr)
if (class(dt8$air_temp)=="factor") dt8$air_temp <-as.numeric(levels(dt8$air_temp))[as.integer(dt8$air_temp) ]               
if (class(dt8$air_temp)=="character") dt8$air_temp <-as.numeric(dt8$air_temp)
if (class(dt8$water_temp)=="factor") dt8$water_temp <-as.numeric(levels(dt8$water_temp))[as.integer(dt8$water_temp) ]               
if (class(dt8$water_temp)=="character") dt8$water_temp <-as.numeric(dt8$water_temp)
if (class(dt8$wind)!="factor") dt8$wind<- as.factor(dt8$wind)
if (class(dt8$sun)!="factor") dt8$sun<- as.factor(dt8$sun)
if (class(dt8$lake_max_depth)=="factor") dt8$lake_max_depth <-as.numeric(levels(dt8$lake_max_depth))[as.integer(dt8$lake_max_depth) ]               
if (class(dt8$lake_max_depth)=="character") dt8$lake_max_depth <-as.numeric(dt8$lake_max_depth)
if (class(dt8$zoo_sample_ind)!="factor") dt8$zoo_sample_ind<- as.factor(dt8$zoo_sample_ind)
if (class(dt8$zoo_tow_number)=="factor") dt8$zoo_tow_number <-as.numeric(levels(dt8$zoo_tow_number))[as.integer(dt8$zoo_tow_number) ]               
if (class(dt8$zoo_tow_number)=="character") dt8$zoo_tow_number <-as.numeric(dt8$zoo_tow_number)
if (class(dt8$zoo_tow_type)!="factor") dt8$zoo_tow_type<- as.factor(dt8$zoo_tow_type)
if (class(dt8$zoo_tow_depth)=="factor") dt8$zoo_tow_depth <-as.numeric(levels(dt8$zoo_tow_depth))[as.integer(dt8$zoo_tow_depth) ]               
if (class(dt8$zoo_tow_depth)=="character") dt8$zoo_tow_depth <-as.numeric(dt8$zoo_tow_depth)
if (class(dt8$benthic_sample_ind)!="factor") dt8$benthic_sample_ind<- as.factor(dt8$benthic_sample_ind)
if (class(dt8$benthic_sample_percent)=="factor") dt8$benthic_sample_percent <-as.numeric(levels(dt8$benthic_sample_percent))[as.integer(dt8$benthic_sample_percent) ]               
if (class(dt8$benthic_sample_percent)=="character") dt8$benthic_sample_percent <-as.numeric(dt8$benthic_sample_percent)
if (class(dt8$nbr_benthic_sweeps)=="factor") dt8$nbr_benthic_sweeps <-as.numeric(levels(dt8$nbr_benthic_sweeps))[as.integer(dt8$nbr_benthic_sweeps) ]               
if (class(dt8$nbr_benthic_sweeps)=="character") dt8$nbr_benthic_sweeps <-as.numeric(dt8$nbr_benthic_sweeps)
if (class(dt8$lake_fairy_shrimp_ind)!="factor") dt8$lake_fairy_shrimp_ind<- as.factor(dt8$lake_fairy_shrimp_ind)
if (class(dt8$lake_shrimp_collection)!="factor") dt8$lake_shrimp_collection<- as.factor(dt8$lake_shrimp_collection)
if (class(dt8$pool_fairy_shrimp_ind)!="factor") dt8$pool_fairy_shrimp_ind<- as.factor(dt8$pool_fairy_shrimp_ind)
if (class(dt8$pool_shrimp_collection)!="factor") dt8$pool_shrimp_collection<- as.factor(dt8$pool_shrimp_collection)
if (class(dt8$amphib_survey_duration)=="factor") dt8$amphib_survey_duration <-as.numeric(levels(dt8$amphib_survey_duration))[as.integer(dt8$amphib_survey_duration) ]               
if (class(dt8$amphib_survey_duration)=="character") dt8$amphib_survey_duration <-as.numeric(dt8$amphib_survey_duration)
if (class(dt8$amphib_survey_desc)!="factor") dt8$amphib_survey_desc<- as.factor(dt8$amphib_survey_desc)
if (class(dt8$amphib_survey_fish_presence)!="factor") dt8$amphib_survey_fish_presence<- as.factor(dt8$amphib_survey_fish_presence)
if (class(dt8$actual_fish_presence)!="factor") dt8$actual_fish_presence<- as.factor(dt8$actual_fish_presence)
if (class(dt8$fish_survey_type)!="factor") dt8$fish_survey_type<- as.factor(dt8$fish_survey_type)
if (class(dt8$fish_net_location_type)!="factor") dt8$fish_net_location_type<- as.factor(dt8$fish_net_location_type)                                   
# attempting to convert dt8$fish_net_set_datetime dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d %H:%M:%S â 07" 
tmp8fish_net_set_datetime<-as.POSIXct(dt8$fish_net_set_datetime,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp8fish_net_set_datetime) == length(tmp8fish_net_set_datetime[!is.na(tmp8fish_net_set_datetime)])){dt8$fish_net_set_datetime <- tmp8fish_net_set_datetime } else {print("Date conversion failed for dt8$fish_net_set_datetime. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp8fish_net_set_datetime)                                    
# attempting to convert dt8$fish_net_pull_datetime dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d %H:%M:%S â 07" 
tmp8fish_net_pull_datetime<-as.POSIXct(dt8$fish_net_pull_datetime,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp8fish_net_pull_datetime) == length(tmp8fish_net_pull_datetime[!is.na(tmp8fish_net_pull_datetime)])){dt8$fish_net_pull_datetime <- tmp8fish_net_pull_datetime } else {print("Date conversion failed for dt8$fish_net_pull_datetime. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp8fish_net_pull_datetime) 

# Convert Missing Values to NA for non-dates

dt8$air_temp <- ifelse((trimws(as.character(dt8$air_temp))==trimws("NA")),NA,dt8$air_temp)               
suppressWarnings(dt8$air_temp <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt8$air_temp))==as.character(as.numeric("NA"))),NA,dt8$air_temp))
dt8$water_temp <- ifelse((trimws(as.character(dt8$water_temp))==trimws("NA")),NA,dt8$water_temp)               
suppressWarnings(dt8$water_temp <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt8$water_temp))==as.character(as.numeric("NA"))),NA,dt8$water_temp))
dt8$wind <- as.factor(ifelse((trimws(as.character(dt8$wind))==trimws("NA")),NA,as.character(dt8$wind)))
dt8$sun <- as.factor(ifelse((trimws(as.character(dt8$sun))==trimws("NA")),NA,as.character(dt8$sun)))
dt8$lake_max_depth <- ifelse((trimws(as.character(dt8$lake_max_depth))==trimws("NA")),NA,dt8$lake_max_depth)               
suppressWarnings(dt8$lake_max_depth <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt8$lake_max_depth))==as.character(as.numeric("NA"))),NA,dt8$lake_max_depth))
dt8$zoo_sample_ind <- as.factor(ifelse((trimws(as.character(dt8$zoo_sample_ind))==trimws("NA")),NA,as.character(dt8$zoo_sample_ind)))
dt8$zoo_tow_number <- ifelse((trimws(as.character(dt8$zoo_tow_number))==trimws("NA")),NA,dt8$zoo_tow_number)               
suppressWarnings(dt8$zoo_tow_number <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt8$zoo_tow_number))==as.character(as.numeric("NA"))),NA,dt8$zoo_tow_number))
dt8$zoo_tow_type <- as.factor(ifelse((trimws(as.character(dt8$zoo_tow_type))==trimws("NA")),NA,as.character(dt8$zoo_tow_type)))
dt8$zoo_tow_depth <- ifelse((trimws(as.character(dt8$zoo_tow_depth))==trimws("NA")),NA,dt8$zoo_tow_depth)               
suppressWarnings(dt8$zoo_tow_depth <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt8$zoo_tow_depth))==as.character(as.numeric("NA"))),NA,dt8$zoo_tow_depth))
dt8$benthic_sample_ind <- as.factor(ifelse((trimws(as.character(dt8$benthic_sample_ind))==trimws("NA")),NA,as.character(dt8$benthic_sample_ind)))
dt8$benthic_sample_percent <- ifelse((trimws(as.character(dt8$benthic_sample_percent))==trimws("NA")),NA,dt8$benthic_sample_percent)               
suppressWarnings(dt8$benthic_sample_percent <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt8$benthic_sample_percent))==as.character(as.numeric("NA"))),NA,dt8$benthic_sample_percent))
dt8$nbr_benthic_sweeps <- ifelse((trimws(as.character(dt8$nbr_benthic_sweeps))==trimws("NA")),NA,dt8$nbr_benthic_sweeps)               
suppressWarnings(dt8$nbr_benthic_sweeps <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt8$nbr_benthic_sweeps))==as.character(as.numeric("NA"))),NA,dt8$nbr_benthic_sweeps))
dt8$lake_fairy_shrimp_ind <- as.factor(ifelse((trimws(as.character(dt8$lake_fairy_shrimp_ind))==trimws("NA")),NA,as.character(dt8$lake_fairy_shrimp_ind)))
dt8$pool_fairy_shrimp_ind <- as.factor(ifelse((trimws(as.character(dt8$pool_fairy_shrimp_ind))==trimws("NA")),NA,as.character(dt8$pool_fairy_shrimp_ind)))
dt8$amphib_survey_duration <- ifelse((trimws(as.character(dt8$amphib_survey_duration))==trimws("NA")),NA,dt8$amphib_survey_duration)               
suppressWarnings(dt8$amphib_survey_duration <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt8$amphib_survey_duration))==as.character(as.numeric("NA"))),NA,dt8$amphib_survey_duration))
dt8$amphib_survey_fish_presence <- as.factor(ifelse((trimws(as.character(dt8$amphib_survey_fish_presence))==trimws("NA")),NA,as.character(dt8$amphib_survey_fish_presence)))
dt8$actual_fish_presence <- as.factor(ifelse((trimws(as.character(dt8$actual_fish_presence))==trimws("NA")),NA,as.character(dt8$actual_fish_presence)))
dt8$fish_survey_type <- as.factor(ifelse((trimws(as.character(dt8$fish_survey_type))==trimws("NA")),NA,as.character(dt8$fish_survey_type)))
dt8$fish_net_location_type <- as.factor(ifelse((trimws(as.character(dt8$fish_net_location_type))==trimws("NA")),NA,as.character(dt8$fish_net_location_type)))


# Here is the structure of the input data frame:
str(dt8)                            
attach(dt8)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(site_code_nbr)
summary(air_temp)
summary(air_temp_time)
summary(water_temp)
summary(water_temp_time)
summary(wind)
summary(sun)
summary(lake_max_depth)
summary(zoo_sample_ind)
summary(zoo_sample_time)
summary(zoo_tow_number)
summary(zoo_tow_type)
summary(zoo_tow_depth)
summary(benthic_sample_ind)
summary(benthic_sample_percent)
summary(nbr_benthic_sweeps)
summary(lake_fairy_shrimp_ind)
summary(lake_shrimp_collection)
summary(pool_fairy_shrimp_ind)
summary(pool_shrimp_collection)
summary(amphib_survey_starttime)
summary(amphib_survey_endtime)
summary(amphib_survey_duration)
summary(amphib_survey_desc)
summary(amphib_survey_fish_presence)
summary(actual_fish_presence)
summary(fish_survey_type)
summary(fish_net_location_type)
summary(fish_net_set_datetime)
summary(fish_net_pull_datetime) 
# Get more details on character variables

summary(as.factor(dt8$site_code_nbr)) 
summary(as.factor(dt8$wind)) 
summary(as.factor(dt8$sun)) 
summary(as.factor(dt8$zoo_sample_ind)) 
summary(as.factor(dt8$zoo_tow_type)) 
summary(as.factor(dt8$benthic_sample_ind)) 
summary(as.factor(dt8$lake_fairy_shrimp_ind)) 
summary(as.factor(dt8$lake_shrimp_collection)) 
summary(as.factor(dt8$pool_fairy_shrimp_ind)) 
summary(as.factor(dt8$pool_shrimp_collection)) 
summary(as.factor(dt8$amphib_survey_desc)) 
summary(as.factor(dt8$amphib_survey_fish_presence)) 
summary(as.factor(dt8$actual_fish_presence)) 
summary(as.factor(dt8$fish_survey_type)) 
summary(as.factor(dt8$fish_net_location_type))
detach(dt8)               


inUrl9  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/a343e504c5c40a2a6c9cfc989cdcb829" 
infile9 <- tempfile()
try(download.file(inUrl9,infile9,method="curl"))
if (is.na(file.size(infile9))) download.file(inUrl9,infile9,method="auto")


dt9 <-read.csv(infile9,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "lake_id",     
                 "survey_date",     
                 "survey_type",     
                 "crew_id"    ), check.names=TRUE)

unlink(infile9)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt9$lake_id)=="factor") dt9$lake_id <-as.numeric(levels(dt9$lake_id))[as.integer(dt9$lake_id) ]               
if (class(dt9$lake_id)=="character") dt9$lake_id <-as.numeric(dt9$lake_id)                                   
# attempting to convert dt9$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp9survey_date<-as.Date(dt9$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp9survey_date) == length(tmp9survey_date[!is.na(tmp9survey_date)])){dt9$survey_date <- tmp9survey_date } else {print("Date conversion failed for dt9$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp9survey_date) 
if (class(dt9$survey_type)!="factor") dt9$survey_type<- as.factor(dt9$survey_type)
if (class(dt9$crew_id)=="factor") dt9$crew_id <-as.numeric(levels(dt9$crew_id))[as.integer(dt9$crew_id) ]               
if (class(dt9$crew_id)=="character") dt9$crew_id <-as.numeric(dt9$crew_id)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt9)                            
attach(dt9)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(survey_type)
summary(crew_id) 
# Get more details on character variables

summary(as.factor(dt9$survey_type))
detach(dt9)               


inUrl10  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/227238b58704514f38cfa290ba901407" 
infile10 <- tempfile()
try(download.file(inUrl10,infile10,method="curl"))
if (is.na(file.size(infile10))) download.file(inUrl10,infile10,method="auto")


dt10 <-read.csv(infile10,header=F 
                ,skip=1
                ,sep=","  
                ,quot='"' 
                , col.names=c(
                  "lake_id",     
                  "survey_date",     
                  "Subsample",     
                  "SpeciesID",     
                  "Number"    ), check.names=TRUE)

unlink(infile10)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt10$lake_id)=="factor") dt10$lake_id <-as.numeric(levels(dt10$lake_id))[as.integer(dt10$lake_id) ]               
if (class(dt10$lake_id)=="character") dt10$lake_id <-as.numeric(dt10$lake_id)                                   
# attempting to convert dt10$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp10survey_date<-as.Date(dt10$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp10survey_date) == length(tmp10survey_date[!is.na(tmp10survey_date)])){dt10$survey_date <- tmp10survey_date } else {print("Date conversion failed for dt10$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp10survey_date) 
if (class(dt10$Subsample)=="factor") dt10$Subsample <-as.numeric(levels(dt10$Subsample))[as.integer(dt10$Subsample) ]               
if (class(dt10$Subsample)=="character") dt10$Subsample <-as.numeric(dt10$Subsample)
if (class(dt10$SpeciesID)=="factor") dt10$SpeciesID <-as.numeric(levels(dt10$SpeciesID))[as.integer(dt10$SpeciesID) ]               
if (class(dt10$SpeciesID)=="character") dt10$SpeciesID <-as.numeric(dt10$SpeciesID)
if (class(dt10$Number)=="factor") dt10$Number <-as.numeric(levels(dt10$Number))[as.integer(dt10$Number) ]               
if (class(dt10$Number)=="character") dt10$Number <-as.numeric(dt10$Number)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt10)                            
attach(dt10)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(Subsample)
summary(SpeciesID)
summary(Number) 
# Get more details on character variables

detach(dt10)               


inUrl11  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/d4d60ccd5034909faabadd8cccb14603" 
infile11 <- tempfile()
try(download.file(inUrl11,infile11,method="curl"))
if (is.na(file.size(infile11))) download.file(inUrl11,infile11,method="auto")


dt11 <-read.csv(infile11,header=F 
                ,skip=1
                ,sep=","  
                ,quot='"' 
                , col.names=c(
                  "lake_id",     
                  "survey_date",     
                  "ID",     
                  "SpeciesID",     
                  "Length"    ), check.names=TRUE)

unlink(infile11)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt11$lake_id)=="factor") dt11$lake_id <-as.numeric(levels(dt11$lake_id))[as.integer(dt11$lake_id) ]               
if (class(dt11$lake_id)=="character") dt11$lake_id <-as.numeric(dt11$lake_id)                                   
# attempting to convert dt11$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp11survey_date<-as.Date(dt11$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp11survey_date) == length(tmp11survey_date[!is.na(tmp11survey_date)])){dt11$survey_date <- tmp11survey_date } else {print("Date conversion failed for dt11$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp11survey_date) 
if (class(dt11$ID)=="factor") dt11$ID <-as.numeric(levels(dt11$ID))[as.integer(dt11$ID) ]               
if (class(dt11$ID)=="character") dt11$ID <-as.numeric(dt11$ID)
if (class(dt11$SpeciesID)=="factor") dt11$SpeciesID <-as.numeric(levels(dt11$SpeciesID))[as.integer(dt11$SpeciesID) ]               
if (class(dt11$SpeciesID)=="character") dt11$SpeciesID <-as.numeric(dt11$SpeciesID)
if (class(dt11$Length)=="factor") dt11$Length <-as.numeric(levels(dt11$Length))[as.integer(dt11$Length) ]               
if (class(dt11$Length)=="character") dt11$Length <-as.numeric(dt11$Length)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt11)                            
attach(dt11)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(ID)
summary(SpeciesID)
summary(Length) 
# Get more details on character variables

detach(dt11)               


inUrl12  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/48967dffec165077f7130f7a13025fba" 
infile12 <- tempfile()
try(download.file(inUrl12,infile12,method="curl"))
if (is.na(file.size(infile12))) download.file(inUrl12,infile12,method="auto")


dt12 <-read.csv(infile12,header=F 
                ,skip=1
                ,sep=","  
                ,quot='"' 
                , col.names=c(
                  "SpeciesID",     
                  "Species_Name"    ), check.names=TRUE)

unlink(infile12)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt12$SpeciesID)=="factor") dt12$SpeciesID <-as.numeric(levels(dt12$SpeciesID))[as.integer(dt12$SpeciesID) ]               
if (class(dt12$SpeciesID)=="character") dt12$SpeciesID <-as.numeric(dt12$SpeciesID)
if (class(dt12$Species_Name)!="factor") dt12$Species_Name<- as.factor(dt12$Species_Name)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt12)                            
attach(dt12)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(SpeciesID)
summary(Species_Name) 
# Get more details on character variables

summary(as.factor(dt12$Species_Name))
detach(dt12)               


inUrl13  <- "https://pasta.lternet.edu/package/data/eml/edi/577/2/2e45eac50efe844d2b1017e74e024627" 
infile13 <- tempfile()
try(download.file(inUrl13,infile13,method="curl"))
if (is.na(file.size(infile13))) download.file(inUrl13,infile13,method="auto")


dt13 <-read.csv(infile13,header=F 
                ,skip=1
                ,sep=","  
                ,quot='"' 
                , col.names=c(
                  "lake_id",     
                  "survey_date",     
                  "collect_date",     
                  "sample_select",     
                  "process_date",     
                  "process_name",     
                  "subsample",     
                  "sample_vol",     
                  "sample_type",     
                  "magnif",     
                  "notes"    ), check.names=TRUE)

unlink(infile13)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt13$lake_id)=="factor") dt13$lake_id <-as.numeric(levels(dt13$lake_id))[as.integer(dt13$lake_id) ]               
if (class(dt13$lake_id)=="character") dt13$lake_id <-as.numeric(dt13$lake_id)                                   
# attempting to convert dt13$survey_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp13survey_date<-as.Date(dt13$survey_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp13survey_date) == length(tmp13survey_date[!is.na(tmp13survey_date)])){dt13$survey_date <- tmp13survey_date } else {print("Date conversion failed for dt13$survey_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp13survey_date)                                    
# attempting to convert dt13$collect_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp13collect_date<-as.Date(dt13$collect_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp13collect_date) == length(tmp13collect_date[!is.na(tmp13collect_date)])){dt13$collect_date <- tmp13collect_date } else {print("Date conversion failed for dt13$collect_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp13collect_date) 
if (class(dt13$sample_select)!="factor") dt13$sample_select<- as.factor(dt13$sample_select)                                   
# attempting to convert dt13$process_date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d %H:%M:%S â 07" 
tmp13process_date<-as.POSIXct(dt13$process_date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp13process_date) == length(tmp13process_date[!is.na(tmp13process_date)])){dt13$process_date <- tmp13process_date } else {print("Date conversion failed for dt13$process_date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp13process_date) 
if (class(dt13$process_name)=="factor") dt13$process_name <-as.numeric(levels(dt13$process_name))[as.integer(dt13$process_name) ]               
if (class(dt13$process_name)=="character") dt13$process_name <-as.numeric(dt13$process_name)
if (class(dt13$subsample)=="factor") dt13$subsample <-as.numeric(levels(dt13$subsample))[as.integer(dt13$subsample) ]               
if (class(dt13$subsample)=="character") dt13$subsample <-as.numeric(dt13$subsample)
if (class(dt13$sample_vol)=="factor") dt13$sample_vol <-as.numeric(levels(dt13$sample_vol))[as.integer(dt13$sample_vol) ]               
if (class(dt13$sample_vol)=="character") dt13$sample_vol <-as.numeric(dt13$sample_vol)
if (class(dt13$sample_type)!="factor") dt13$sample_type<- as.factor(dt13$sample_type)
if (class(dt13$magnif)!="factor") dt13$magnif<- as.factor(dt13$magnif)
if (class(dt13$notes)!="factor") dt13$notes<- as.factor(dt13$notes)

# Convert Missing Values to NA for non-dates

dt13$sample_type <- as.factor(ifelse((trimws(as.character(dt13$sample_type))==trimws("NA")),NA,as.character(dt13$sample_type)))
dt13$magnif <- as.factor(ifelse((trimws(as.character(dt13$magnif))==trimws("NA")),NA,as.character(dt13$magnif)))


# Here is the structure of the input data frame:
str(dt13)                            
attach(dt13)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(lake_id)
summary(survey_date)
summary(collect_date)
summary(sample_select)
summary(process_date)
summary(process_name)
summary(subsample)
summary(sample_vol)
summary(sample_type)
summary(magnif)
summary(notes) 
# Get more details on character variables

summary(as.factor(dt13$sample_select)) 
summary(as.factor(dt13$sample_type)) 
summary(as.factor(dt13$magnif)) 
summary(as.factor(dt13$notes))
detach(dt13)               





