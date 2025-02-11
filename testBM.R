# Load species occurrences (6 species available)
data(DataSpecies)
head(DataSpecies)

# Select the name of the studied species
myRespName <- 'GuloGulo'

# Get corresponding presence/absence data
myResp <- as.numeric(DataSpecies[, myRespName])

# Get corresponding XY coordinates
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]

# Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
data(bioclim_current)
myExpl <- terra::rast(bioclim_current)



# ---------------------------------------------------------------#
# Format Data with true absences
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
myBiomodData
summary(myBiomodData)
plot(myBiomodData)


# ---------------------------------------------------------------#
# # Transform true absences into potential pseudo-absences
# myResp.PA <- ifelse(myResp == 1, 1, NA)
# 
# # Format Data with pseudo-absences : random method
# myBiomodData.r <- BIOMOD_FormatingData(resp.var = myResp.PA,
#                                        expl.var = myExpl,
#                                        resp.xy = myRespXY,
#                                        resp.name = myRespName,
#                                        PA.nb.rep = 4,
#                                        PA.nb.absences = 1000,
#                                        PA.strategy = 'random')
# 
# # Format Data with pseudo-absences : disk method
# myBiomodData.d <- BIOMOD_FormatingData(resp.var = myResp.PA,
#                                        expl.var = myExpl,
#                                        resp.xy = myRespXY,
#                                        resp.name = myRespName,
#                                        PA.nb.rep = 4,
#                                        PA.nb.absences = 500,
#                                        PA.strategy = 'disk',
#                                        PA.dist.min = 5,
#                                        PA.dist.max = 35)
# 
# # Format Data with pseudo-absences : SRE method
# myBiomodData.s <- BIOMOD_FormatingData(resp.var = myResp.PA,
#                                        expl.var = myExpl,
#                                        resp.xy = myRespXY,
#                                        resp.name = myRespName,
#                                        PA.nb.rep = 4,
#                                        PA.nb.absences = 1000,
#                                        PA.strategy = 'sre',
#                                        PA.sre.quant = 0.025)
# 
# # Format Data with pseudo-absences : user.defined method
myPAtable <- data.frame(PA1 = ifelse(myResp == 1, TRUE, FALSE),
                         PA2 = ifelse(myResp == 1, TRUE, FALSE))
for (i in 1:ncol(myPAtable)) myPAtable[sample(which(myPAtable[, i] == FALSE), 500), i] = TRUE
myBiomodData.u <- BIOMOD_FormatingData(resp.var = myResp.PA,
                                        expl.var = myExpl,
                                        resp.xy = myRespXY,
                                        resp.name = myRespName,
                                        PA.strategy = 'user.defined',
                                        PA.user.table = myPAtable)
 
# Documentation Example
class(myResp.PA) #numeric
length(myResp.PA) #2488

class(myExpl) #SpatRaster
dim(myExpl) #47x120x5

class(myRespXY) #data.frame
dim(myRespXY) #2488x2

class(myPAtable) #data.frame
dim(myPAtable) #2488x2 

# My Data
class(resp.var) #numeric
length(resp.var) #5433
print(resp.var)

class(resp.xy) #data.frame
dim(resp.xy) #5433x2
print(resp.xy)

class(PA.user.table) #data.frame
dim(PA.user.table) #5443x10
print(PA.user.table)




