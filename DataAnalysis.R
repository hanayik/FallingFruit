# Data analysis for Falling Fruit task
library(Hmisc)
# intialize grouping variables arrays
allConHitAcc = c()
allConTotAcc = c()
allConTotNumStim = c()
allConHits = c()
allConTargsHit = c()
allConNumTargsShown = c()
allConTargetAcc = c()

allCussHitAcc = c()
allCussTotAcc = c()
allCussTotNumStim = c()
allCussHits = c()
allCussTargsHit = c()
allCussNumTargsShown = c()
allCussTargetAcc = c()
######################  con group   ######################################
baseFolder = "/Users/rorden/Documents/MATLAB/FallingFruit-master (3)/data/con"
fileList = list.files(baseFolder, pattern = ".csv", full.names = TRUE)
numberOfFiles = length(fileList)
for (i in 1:numberOfFiles) {
  file = fileList[i]
  print(file)
  data = read.csv(file)
  # get number of targets seen
  onlyTargs = length(data$isTarget[data$isTarget == 1])
  # get number of hits
  onlyHits = length(data$hit[data$hit == 1])
  targsHit = length(data$hit[data$hit==1 & data$isTarget==1])
  onlyApples = length(data$word[data$word == "apple"])
  hitAcc = round((targsHit/onlyHits)*100, digits = 2)
  #hit accuracy
  totalAcc = round((targsHit/onlyApples)*100, digits = 2)
  totalNumOfStimuli = length(data$hit)
  # get left side accuracy
  leftOnlyHits = length(data$hit[data$hit == 1 & data$scrXpos < 1422])
  leftTargsHit = length(data$hit[data$hit==1 & data$isTarget==1 & data$scrXpos < 1422])
  leftOnlyApples = length(data$word[data$word == "apple" & data$scrXpos < 1422])
  leftAcc = round((leftTargsHit/leftOnlyHits)*100, digits = 2)
  leftTotalAcc = round((leftTargsHit/leftOnlyApples)*100, digits = 2)
  # get right side accurcy
  rightOnlyHits = length(data$hit[data$hit == 1 & data$scrXpos > 1422])
  rightTargsHit = length(data$hit[data$hit==1 & data$isTarget==1 & data$scrXpos > 1422])
  rightOnlyApples = length(data$word[data$word == "apple" & data$scrXpos > 1422])
  rightAcc = round((rightTargsHit/rightOnlyHits)*100, digits = 2)
  rightTotalAcc = round((rightTargsHit/rightOnlyApples)*100, digits = 2)
  
  allConHitAcc = append(allConHitAcc, hitAcc)
  allConTotAcc = append(allConTotAcc, totalAcc)
  allConTotNumStim = append(allConTotNumStim, totalNumOfStimuli)
  allConHits = append(allConHits, onlyHits)
  allConTargsHit = append(allConTargsHit, targsHit)
  allConNumTargsShown = append(allConNumTargsShown, onlyTargs)
  allConTargetAcc = append(allConTargetAcc, targsHit/onlyTargs)
  
  
}

dataForCard_zConFruit = scale(allConTotAcc)

######################  cuss group   ######################################
# Data analysis for Falling Fruit task
baseFolder = "/Users/rorden/Documents/MATLAB/FallingFruit-master (3)/data/cuss"
fileList = list.files(baseFolder, pattern = ".csv", full.names = TRUE)
numberOfFiles = length(fileList)

for (i in 1:numberOfFiles) {

  file = fileList[i]
  data = read.csv(file)
  # get only targs shown
  onlyTargs = length(data$isTarget[data$isTarget == 1])
  # get number of hits
  onlyHits = length(data$hit[data$hit == 1])
  targsHit = length(data$hit[data$hit==1 & data$isTarget==1])
  onlyApples = length(data$word[data$word == "apple"])
  hitAcc = round((targsHit/onlyHits)*100, digits = 2)
  #hit accuracy
  totalAcc = round((targsHit/onlyApples)*100, digits = 2)
  totalNumOfStimuli = length(data$hit)
  # get left side accuracy
  leftOnlyHits = length(data$hit[data$hit == 1 & data$scrXpos < 1422])
  leftTargsHit = length(data$hit[data$hit==1 & data$isTarget==1 & data$scrXpos < 1422])
  leftOnlyApples = length(data$word[data$word == "apple" & data$scrXpos < 1422])
  leftAcc = round((leftTargsHit/leftOnlyHits)*100, digits = 2)
  leftTotalAcc = round((leftTargsHit/leftOnlyApples)*100, digits = 2)
  # get right side accurcy
  rightOnlyHits = length(data$hit[data$hit == 1 & data$scrXpos > 1422])
  rightTargsHit = length(data$hit[data$hit==1 & data$isTarget==1 & data$scrXpos > 1422])
  rightOnlyApples = length(data$word[data$word == "apple" & data$scrXpos > 1422])
  rightAcc = round((rightTargsHit/rightOnlyHits)*100, digits = 2)
  rightTotalAcc = round((rightTargsHit/rightOnlyApples)*100, digits = 2)
 
  allCussHitAcc = append(allCussHitAcc, hitAcc)
  allCussTotAcc = append(allCussTotAcc, totalAcc)
  allCussTotNumStim = append(allCussTotNumStim, totalNumOfStimuli)
  # collect number of hits and targets hit
  allCussHits = append(allCussHits, onlyHits)
  allCussTargsHit = append(allCussTargsHit, targsHit)
  allCussNumTargsShown = append(allCussNumTargsShown, onlyTargs)
  allCussTargetAcc = append(allCussTargetAcc, targsHit/onlyTargs)
  
  
}

dataForCard_zCussFruit = scale(allCussTotAcc)

t = t.test(allConTotAcc, allCussTotAcc, var.equal = TRUE)
print("TotAcc")
t

t = t.test(allConHitAcc, allCussHitAcc, var.equal = TRUE)
print("HitAcc")
t

t = t.test(allConTotNumStim, allCussTotNumStim, var.equal = TRUE)
print("Number of stims seen")
t

t = t.test(allConTargsHit, allCussTargsHit, var.equal = TRUE)
print("Targets Hit")
t

t = t.test(allConHits, allCussHits, var.equal = TRUE)
print("Just Hits")
t

t = t.test(allConNumTargsShown, allCussNumTargsShown, var.equal = TRUE)
print("Targets shown")
t

t = t.test(allConTargetAcc, allCussTargetAcc, var.equal = TRUE)
print("Targets hit out of targets shown")
t

mat4corr = matrix(c(allConNumTargsShown, allConTargsHit),nrow=length(allConNumTargsShown))
print("Correlation of targets seen and targets hit")
rcorr(mat4corr, type="pearson")

mat4corr = matrix(c(allCussNumTargsShown, allCussTargsHit),nrow=length(allCussNumTargsShown))
print("Correlation of targets seen and targets hit")
rcorr(mat4corr, type="pearson")




