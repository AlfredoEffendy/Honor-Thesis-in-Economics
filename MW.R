library(foreign)
library(plm)
setwd("C:/Education/University of California Irvine/Winter 2018/ECON 190BW/Data")
borderpairs = read.dta("borderpairs.dta")
contigstatepairs = read.dta("contigstatepairs.dta")
CPSdata = read.dta("CPSdata.dta")
QCEWindustry_minwage_all = read.dta("QCEWindustry_minwage_all.dta")
QCEWindustry_minwage_contig = read.dta("QCEWindustry_minwage_contig.dta")

time = c(1:86)
time = rep(time,51)
CPSdata$date = time

teenOLS = lm(log(epteen) ~ log(minwage) + urall + rpteen, data=CPSdata)
teenFE = lm(log(epteen) ~ log(minwage) + urall + rpteen + factor(state) +
              factor(date), data=CPSdata)

teenFE2 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=CPSdata,
              index=c("state","date"), model = "within")
teenRE = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=CPSdata,
              index=c("state","date"), model = "random")

teenCCE = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=CPSdata,
               index=c("state","date"), model = "mg")

teenCCEP = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=CPSdata,
               index=c("state","date"), model = "p")

urallbar = rep(0,86)
for(t in (1:86)){
  urall = rep(0,51)
  for(j in (1:51)){
    urall[j] = CPSdata$urall[t+86*(j-1)]
  }
  urallbar[t] = mean(urall)
}

rpteenbar = rep(0,86)
for(t in (1:86)){
  rpteen = rep(0,51)
  for(j in (1:51)){
    rpteen[j] = CPSdata$rpteen[t+86*(j-1)]
  }
  rpteenbar[t] = mean(rpteen)
}

epteenbar = rep(0,86)
for(t in (1:86)){
  epteen = rep(0,51)
  for(j in (1:51)){
    epteen[j] = CPSdata$epteen[t+86*(j-1)]
  }
  epteenbar[t] = mean(epteen)
}

minwagebar = rep(0,86)
for(t in (1:86)){
  minwage = rep(0,51)
  for(j in (1:51)){
    minwage[j] = CPSdata$minwage[t+86*(j-1)]
  }
  minwagebar[t] = mean(minwage)
}

urallbar = rep(urallbar,51)
rpteenbar = rep(rpteenbar,51)
epteenbar = rep(epteenbar,51)
minwagebar = rep(minwagebar,51)

CPSdata$urallbar = urallbar
CPSdata$rpteenbar = rpteenbar
CPSdata$epteenbar = epteenbar
CPSdata$minwagebar = minwagebar

#Half

firsthalf = CPSdata[1:43,]
for(i in (1:50)){
  initial = i*86 + 1
  end = i*86 + 43
  newentry = CPSdata[initial:end,]
  firsthalf = rbind(firsthalf,newentry)
}

teenOLS1 = lm(log(epteen) ~ log(minwage) + urall + rpteen, data=firsthalf)
teenFE1 = lm(log(epteen) ~ log(minwage) + urall + rpteen + factor(state) +
              factor(date), data=firsthalf)
teenFE21 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=firsthalf,
              index=c("state","date"), model = "within")
teenRE1 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=firsthalf,
             index=c("state","date"), model = "random")

teenCCE1 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=firsthalf,
               index=c("state","date"), model = "mg")

teenCCEP1 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=firsthalf,
                index=c("state","date"), model = "p")

secondhalf = CPSdata[44:86,]
for(i in (1:50)){
  initial = i*86 + 44
  end = i*86 + 86
  newentry = CPSdata[initial:end,]
  secondhalf = rbind(secondhalf,newentry)
}

teenOLS2 = lm(log(epteen) ~ log(minwage) + urall + rpteen, data=secondhalf)
teenFE2 = lm(log(epteen) ~ log(minwage) + urall + rpteen + factor(state) +
               factor(date), data=secondhalf)
teenFE22 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=secondhalf,
               index=c("state","date"), model = "within")
teenRE2 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=secondhalf,
              index=c("state","date"), model = "random")

teenCCE2 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=secondhalf,
                index=c("state","date"), model = "mg")

teenCCEP2 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=secondhalf,
                 index=c("state","date"), model = "p")

#Quarter

firstquarter = CPSdata[1:22,]
for(i in (1:50)){
  initial = i*86 + 1
  end = i*86 + 22
  newentry = CPSdata[initial:end,]
  firstquarter = rbind(firstquarter,newentry)
}

teenOLS3 = lm(log(epteen) ~ log(minwage) + urall + rpteen, data=firstquarter)
teenFE3 = lm(log(epteen) ~ log(minwage) + urall + rpteen + factor(state) +
               factor(date), data=firstquarter)
teenFE23 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=firstquarter,
               index=c("state","date"), model = "within")
teenRE3 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=firstquarter,
              index=c("state","date"), model = "random")

teenCCE3 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=firstquarter,
                index=c("state","date"), model = "mg")

teenCCEP3 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=firstquarter,
                 index=c("state","date"), model = "p")


secondquarter = CPSdata[23:43,]
for(i in (1:50)){
  initial = i*86 + 23
  end = i*86 + 43
  newentry = CPSdata[initial:end,]
  secondquarter = rbind(secondquarter,newentry)
}

teenOLS4 = lm(log(epteen) ~ log(minwage) + urall + rpteen, data=secondquarter)
teenFE4 = lm(log(epteen) ~ log(minwage) + urall + rpteen + factor(state) +
               factor(date), data=secondquarter)
teenFE24 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=secondquarter,
               index=c("state","date"), model = "within")
teenRE4 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=secondquarter,
              index=c("state","date"), model = "random")

teenCCE4 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=secondquarter,
                index=c("state","date"), model = "mg")

teenCCEP4 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=secondquarter,
                 index=c("state","date"), model = "p")

thirdquarter = CPSdata[44:65,]
for(i in (1:50)){
  initial = i*86 + 44
  end = i*86 + 65
  newentry = CPSdata[initial:end,]
  thirdquarter = rbind(thirdquarter,newentry)
}

teenOLS5 = lm(log(epteen) ~ log(minwage) + urall + rpteen, data=thirdquarter)
teenFE5 = lm(log(epteen) ~ log(minwage) + urall + rpteen + factor(state) +
               factor(date), data=thirdquarter)
teenFE25 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=thirdquarter,
               index=c("state","date"), model = "within")
teenRE5 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=thirdquarter,
              index=c("state","date"), model = "random")

teenCCE5 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=thirdquarter,
                index=c("state","date"), model = "mg")

teenCCEP5 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=thirdquarter,
                 index=c("state","date"), model = "p")

fourthquarter = CPSdata[66:86,]
for(i in (1:50)){
  initial = i*86 + 66
  end = i*86 + 86
  newentry = CPSdata[initial:end,]
  fourthquarter = rbind(fourthquarter,newentry)
}

teenOLS6 = lm(log(epteen) ~ log(minwage) + urall + rpteen, data=fourthquarter)
teenFE6 = lm(log(epteen) ~ log(minwage) + urall + rpteen + factor(state) +
               factor(date), data=fourthquarter)
teenFE26 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=fourthquarter,
               index=c("state","date"), model = "within")
teenRE6 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=fourthquarter,
              index=c("state","date"), model = "random")

teenCCE6 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=fourthquarter,
                index=c("state","date"), model = "mg")

teenCCEP6 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=fourthquarter,
                 index=c("state","date"), model = "p")

last5 = CPSdata[82:86,]
for(i in (1:50)){
  initial = i*86 + 82
  end = i*86 + 86
  newentry = CPSdata[initial:end,]
  last5 = rbind(last5,newentry)
}

teenOLS7 = lm(log(epteen) ~ log(minwage) + urall + rpteen, data=last5)
teenFE7 = lm(log(epteen) ~ log(minwage) + urall + rpteen + factor(state) +
               factor(date), data=last5)
teenFE27 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=last5,
               index=c("state","date"), model = "within")
teenRE7 = plm(log(epteen) ~ log(minwage) + urall + rpteen, data=last5,
              index=c("state","date"), model = "random")

teenCCE7 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=last5,
                index=c("state","date"), model = "mg")

teenCCEP7 = pcce(log(epteen) ~ log(minwage) + urall + rpteen, data=last5,
                 index=c("state","date"), model = "p")

