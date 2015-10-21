tbc = read.csv("rawmilk/data/tbc_data.csv")

tbc$Farm = factor(tbc$Farm)
tbc$Month = factor(tbc$Month)
tbc$Volume = tbc$Volume / 1000

# omit data that can't be used in the full model
tbc = na.omit(tbc)

# The below seems to give the table, with the exception of log(Volume) which seems wrong (maybe as offset is not included)?
mod1 = glmer(TBC ~ log(Volume) + log(CC+1) + log(SCC) + as.factor(Month) + (1|Farm), family="poisson", data=tbc)
mod2 = glmer(TBC ~ log(Volume) + log(CC+1) + as.factor(Month) + (1|Farm), family="poisson", data=tbc)
mod3 = glmer(TBC ~ log(Volume) + log(SCC) + as.factor(Month) + (1|Farm), family="poisson", data=tbc)

y   = tbc$TBC
mu1 = predict(mod1, type="response")
mu2 = predict(mod2, type="response")
mu3 = predict(mod3, type="response")

xlogx = function(x)
{
  ifelse(x==0, 0, x*log(x))
}

R2_scc  = 1 - sum(xlogx(y)-y - (y*log(mu1)-mu1)) / sum(xlogx(y)-y - (y*log(mu2)-mu2))
R2_coli = 1 - sum(xlogx(y)-y - (y*log(mu1)-mu1)) / sum(xlogx(y)-y - (y*log(mu3)-mu3))
