library(TeachingDemos)

# R CMD BATCH --no-save --no-restore '--args nbcores=2'
args = commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}
print(42)
print(Sys.time())
print(nbcores)


