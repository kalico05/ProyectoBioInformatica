library(seqinr)

trad =    c(UUU="F", UUC="F", UUA="L", UUG="L",
            UCU="S", UCC="S", UCA="S", UCG="S",
            UAU="Y", UAC="Y", UAA="STOP", UAG="STOP",
            UGU="C", UGC="C", UGA="STOP", UGG="W",
            CUU="L", CUC="L", CUA="L", CUG="L",
            CCU="P", CCC="P", CCA="P", CCG="P",
            CAU="H", CAC="H", CAA="Q", CAG="Q",
            CGU="R", CGC="R", CGA="R", CGG="R",
            AUU="I", AUC="I", AUA="I", AUG="M",
            ACU="T", ACC="T", ACA="T", ACG="T",
            AAU="N", AAC="N", AAA="K", AAG="K",
            AGU="S", AGC="S", AGA="R", AGG="R",
            GUU="V", GUC="V", GUA="V", GUG="V",
            GCU="A", GCC="A", GCA="A", GCG="A",
            GAU="D", GAC="D", GAA="E", GAG="E",
            GGU="G", GGC="G", GGA="G", GGG="G")

ToARN = function(input) {
  newARN = as.vector(input)
  newARN[which(newARN=="t")] = "u"
  newARN = toupper(newARN)
  return (newARN)
}

df = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

fRef = read.fasta("sequence_wuhan.txt")
length(fRef)

fB117 = read.fasta("sequencesP1.fasta")
length(fB117)

cat("Procesando ", as.integer(length(fB117)/12), " genomas \n")

nObs = 1
for (i in seq(1,length(fRef),1)){
  if (i==2) next
  anotaciones = attr(fRef[[i]], "Annot") 
  atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
  geneName = atributos[which(atributos=="gene")+1] 
  genRef = ToARN( fRef[[i]] )  
  cat("#",geneName)
  for (k in seq(i, length(fB117), 12)){
    
     genfB117 = ToARN( fB117[[k]] )  
     cat(i, k, length(genRef), length(genfB117), "\n")
     
    if (length(genRef) == length(genfB117)){
      dif = which(genRef != genfB117)
      cat("length",length(dif))
      print(dif)
    
      if (length(dif) > 0){
        for (x in dif){
          cat(i,k)
          muta = paste(genRef[x],"to",genfB117[x], sep="") 
          inicioCodon = x - (x-1)%%3 
          numCodon = as.integer((x-1)/3+1) 
          codonOri = paste(genRef[inicioCodon], genRef[inicioCodon+1], genRef[inicioCodon+2],sep="")
          codonfB117 = paste(genfB117[inicioCodon], genfB117[inicioCodon+1], genfB117[inicioCodon+2],sep="")
          codonChange = paste(codonOri,"to",codonfB117, sep="")
          aminoChange = paste(trad[codonOri],numCodon,trad[codonfB117], sep="")
          cat(i,k,geneName, codonChange, aminoChange)
          
          if (!is.na(trad[codonfB117]) && trad[codonOri]!=trad[codonfB117]){
            obs = list(muta,codonChange,aminoChange,geneName)
            df[nObs,] = obs 
            nObs = nObs+1
          }
        }
      }
    }else{
    }
  }
}

head(df)
nrow(df)

library(ggplot2)

p = ggplot(df)
p = p + aes(x=Mutation, fill=Mutation, label=after_stat(count))
p = p + ggtitle("Mutaciones de sustitución")
p = p + labs(x="Mutation", y="Frecuencia", fill="Frecuencia")
p = p + geom_bar(stat = "count")
p = p + geom_text(stat = "count", vjust=1.5)
#p = p + facet_grid(~Gene)
p

library(dplyr)
dfgraph = filter(
  summarise(
    select(
      group_by(df, Amino),
      Mutation:Gene
    ),
    Mutation = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta>1
)

dfgraph = dfgraph[order(-dfgraph$Cuenta), ]
dfgraph = dfgraph[1:20, ]
head(dfgraph)
nrow(dfgraph)

p2 = ggplot(dfgraph)
p2 = p2 + aes(x=Amino, y=Cuenta, fill=Amino, label=Cuenta)
p2 = p2 + ggtitle("Cambio de Aminoácidos")
p2 = p2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
p2 = p2 + geom_bar(stat = "identity")
p2 = p2 + geom_text(stat = "identity", vjust=1.5)
p2 = p2 + facet_grid(~Gene, scales="free", space="free_x")
p2



