library(seqinr)
library(dplyr)
library(ggplot2)


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

df_wuhan_vs_B117 = data.frame(
  Mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

df_wuhan_vs_P1 = data.frame(
  mutation2 = character(),
  Codon2 = character(),
  Amino2 = character(),
  Gene2 = character(),
  stringsAsFactors = FALSE
)

df_B117_vs_P1 = data.frame(
  mutation = character(),
  Codon = character(),
  Amino = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

fRef = read.fasta("sequence_wuhan.txt")
length(fRef)

fB117 = read.fasta("sequencesB117.fasta")
length(fB117)

fP1 = read.fasta("sequencesP1.fasta")
length(fP1)

cat("Procesando ", as.integer(length(fB117)/12), " genomas de B117 y P1 \n")

#CICLOS PARA B117

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
      
      if (length(dif) > 0){
        for (x in dif){
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
            df_wuhan_vs_B117[nObs,] = obs 
            nObs = nObs+1
          }
        }
      }
    }else{
    }
  }
}

head(df_wuhan_vs_B117)
nrow(df_wuhan_vs_B117)


p = ggplot(df_wuhan_vs_B117)
p = p + aes(x=Mutation, fill=Mutation, label=after_stat(count))
p = p + ggtitle("Mutaciones de sustitución Wuhan - B117")
p = p + labs(x="Mutation", y="Frecuencia", fill="Frecuencia")
p = p + geom_bar(stat = "count")
p = p + geom_text(stat = "count", vjust=1.5)
p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

df_wuhan_vs_B117graph = filter(
  summarise(
    select(
      group_by(df_wuhan_vs_B117, Amino),
      Mutation:Gene
    ),
    Mutation = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta>1
)

df_wuhan_vs_B117graph <- df_wuhan_vs_B117 %>%
  filter(!is.na(Amino)) %>%
  group_by(Amino) %>%
  summarise(
    mutation2 = first(Mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ) %>%
  filter(Cuenta > 1) %>%
  arrange(desc(Cuenta)) %>%
  slice(1:20)

head(df_wuhan_vs_B117graph)
nrow(df_wuhan_vs_B117graph)

p2 <- ggplot(df_wuhan_vs_B117graph)
p2 <- p2 + aes(x = Amino, y = Cuenta, fill = Amino, label = Cuenta)
p2 <- p2 + ggtitle("Cambio de Aminoácidos Wuhan - B117")
p2 <- p2 + labs(x = "Amino", y = "Frecuencia", fill = "Frecuencia")
p2 <- p2 + geom_bar(stat = "identity")
p2 <- p2 + geom_text(stat = "identity", vjust = 1.5)
p2 <- p2 + facet_grid(~ Gene, scales = "free", space = "free_x")
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2

# CICLOS PARA P1

nObs2 = 1
for (i in seq(1, length(fRef), 1)) {
  if (i == 2) next
  anotaciones2 = attr(fRef[[i]], "Annot") 
  atributos2 = unlist(strsplit(anotaciones2, "\\[|\\]|:|=|\\.|\\(")); 
  geneName2 = atributos2[which(atributos2 == "gene") + 1] 
  genRef2 = ToARN(fRef[[i]])  
  cat("#", geneName2)
  for (k in seq(i, length(fP1), 12)) {
    genfP1 = ToARN(fP1[[k]])  
    cat(i, k, length(genRef2), length(genfP1), "\n")
    
    if (length(genRef2) == length(genfP1)) {
      dif = which(genRef2 != genfP1)
      cat("length", length(dif))
      
      if (length(dif) > 0) {
        for (x in dif) {
          muta2 = paste(genRef2[x], "to", genfP1[x], sep = "") 
          inicioCodon2 = x - (x - 1) %% 3 
          numCodon2 = as.integer((x - 1) / 3 + 1) 
          codonOri2 = paste(genRef2[inicioCodon2], genRef2[inicioCodon2 + 1], genRef2[inicioCodon2 + 2], sep = "")
          codonfP1 = paste(genfP1[inicioCodon2], genfP1[inicioCodon2 + 1], genfP1[inicioCodon2 + 2], sep = "")
          codonChange2 = paste(codonOri2, "to", codonfP1, sep = "")
          aminoChange2 = paste(trad[codonOri2], numCodon2, trad[codonfP1], sep = "")
          cat(i, k, geneName2, codonChange2, aminoChange2)
          
          if (!is.na(trad[codonfP1]) && trad[codonOri2] != trad[codonfP1]) {
            obs2 = list(muta2, codonChange2, aminoChange2, geneName2)
            df_wuhan_vs_P1[nObs2, ] <- obs2 
            nObs2 = nObs2 + 1
          }
        }
      }
    } else {
    }
  }
}

head(df_wuhan_vs_P1)
nrow(df_wuhan_vs_P1)


p12 = ggplot(df_wuhan_vs_P1)
p12 = p12 + aes(x = mutation2, fill = mutation2, label = after_stat(count))
p12 = p12 + ggtitle("Mutaciones de sustitución Wuhan - P1")
p12 = p12 + labs(x = "Mutation", y = "Frecuencia", fill = "Frecuencia")
p12 = p12 + geom_bar(stat = "count")
p12 = p12 + geom_text(stat = "count", vjust = 1.5)
p12 = p12 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p12

dfgraph2 = filter(
  summarise(
    select(
      group_by(df_wuhan_vs_P1, Amino2),
      mutation2:Gene2
    ),
    mutation2 = first(mutation2),
    Codon = first(Codon2),
    Gene = first(Gene2),
    Cuenta = n()
  ),
  Cuenta > 1
)

dfgraph2 = dfgraph2[order(-dfgraph2$Cuenta), ]
dfgraph2 = dfgraph2[1:20, ]
head(dfgraph2)
nrow(dfgraph2)

p22 = ggplot(dfgraph2)
p22 = p22 + aes(x = Amino2, y = Cuenta, fill = Amino2, label = Cuenta)
p22 = p22 + ggtitle("Cambio de Aminoácidos Wuhan - P1")
p22 = p22 + labs(x = "Amino", y = "Frecuencia", fill = "Frecuencia")
p22 = p22 + geom_bar(stat = "identity")
p22 = p22 + geom_text(stat = "identity", vjust = 1.5)
p22 = p22 + facet_grid(~ Gene, scales = "free", space = "free_x")
p22 = p22 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p22

# Comparación B117 - P1
for (i in seq_along(fB117)) {
  # Obtener la secuencia de B117 en la posición i
  seqB117 = ToARN(fB117[[i]])
  
  # Obtener la secuencia correspondiente de P1
  seqP1 = ToARN(fP1[[i]])
  
  # Verificar que las longitudes de las secuencias sean iguales
  if (length(seqB117) == length(seqP1)) {
    # Iterar sobre los elementos de las secuencias
    for (j in seq_along(seqB117)) {
      # Comparar los elementos correspondientes de B117 y P1
      if (seqB117[j] != seqP1[j]) {
        # Encontrar la mutación y guardarla en df_B117_vs_P1
        muta = paste(toupper(seqB117[j]), "to", toupper(seqP1[j]), sep="")
        inicioCodon = j - (j-1) %% 3
        numCodon = as.integer((j-1) / 3 + 1)
        codonOri = paste(seqB117[inicioCodon:(inicioCodon + 2)], collapse = "")
        codonP1 = paste(seqP1[inicioCodon:(inicioCodon + 2)], collapse = "")
        codonChange = paste(codonOri, "to", codonP1, sep="")
        
        aminoOri = trad[codonOri]
        aminoP1 = trad[codonP1]
        cat(codonOri, codonP1, "\n")
        
        # Combinar los aminoácidos de manera adecuada
        aminoChange = paste(aminoOri, numCodon, aminoP1, sep = "")
        
        # Verificar si aminoOri y aminoP1 son NA o iguales
        if (!is.na(aminoOri) && !is.na(aminoP1) && aminoOri != aminoP1) {
          obs = list(mutation = muta, Codon = codonChange, Amino = aminoChange, Gene = "B117_P1")
          df_B117_vs_P1 = rbind(df_B117_vs_P1, obs)
        }
      }
    }
  } else {
    cat("Las longitudes de las secuencias no son iguales en la posición", i, "\n")
  }
}

head(df_B117_vs_P1)
nrow(df_B117_vs_P1)

p13 = ggplot(df_B117_vs_P1)
p13 = p13 + aes(x = mutation, fill = mutation, label = after_stat(count))
p13 = p13 + ggtitle("Mutaciones de sustitución B117 - P1")
p13 = p13 + labs(x = "Mutation", y = "Frecuencia", fill = "Frecuencia")
p13 = p13 + geom_bar(stat = "count")
p13 = p13 + geom_text(stat = "count", vjust = 1.5)
p13 = p13 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p13

dfgraph3 = filter(
  summarise(
    select(
      group_by(df_B117_vs_P1, Amino),
      mutation:Gene
    ),
    mutation = first(mutation),
    Codon = first(Codon),
    Gene = first(Gene),
    Cuenta = n()
  ),
  Cuenta > 1
)

dfgraph3 = dfgraph3[order(-dfgraph3$Cuenta), ]
dfgraph3 = dfgraph3[1:20, ]
head(dfgraph3)
nrow(dfgraph3)

p23 = ggplot(dfgraph3)
p23 = p23 + aes(x = Amino, y = Cuenta, fill = Amino, label = Cuenta)
p23 = p23 + ggtitle("Cambio de Aminoácidos B117 - P1")
p23 = p23 + labs(x = "Amino", y = "Frecuencia", fill = "Frecuencia")
p23 = p23 + geom_bar(stat = "identity")
p23 = p23 + geom_text(stat = "identity", vjust = 1.5)
p23 = p23 + facet_grid(~ Gene, scales = "free", space = "free_x")
p23 = p23 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p23

# Filtrar mutaciones comunes entre Wuhan vs B117 y Wuhan vs P1
common_mutations <- intersect(df_wuhan_vs_B117$Mutation, df_wuhan_vs_P1$mutation2)

df_wuhan_vs_B117_common <- df_wuhan_vs_B117[df_wuhan_vs_B117$Mutation %in% common_mutations, ]
df_wuhan_vs_P1_common <- df_wuhan_vs_P1[df_wuhan_vs_P1$mutation2 %in% common_mutations, ]

# Gráfica de mutaciones comunes entre Wuhan vs B117
p_common_wuhan_b117 <- ggplot(df_wuhan_vs_B117_common) +
  aes(x = Mutation, fill = Mutation, label = after_stat(count)) +
  ggtitle("Mutaciones comunes Wuhan - B117") +
  labs(x = "Mutación", y = "Frecuencia", fill = "Frecuencia") +
  geom_bar(stat = "count") +
  geom_text(stat = "count", vjust = 1.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gráfica de mutaciones comunes entre Wuhan vs P1
p_common_wuhan_p1 <- ggplot(df_wuhan_vs_P1_common) +
  aes(x = mutation2, fill = mutation2, label = after_stat(count)) +
  ggtitle("Mutaciones comunes Wuhan - P1") +
  labs(x = "Mutación", y = "Frecuencia", fill = "Frecuencia") +
  geom_bar(stat = "count") +
  geom_text(stat = "count", vjust = 1.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Mostrar las gráficas
p_common_wuhan_b117
p_common_wuhan_p1

require(gridExtra)
F1=grid.arrange(p,p12,p13,nrow=2)
F1
require(gridExtra)
F2=grid.arrange(p_common_wuhan_b117,p_common_wuhan_p1, nrow=1)
F2
require(gridExtra)
F3=grid.arrange(p2,p22,p23, nrow=3)
F3
