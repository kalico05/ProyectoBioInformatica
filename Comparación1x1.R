library(seqinr)
library(ggplot2)
library(dplyr)

trad = c(
  UUU = "F", UUC = "F", UUA = "L", UUG = "L",
  UCU = "S", UCC = "S", UCA = "S", UCG = "S",
  UAU = "Y", UAC = "Y", UAA = "STOP", UAG = "STOP",
  UGU = "C", UGC = "C", UGA = "STOP", UGG = "W",
  CUU = "L", CUC = "L", CUA = "L", CUG = "L",
  CCU = "P", CCC = "P", CCA = "P", CCG = "P",
  CAU = "H", CAC = "H", CAA = "Q", CAG = "Q",
  CGU = "R", CGC = "R", CGA = "R", CGG = "R",
  AUU = "I", AUC = "I", AUA = "I", AUG = "M",
  ACU = "T", ACC = "T", ACA = "T", ACG = "T",
  AAU = "N", AAC = "N", AAA = "K", AAG = "K",
  AGU = "S", AGC = "S", AGA = "R", AGG = "R",
  GUU = "V", GUC = "V", GUA = "V", GUG = "V",
  GCU = "A", GCC = "A", GCA = "A", GCG = "A",
  GAU = "D", GAC = "D", GAA = "E", GAG = "E",
  GGU = "G", GGC = "G", GGA = "G", GGG = "G"
)

# Función para convertir a ARN
FuncionConvertirARN = function(entrada) {
  nuevoARN = as.vector(entrada)
  nuevoARN[which(nuevoARN == "t")] = "u"
  nuevoARN = toupper(nuevoARN)
  return(nuevoARN)
}

# Leer secuencias
archivoReferen = read.fasta("sequence_wuhan.txt")
archivoMexico = read.fasta("sequencesB1.fasta")

# Inicializar dataframe
datosFrame = data.frame(
  Mutacion = character(),
  Codon = character(),
  Aminoacido = character(),
  Posicion = character(),
  Gen = character()
)

# Procesar secuencias
for (i in seq(1, length(archivoReferen), 1)) {
  anotaciones = attr(archivoReferen[[i]], "Annot")
  atributos = strsplit(anotaciones, "=")[[1]]
  nombreGen = atributos[which(grepl("gene", atributos)) + 1]
  genRef = FuncionConvertirARN(archivoReferen[[i]])
  for(k in seq(i, length(archivoMexico), 12)) {
    genMex = FuncionConvertirARN(archivoMexico[[k]])
    if (length(genRef) == length(genMex)) {
      diferencia = which(genRef != genMex)
      if (length(diferencia) > 0) {
        for (x in diferencia) {
          mutacion = paste(genRef[x], "to", genMex[x], sep = "")
          inicioCodigo = x - (x - 1) %% 3
          numeroCodigo = as.integer((x - 1) / 3 + 1)
          codonOri = paste(genRef[inicioCodigo], genRef[inicioCodigo + 1], genRef[inicioCodigo + 2], sep = "")
          codonMex = paste(genMex[inicioCodigo], genMex[inicioCodigo + 1], genMex[inicioCodigo + 2], sep = "")
          cambioCodigo = paste(codonOri, "to", codonMex, sep = "")
          cambioAminoAcido = paste(trad[codonOri], numeroCodigo, trad[codonMex], sep = "")
          if (!is.na(trad[codonMex])) {
            datosFrame <- rbind(datosFrame, data.frame(
              Mutacion = mutacion,
              Codon = cambioCodigo,
              Aminoacido = cambioAminoAcido,
              Posicion = x,
              Gen = nombreGen
            ))
          }
        }
      }
    } else {
      cat("Mutación de deleción o inserción\n")
    }
  }
}

# Gráfica de cambio de nucleótido
ggplot(datosFrame, aes(x = Mutacion, fill = Mutacion)) +
  geom_bar() +
  xlab("Cambio de nucleótido") +
  ylab("Frecuencia") +
  ggtitle("Gráfica de cambio de nucleótido")

# Gráfica de cambio de aminoácido
ggplot(datosFrame, aes(x = Aminoacido, fill = Aminoacido)) +
  geom_bar() +
  xlab("Cambio de aminoácido") +
  ylab("Frecuencia") +
  ggtitle("Gráfica de cambio de aminoácido")