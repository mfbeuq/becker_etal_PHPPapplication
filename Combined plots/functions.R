rename.ancombc.output <- function(tax.clean)
{
  for(i in 1:nrow(tax.clean)){
    for (j in 1:5){
      tax.clean[,j][is.na(tax.clean[,j])] <- FALSE
    }
    if(tax.clean$Family[i] == as.character(tax.clean$Order[i])) {
      tax.clean$Name[i] <- paste("Unclassified", tax.clean$Name[i], sep = " ")
    } else if (tax.clean$Family[i] == as.character(tax.clean$Class[i])) {
      tax.clean$Name[i] <- paste("Unclassified", tax.clean$Name[i], sep = " ")
    } else if (tax.clean$Family[i] == as.character(tax.clean$Phylum[i])) {
      tax.clean$Name[i] <- paste("Unclassified", tax.clean$Name[i], sep = " ")
    } else (print = "Hi")
    if(tax.clean$Name[i] == "Ellin6067") {
      tax.clean$Name[i] <- paste("Nitrosomonadaceae", tax.clean$Name[i], sep = " ")
    }
    if(tax.clean$Name[i] == "MND1") {
      tax.clean$Name[i] <- paste("Nitrosomonadaceae", tax.clean$Name[i], sep = " ")
    }
    if(tax.clean$Name[i] == "Unclassified JG36-TzT-191") {
      tax.clean$Name[i] <- "Unclassified Gammaproteo. JG36-TzT-191"
    }
    if(tax.clean$Name[i] == "SC-I-84") {
      tax.clean$Name[i] <- "Unclassified Burkholderiales SC-I-84"
    }
    if(tax.clean$Name[i] == "Vermiphilaceae") {
      tax.clean$Name[i] <- "Unclassified Vermiphilaceae"
    }
    if(tax.clean$Name[i] == "AKYH767") {
      tax.clean$Name[i] <- "Unclassified Sphingobacteriales AKYH767"
    }
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Saprospiraceae") {
        tax.clean$Name[i] <- "Unclassified Saprospiraceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Comamonadaceae") {
        tax.clean$Name[i] <- "Unclassified Comamonadaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Xanthomonadaceae") {
        tax.clean$Name[i] <- "Unclassified Xanthomonadaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Halieaceae") {
        tax.clean$Name[i] <- "Unclassified Halieaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "RhizobialesIncertaeSedis") {
        tax.clean$Name[i] <- "Unclassified Rhizobiales Incertae Sedis"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Rhodocyclaceae") {
        tax.clean$Name[i] <- "Unclassified Rhodocyclaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Steroidobacteraceae") {
        tax.clean$Name[i] <- "Unclassified Steroidobacteraceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Enterobacteriaceae") {
        tax.clean$Name[i] <- "Unclassified Enterobacteriaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Solimonadaceae") {
        tax.clean$Name[i] <- "Unclassified Solimonadaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Sphingomonadaceae") {
        tax.clean$Name[i] <- "Unclassified Sphingomonadaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Chitinophagaceae") {
        tax.clean$Name[i] <- "Unclassified Chitinophagaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Oxalobacteraceae") {
        tax.clean$Name[i] <- "Unclassified Oxalobacteraceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Rhodanobacteraceae") {
        tax.clean$Name[i] <- "Unclassified Rhodanobacteraceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Nitrosomonadaceae") {
        tax.clean$Name[i] <- "Unclassified Nitrosomonadaceae"
      }}
    if(tax.clean$Genus[i] == "Unclassified") {
      if(tax.clean$Family[i] == "Ilumatobacteraceae") {
        tax.clean$Name[i] <- "Unclassified Ilumatobacteraceae"
      }}
    if(tax.clean$Name[i] == "Unclassified Unclassified") {
      if(tax.clean$Class[i] == "Actinobacteria") {
        tax.clean$Name[i] <- "Unclassified Actinobacteria"
      }}
    if(tax.clean$Order[i] == "Microtrichales") {
      if(tax.clean$Order[i] == "uncultured") {
        tax.clean$Name[i] <- "Unclassified Microtrichales"
      } else if(tax.clean$Order[i] == "FALSE") {
        tax.clean$Name[i] <- "Unclassified Microtrichales"
      }}
    if(tax.clean$Name[i] == "TM7") {
      tax.clean$Name[i] <- "Unclassified Saccharimonadales TM7"
    }
    if(tax.clean$Name[i] == "A21b") {
      tax.clean$Name[i] <- "Unclassified Burkholderiales A21b"
    }
    if(tax.clean$Name[i] == "BD1-7clade") {
      tax.clean$Name[i] <- "SpongiibacteraceaeBD1-7 Clade"
    }
    if(tax.clean$Name[i] == "mle1-7") {
      tax.clean$Name[i] <- "Nitrosomonadaceae mle1-7"
    }
    if(tax.clean$Name[i] == "Unclassified PLTA13") {
      tax.clean$Name[i] <- "Unclassified Gammaproteo. PLTA13"
    }
    if(tax.clean$Name[i] == "Unclassified Subgroup5") {
      tax.clean$Name[i] <- "Unclassified Acidobacteriota Subgroup 5"
    }
    if(tax.clean$Name[i] == "Unclassified Subgroup7") {
      tax.clean$Name[i] <- "Unclassified Holophagae Subgroup 7"
    }
    if(tax.clean$Name[i] == "Unclassified Subgroup22") {
      tax.clean$Name[i] <- "Unclassified Acidobacteriota Subgroup 22"
    }
    if(tax.clean$Name[i] == "Unclassified IMCC26256") {
      tax.clean$Name[i] <- "Unclassified Acidimicrobiia IMCC26256"
    }
    if(tax.clean$Name[i] == "CL500-29marinegroup") {
      tax.clean$Name[i] <- "CL500-29 Marine Group"
    }
    if(tax.clean$Name[i] == "Unclassified SJA-28") {
      tax.clean$Name[i] <- "Unclassified Bacteroidota SJA-28"
    }
    if(tax.clean$Name[i] == "BSV26") {
      tax.clean$Name[i] <- "Unclassified Kryptoniales BSV26"
    }
    if(tax.clean$Name[i] == "Subgroup10") {
      tax.clean$Name[i] <- "Thermoanaerobaculaceae Subgroup 10"
    }
    if(tax.clean$Name[i] == "TRA3-20") {
      tax.clean$Name[i] <- "Unclassified Burkholderiales TRA3-20"
    } else(print="Hi")
    tax.clean[] <- lapply(tax.clean, gsub, pattern = "_", replacement = " ", fixed = TRUE)
  } 
  return(tax.clean)
  print(tax.clean$Name)
}