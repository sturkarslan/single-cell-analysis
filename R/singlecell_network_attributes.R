# This script creates network attributes for single cell mutation nodes

#load early geenration mutations
ua3.b = read.delim("~/Google Drive/Single-Cell-Genomics/variants/EPD/UA3_variants.FINAL-dvh.txt", sep="\t", header=F, stringsAsFactors = F)
ua3.b$line = "ua3.b"
ua3.03 = read.delim("~/Google Drive/Single-Cell-Genomics/variants/EPD/UA3_03_variants.FINAL-dvh.txt", sep="\t", header=F, stringsAsFactors = F)
ua3.03$line = "ua3.03"
ua3.09 = read.delim("~/Google Drive/Single-Cell-Genomics/variants/EPD/UA3_09_variants.FINAL-dvh.txt", sep="\t", header=F, stringsAsFactors = F)
ua3.09$line = "ua3.09"
ua3.10 = read.delim("~/Google Drive/Single-Cell-Genomics/variants/after_300/UA3-10_variants.FINAL-dvh.txt", sep="\t", header=F, stringsAsFactors = F)
ua3.10$line = "ua3.10"
ua3.15 = read.delim("~/Google Drive/Single-Cell-Genomics/variants/after_300/UA3-15_variants.FINAL-dvh.txt", sep="\t", header=F, stringsAsFactors = F)
ua3.15$line = "ua3.15"
ua3.45 = read.delim("~/Google Drive/Single-Cell-Genomics/variants/after_300/UA3-45_variants.FINAL-dvh.txt", sep="\t", header=F, stringsAsFactors = F)
ua3.45$line = "ua3.45"
ua3.76 = read.delim("~/Google Drive/Single-Cell-Genomics/variants/after_300/UA3-76_variants.FINAL-dvh.txt", sep="\t", header=F, stringsAsFactors = F)
ua3.76$line = "ua3.76"
ua3.118 = read.delim("~/Google Drive/Single-Cell-Genomics/variants/after_300/UA3-118_variants.FINAL-dvh.txt", sep="\t", header=F, stringsAsFactors = F)
ua3.118$line = "ua3.118"

# Load dvh genome data files for Dvh-UA3-152-03</h3>
essential = read.delim("~/Google Drive/Portals/Snytrophy_Portal/dvh-essentiality-data.txt", header=T, sep="\t", stringsAsFactors=F)
genome = read.delim("~/Google Drive/Portals/Snytrophy_Portal/dvu_genomeInfo.txt", header=T, sep="\t", stringsAsFactors=F)

# read original mutation data
dvh03.original = read.delim("~/Google Drive/Single-Cell-Genomics/variants/singlecell/dvh-03-single-cell-Final-Merged-Variant-Filtered2.txt", sep="\t", header=F)

## function to format names
generation.names = list(ua3.b,ua3.03, ua3.09, ua3.10, ua3.15, ua3.45, ua3.76, ua3.118)
generations = data.frame()
for(generation in generation.names){
  line = unique(generation$line)
  for(row in 1:length(generation$V1)){

    if("" %in% generation[row,"V11"]){
      name = sub("Chromosome", "", paste(paste("IG", generation[row, "V2"], sep="_"), sub("pDV", "p", generation[row, "V1"] ), sep = ""))
    } else {
      name =sub("Chromosome", "", paste(sub("DVU_", "DVU", paste(generation[row, "V11"], generation[row, "V2"], sep="_")), sub("pDV", "p", generation[row, "V1"] ), sep = ""))
    }
    generations = rbind(generations, cbind(generation = line, names = name))
  }
  
}
# convert data frame into a list
generation.list = list()
for(line in unique(generations$generation)){
  generation.list[[line]] = as.vector(generations[which(generations$generation == line),"names"])
}

# append unique mutation id to each row
for(row in 1:length(dvh03.original$V1)){
  if("" %in% dvh03.original[row,"V11"]){
    dvh03.original[row,"ID"] = paste("IG", dvh03.original[row, "V2"], sep="_")
  } else {
    dvh03.original[row,"ID"] = sub("DVU_", "DVU", paste(dvh03.original[row, "V11"], dvh03.original[row, "V2"], sep="_"))
  }
}


# read mutations bedfile
dvh03.bedfile = read.delim("/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/dvh-UA3-152-03-singlecell-variants-2callers-80percent-2cells_noan-bed.txt", sep="\t", header=F)
for(row in 1:length(dvh03.bedfile$V1)){
  dvh03.bedfile[row,"ID"] = sub("DVU_", "DVU", paste(dvh03.bedfile[row,"V7"], dvh03.bedfile[row,"V2"], sep="_"))
}

#load clonal isolate naming file
clonal.isolate.names = read.delim("/Volumes/omics4tb/sturkarslan/clonal-isolates/clonal_isolate_names.txt", sep="\t", header=T, stringsAsFactors = F)
# read clonal isolate mutations
dvhUA3.isolates = read.delim("/Volumes/omics4tb/sturkarslan/clonal-isolates/results/dvh/clonal-isolates_2callers-filtered-variants.txt", sep="\t", header=F, stringsAsFactors = F)
dvhUA3.isolates$cloneno = sapply(dvhUA3.isolates$V18, function(x) strsplit(x, split = "_")[[1]][1])


for(i in 1:length(dvhUA3.isolates$V1)){
  id = dvhUA3.isolates[i,"cloneno"]
  dvhUA3.isolates[i,"clonename"] = clonename = clonal.isolate.names[which(clonal.isolate.names$pair == id),"isolate"][1]
  dvhUA3.isolates[i,"line"] = strsplit(clonename, split=".", fixed=T)[[1]][1]
  dvhUA3.isolates[i,"epd"] = strsplit(clonename, split=".", fixed=T)[[1]][3]
  dvhUA3.isolates[i,"clone"] = line = strsplit(clonename, split=".", fixed=T)[[1]][4]
  
  if("" %in% dvhUA3.isolates[i,"V11"]){
    dvhUA3.isolates[i,"ID"] = sub("Chromosome", "", paste(paste("IG", dvhUA3.isolates[i, "V2"], sep="_"), sub("pDV", "p", dvhUA3.isolates[i, "V1"] ), sep = ""))
  } else {
    dvhUA3.isolates[i,"ID"] = sub("Chromosome", "", paste(sub("DVU_", "DVU", paste(dvhUA3.isolates[i, "V11"], dvhUA3.isolates[i, "V2"], sep="_")), sub("pDV", "p", dvhUA3.isolates[i, "V1"] ), sep = ""))
  }
  
}

# selet only UA3/03 epd line variants
dvh03.isolates = dvhUA3.isolates[which(dvhUA3.isolates$line == "UA3" & dvhUA3.isolates$epd == "03"),]

# read mutation matrix
dvh03.matrix = read.delim("/Volumes/omics4tb/sturkarslan/scite_single_cells_syntrophy/dvh-UA3-152-03_noan_mutation_counts_verified_5cells_50nas_mutation_matrix.txt", sep="\t", header=F, stringsAsFactors = F)

# read mutation names
dvh03.names = read.delim("/Volumes/omics4tb/sturkarslan/scite_single_cells_syntrophy/dvh-UA3-152-03_noan_mutation_counts_verified_5cells_50nas_mutation_names.txt", sep="\t", header=F)



# Attach attributes to mutations</h3>
dvh03.features = data.frame()
for(name in c(dvh03.names)[[1]]){
  cat("Now analyzing ", name, "...\n")
  cat("\n")
  locus = strsplit(name, split = "_", fixed = T)[[1]][1]
  name.us = sub("p", "", name)
  # name.us = sub("DVUA", "DVXA", name.us)
  # name.us = sub("DVU", "DVU_", name.us)
  # name.us = sub("DVXA", "DVUA", name.us)
  cat(name.us,"\n")
  
  if(name.us %in% dvh03.original$ID){
    # mutations
    type = as.character(dvh03.original[which(dvh03.original$ID == name.us), "V6"])
    impact = as.character(dvh03.original[which(dvh03.original$ID == name.us), "V7"])
    effect = as.character(dvh03.original[which(dvh03.original$ID == name.us), "V8"])
    # info
    
    if(locus == "IG"){
      gene.name = ""
      gene.desc = ""
    } else if(length(genome[which(genome$sysName == locus), "name"]) == 0){
      gene.name = ""
      gene.desc = ""
    } else {
      gene.name = genome[which(genome$sysName == locus), "name"]
      gene.desc = genome[which(genome$sysName == locus), "desc"]
    }
    
    # essentialiy
    moyls4.1 = essential[which(essential$locus_tag == locus),"WT.MOYLS4.1"]
    mols4 = essential[which(essential$locus_tag == locus),"WT.MOLS4"]
    
    if(length(moyls4.1) == 0){
      moyls4.1 = ""
    } else {
      moyls4.1 = moyls4.1
    }
    if(length(mols4) == 0){
      mols4 = ""
    } else {
      mols4 = mols4
    }
    
    # go terms
    go = genome[which(genome$sysName == locus), "GO"]
    if (length(go) == 0){
      go = ""
    }
    if(go == ""){
      go = ""
    } else {
      go = paste(go)
    }
    
    # COGFun
    cog = genome[which(genome$sysName == locus), "COGFun"]
    if(length(cog) == 0){
      cog = ""
    }
    if(cog == ""){
      cog = ""
    } else {
      cog.length = length(strsplit(cog, split = "")[[1]])
      if(cog.length > 1){
        cog = sapply(cog, function(x) paste(strsplit(x, split = "")[[1]][1], strsplit(x, split = "")[[1]][2], sep = ":" ))
      } else {
        cog = cog
      }
      
    }
    
    # accession
    accession = genome[which(genome$sysName == locus), "accession"]
    if(length(accession) == 0){
      accession = ""
    }
    if(accession == ""){
      accession = ""
    } else {
      accession = accession
    }
    
    # GI
    GI = as.character(genome[which(genome$sysName == locus), "GI"])
    if(length(GI) == 0){
      GI = ""
    } else {
      GI = GI
    }
    
    # number of cells w/ mutation, w/o and NA
    rownumber = grep(name, dvh03.names$V1)
    in.cells = length(grep(1, dvh03.matrix[rownumber,]))
    notin.cells = length(grep(0, dvh03.matrix[rownumber,]))
    NA.cells = length(grep(3, dvh03.matrix[rownumber,]))
    
  } else {
    cat(name.us, "Not found..\n")
  }
  ## check the clonal isolates file
  if(name %in% dvh03.isolates$ID){
    clones.1 = dvh03.isolates[which(dvh03.isolates$ID == name),"clone"]
    clones = paste(clones.1, collapse = ":", sep="")
    clone.count = length(clones.1)
  } else {
    cat(name, "Not found in clonal isolates..\n")
    clones = ""
    clone.count = 0
  }
  
  ## check if mutation is in early generations
  early.gen = paste(names(generation.list)[grep(name, generation.list)], sep = "", collapse = ":")
  gen.names = c("ua3.b","ua3.03", "ua3.09", "ua3.10", "ua3.15", "ua3.45", "ua3.76", "ua3.118")
  
  k.list = list()
  for(k in gen.names){
        k.list[[k]] = length(grep(k, early.gen))
  }
  
  count.ua3.b = k.list[["ua3.b"]]
  count.ua3.03 = k.list[["ua3.03"]]
  count.ua3.09 = k.list[["ua3.09"]]
  count.ua3.10 = k.list[["ua3.10"]]
  count.ua3.15 = k.list[["ua3.15"]]
  count.ua3.45 = k.list[["ua3.45"]]
  count.ua3.76 = k.list[["ua3.76"]]
  count.ua3.118 = k.list[["ua3.118"]]
  
  

  dvh03.features = rbind(dvh03.features, 
                         cbind(name = name,
                               locus = locus, gene.name = gene.name, gene.desc = gene.desc,
                               type = type, impact = impact, effect=effect,
                               moyls4 = moyls4.1, mols4=mols4,
                               go = go, cog = cog, accession = accession, GI = GI,
                               in.cells = in.cells, notin.cells = notin.cells, NA.cells = NA.cells, 
                               clones = clones, clone.count = clone.count, early.gen = early.gen,
                               count.152 = count.ua3.b, count.03 = count.ua3.03, count.09 = count.ua3.09, count.10 = count.ua3.10,
                               count.15 = count.ua3.15, count.45 = count.ua3.45, count.76 = count.ua3.76, count.118 = count.ua3.118))
  cat(name,"|",locus,"|",gene.name,"|",gene.desc,"|",name.us,"|",type,"|",impact,"|",effect,"|",mols4,"|",moyls4.1,"|",go,"|",cog,"|",accession,"|",GI,"|",in.cells,"|",notin.cells,"|",NA.cells,"\n")
  
}

write.table(dvh03.features, file="/Volumes/omics4tb/sturkarslan/scite_single_cells_syntrophy/dvh-UA3-152-03_noan_5cells_50nas-network-attributes.txt", sep="\t", row.names=F, quote=F)
