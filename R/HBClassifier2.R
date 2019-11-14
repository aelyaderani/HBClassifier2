#' classification function
#'
#' this function classifies scRNA Human brain.
#' @param SampleSet
#' @param nucIsoType
#' @return SampleSet
#' @export
#' HBClassifier2()
#'
HBClassifier2 <- function(SampleSet, nucIsoType)
{require(Matrix);require(devtools);require(Seurat);require(repmis);require(dplyr);require(RCurl);require(repmis)

  SampleSet <- RunPCA(object = SampleSet, verbose = FALSE)
  SampleSet <- ProjectDim(object = SampleSet)

  v <- SampleSet@reductions$pca@stdev^2
  PCA_percentage <-(v/sum(v))*100
  rm(v)
  paste(PCA_percentage,"%", sep = " ")


  PClist <- c()
  TotalPCPer <- 0
  PCNumber <- 1

  while (TotalPCPer <= 85 && PCNumber <= 50)
  {
    TotalPCPer <- TotalPCPer + PCA_percentage[PCNumber]
    PClist <- c(PClist,PCNumber)
    PCNumber <- PCNumber + 1
  }

  PCs <- PClist
  print("Number of PCs Selected:")
  print(NROW(PCs))
  SampleSet <- FindNeighbors(object = SampleSet, dims = PCs)
  SampleSet <- FindClusters(object = SampleSet, resolution = .2)
  SampleSet <- RunUMAP(object = SampleSet, assay = "SCT", dims = PCs)
  DimPlot(object = SampleSet, reduction.use = "umap")

  cluster1.markers <- FindMarkers(object = SampleSet, ident.1 = 1, min.pct = 0.35)

  SampleSet.markers <- FindAllMarkers(object = SampleSet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  class_cluster_list <- SampleSet.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

  reflist <- getURL("https://raw.githubusercontent.com/aelyaderani/HBClassifier2/master/HBClassifier_marker_ref.csv")
  HBClassifier_marker_ref <- read.csv(text = reflist)
  unknown <- "Unknown"

  reflist <- data.frame(Gene_name= HBClassifier_marker_ref$gene, Cell_name= HBClassifier_marker_ref$cell,stringsAsFactors = FALSE)
  samplesetlist <- data.frame(genes= class_cluster_list$gene, cluster= class_cluster_list$cluster,cellname = unknown,stringsAsFactors = FALSE)



  samSize <- nrow(samplesetlist)
  refSize <- nrow(reflist)
  currentSamRow <- 1
  currentRefRow <- 1
  pb <- txtProgressBar(min = 0, max = samSize, style = 3)

  while(currentSamRow <= samSize)
  {
    while (currentRefRow <= refSize)
    {
      if (samplesetlist[currentSamRow,1] == reflist[currentRefRow,1])
      {
        samplesetlist[currentSamRow,3] <- toString(reflist[currentRefRow,2])
      }
      currentRefRow = currentRefRow + 1
    }
    currentRefRow <- 1
    Sys.sleep(0.1)
    #cat("\r",currentSamRow)
    setTxtProgressBar(pb, currentSamRow)
    currentSamRow = currentSamRow + 1
  }

  samplesetlist <- data.frame(Gene_name= samplesetlist$genes, Cluster= factor(samplesetlist$cluster),Cell_name = samplesetlist$cellname ,stringsAsFactors = TRUE)

  clusterSize <- NROW(levels(samplesetlist$Cluster)) - 1
  currentCluster <- 0
  new.cluster.ids <- c()

  while (currentCluster <= clusterSize)
  {
    testinglist <- tail(names(sort(table(samplesetlist$Cell_name[samplesetlist$Cluster==currentCluster]))), 2)

    if (tail(names(sort(table(samplesetlist$Cell_name[samplesetlist$Cluster==currentCluster]))), 1) == unknown)
    {
      new.cluster.ids <- c(new.cluster.ids,testinglist[1])
    }
    else
    {
      new.cluster.ids <- c(new.cluster.ids,testinglist[2])
    }
    currentCluster = currentCluster + 1
  }

  if (nucIsoType == TRUE)
  {
    NumCelltypes <- NROW(new.cluster.ids)
    currentcell <- 1
    while (currentcell <= NumCelltypes)
    {
      if (new.cluster.ids[currentcell] == "Neuronal Stem Cells")
      {new.cluster.ids[currentcell] <- "Mitochondrial Cluster"}
      currentcell = currentcell + 1
    }

  }

  names(x = new.cluster.ids) <- levels(x = SampleSet)
  SampleSet <- RenameIdents(object = SampleSet, new.cluster.ids)
  DimPlot(object = SampleSet, reduction.use = "umap")
  return(SampleSet)

}
