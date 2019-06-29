library(testthat)
assignInNamespace("cedta.override", c(data.table:::cedta.override,"AFM"),"data.table")
assignInNamespace("cedta.override", c(data.table:::cedta.override,"AFMImageAnalyser"),"data.table")

source("./configuration.R")


test_that("generateReport", {

  # exportDirectory=tempdir()
  # #exportDirectory="P:/MAT_DOCUMENTS/AccPlatform/person_display_R/"
<<<<<<< HEAD
  # 
=======
>>>>>>> 20dadf6f59c168d54865c6b5b5e69cce916aa856
  # data("AFMImageOfRegularPeaks")
  # newAFMImage<-AFMImageOfRegularPeaks
  # newAFMImage@fullfilename<-paste(exportDirectory,"AFMImageOfRegularPeaks.txt",sep="/")
  # generateReport(newAFMImage)
  # 
  # #   data("AFMImageOfNormallyDistributedHeights")
  # #   newAFMImage<-AFMImageOfNormallyDistributedHeights
  # #   newAFMImage@fullfilename<-paste(exportDirectory,"AFMImageOfNormallyDistributedHeights.txt",sep="/")
  # #   generateReport(newAFMImage)
  # #
  # #   data("AFMImageOfOnePeak")
  # #   newAFMImage<-AFMImageOfOnePeak
  # #   newAFMImage@fullfilename<-paste(exportDirectory,"AFMImageOfOnePeak.txt",sep="/")
  # #   generateReport(newAFMImage)
  # #
  # data("AFMImageOfAluminiumInterface")
  # newAFMImage<-AFMImageOfAluminiumInterface
  # newAFMImage@fullfilename<-paste(exportDirectory,"AFMImageOfAluminiumInterface.txt",sep="/")
  # newAFMImage<-extractAFMImage(newAFMImage, 100, 100, 32)
  # #generateCheckReport(newAFMImage)
  # generateReport(newAFMImage)
})
