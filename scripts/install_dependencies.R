# Run explicitly; the analysis never installs packages automatically.
install.packages(c("rTensor","readxl","jsonlite","curl","RhpcBLASctl"),
  repos="https://cloud.r-project.org",Ncpus=1L)
