library(rsconnect)
#since we have raw data that's big, only deploy whats needed
deployApp(appFileManifest = "appFileManifestForDeployment.txt", forceUpdate=T)
