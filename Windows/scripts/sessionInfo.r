.libPaths("packrat/packrat_source")
install.packages("packrat/packrat_source/packrat_0.4.4.tar.gz",lib="packrat/packrat_source",repos=NULL,type="source")

print("Packrat installed")
print("Turning packrat on")

packrat::on()

print("Restoring packrat packages")
packrat::restore()

print(sessionInfo())

file.copy("packrat/packrat_source/.Rprofile",".Rprofile",overwrite=TRUE)
