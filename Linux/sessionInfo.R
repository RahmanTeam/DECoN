.libPaths("packrat/packrat_source")
install.packages("packrat/packrat_source/packrat_0.4.4.tar.gz",lib="packrat/packrat_source/",repos=NULL,type="source")

packrat::on()

packrat::restore()

print(sessionInfo())
