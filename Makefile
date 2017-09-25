objects := $(wildcard R/*.R) DESCRIPTION
dir := $(shell pwd)
version := $(shell grep "Version" DESCRIPTION | sed "s/Version: //")
pkg := $(shell grep "Package" DESCRIPTION | sed "s/Package: //")
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log

.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

$(tar): $(objects)
	Rscript -e "library(methods); devtools::document();";
	R CMD build $(dir)

$(checkLog): $(tar)
	R CMD check --as-cran $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateDate
updateDate:
	dt=$$(date +"%Y-%m-%d");\
	sed -i "s/Date: [0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}/Date: $$dt/" DESCRIPTION;

.PHONY: clean
clean:
	rm -rf *~ */*~ *.Rhistroy *.tar.gz *.Rcheck/ .\#*
