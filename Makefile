pkg := $(shell sed -n 's/^Package: *//p' DESCRIPTION)
version := $(shell sed -n 's/^Version: *//p' DESCRIPTION)
tar := $(pkg)_$(version).tar.gz

.PHONY: document build check install

document:
	Rscript -e "devtools::document()"

build:
	R CMD build .

check: build
	R CMD check --as-cran $(tar)

install: build
	R CMD INSTALL $(tar)
