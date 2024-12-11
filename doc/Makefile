TEX = gtsmanual
BIB = gtsmanual
PDF = gtsmanual

DEPS_COMMON = ${TEX}.tex
#ifdef BIB
#	DEPS = ${DEPS_COMMON} ${TEX}.bbl
#else
	DEPS = ${DEPS_COMMON}
#endif

all : ${PDF}.pdf

#${TEX}.bbl : ${BIB}.bib ${TEX}.aux
#	bibtex ${TEX}
#	latex ${TEX}
#
#${TEX}.aux : ${DEPS_COMMON}
#	latex $<

${TEX}.dvi : ${DEPS}
	latex $<

${PDF}.pdf : ${TEX}.dvi
	dvipdf $<
	cp -p ${PDF}.pdf ~/w3_html/


clean :
	rm -f ${PDF}.pdf ${TEX}.dvi ${TEX}.aux ${TEX}.bbl
