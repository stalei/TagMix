EXECS=ExportForTagging
CC=gcc
all:${EXECS}

ExportGalaxies:ExportForTagging.c
	${CC} -o ExportForTagging ExportForTagging.c
	
Clean:
	rm -f ${EXECS}
