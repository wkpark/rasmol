CC = gcc
CFLAGS = -g -O2

prepdoc:	prepdoc.c
		$(CC) $(CFLAGS) -o prepdoc prepdoc.c

rasmol.txt:	prepdoc rasmol.doc
		rm -f rasmol.txt
		prepdoc -text rasmol.doc > rasmol.txt

rasmol.html:	prepdoc rasmol.doc
		rm -f rasmol.html
		prepdoc -html rasmol.doc > rasmol.html

rasmol.hlp:	prepdoc rasmol.doc
		rm -f rasmol.hlp
		prepdoc -help rasmol.doc > rasmol.hlp

raswin.rtf:	prepdoc rasmol.doc
		rm -f raswin.rtf
		prepdoc -rtf rasmol.doc > raswin.rtf

rasmol.man:	prepdoc rasmol.doc
		rm -f rasmol.man
		prepdoc -troff rasmol.doc > rasmol.man
