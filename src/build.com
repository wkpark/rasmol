$! BUILD.COM
$! VMS Compilation Script
$! RasMol v2.3
$!
$ If F$Trnlnm ("X11") .eqs. "" Then Define X11 DECW$Include
$ If F$Trnlnm ("SYS") .eqs. "" Then Define Sys Sys$Share
$!
$ cc/optimize rasmol.c
$ cc/optimize molecule.c
$ cc/optimize transfor.c
$ cc/optimize command.c
$ cc/optimize abstree.c
$ cc/optimize render.c
$ cc/optimize x11win.c
$ cc/optimize pixutils.c
$ cc/optimize outfile.c
$!
$ link /exec=rasmol rasmol.obj, molecule.obj, transfor.obj, command.obj, -
    abstree.obj, render.obj, x11win.obj, pixutils.obj, outfile.obj, rasmol/opt

