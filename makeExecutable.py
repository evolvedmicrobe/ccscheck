import os
# This command may need to run in the terminal ahead of time
# You may also need to add -liconv to the CC command and run it in the terminal if you 
# see something along the lines of:
# undefined reference to `libiconv_open'
## SEE


## MORE NOTES:
## Mono knows what directory it was originally installed in, as it is a string in the original library.  This can cause 
## problems in the re-distributed mono if the configuration file or the other file go missing.  To avoid this, I set
## -config-dir equal to a made up name, and check that MONO_CFG_DIR is null and that MONO_PATH are null on boot.
## Either of these can goof the runtime so that it fails, so we need to ensure there are no conflicts.
## Another option is to edit the temp.c file emitted by mkbundle, os that mono_set_dirs does not point someplace funny.


from sys import platform as _platform


if _platform == "linux" or _platform == "linux2":
	os.system("export PKG_CONFIG_PATH=$HOME/mono64/lib/pkgconfig/")
        # note the dummy config-dir /nothing, just want to prevent a search on main path.
	cmd = "mkbundle --keeptemp --static  --deps --config-dir /nothing --config /home/UNIXHOME/ndelaney/mono64/etc/mono/config -o ccscheck ccscheck.exe Bio.BWA.dll Bio.Core.dll Bio.Desktop.dll Bio.Platform.Helpers.dll"
	os.system(cmd)
	pass
    # linux
elif _platform == "darwin":
	os.system("export PKG_CONFIG_PATH=/Users/nigel/mono64_4.0/lib/pkgconfig/")
	# This will fail due to missing flags
	cmd = "mkbundle --config-dir /nothing --static  --deps -o ccscheck ccscheck.exe Bio.BWA.dll Bio.Core.dll Bio.Desktop.dll Bio.Platform.Helpers.dll"
	os.system(cmd)
	# so we add the flags back in and remake it, like all that corefoundation stuff
	cmd = "cc -o ccscheck  -framework CoreFoundation -lobjc -liconv -Wall `pkg-config --cflags mono-2` temp.c  `pkg-config --libs-only-L mono-2` `pkg-config --variable=libdir mono-2`/libmono-2.0.a `pkg-config --libs-only-l mono-2 | sed -e \"s/\-lmono-2.0 //\"` temp.o"
	print cmd
	os.system(cmd)
elif _platform == "win32":
    # Windows...
    pass


### OLD CODE BELOW, this generates a command to mkbundle, but it needs some manual futzing

#os.system("setenv PKG_CONFIG_PATH /Users/nigel/mono64_4.0/lib/pkgconfig/")
#t=[x for x in os.listdir(os.getcwd()) if x.endswith((".dll",".exe"))]
#cmd="mkbundle --static  --deps -o ccscheck "+" ".join(t)
#above didn't work, order seems to matter
#cmd="mkbundle --deps -o curvefitter --static CurveFitterMonoGUI.exe  ShoOptimizer.dll Microsoft.Solver.Foundation.dll MatrixInterf.dll ShoArray.dll alglibnet2.dll MatrixArrayPlot.dll ZedGraph.dll GrowthCurveLibrary.dll "

# -framework CoreFoundation -lobjc -liconv" <- also needed on mac for CC command
#print cmd
#os.system(cmd)

