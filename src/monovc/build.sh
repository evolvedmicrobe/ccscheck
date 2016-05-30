#!/bin/bash
# Check for the platform first
if [[ "$OSTYPE" == "linux-gnu" ]]; then
        # ubuntu
        export MONO_ROOT=$HOME/mono64/       

elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
        export MONO_ROOT=/Users/nigel/mono64_4.2.1/       
else
	echo "You can't build the release versions here!"
        # Unknown.
fi

rm -rf build
mkdir build

xbuild /p:Configuration=Release ../ccscheck.sln
cp ../EmbeddedCCSCheck/bin/Release/* build/
cp ../../lib/libbwacsharp.so build/

cd build
## PACKAGE INTO EXECUTABLE and SHARED LIBRARY
export PKG_CONFIG_PATH=$MONO_ROOT/lib/pkgconfig/

# FIRST SHARED LIB 
# the C bit
mkbundle -c -o mono_embed_host.c -oo mono_embed_libs.o --static --keeptemp --deps --nomain EmbeddedCCSCheck.exe ccscheck.exe Bio.BWA.dll Bio.Core.dll Bio.Desktop.dll Bio.Platform.Helpers.dll 
clang -g -c -Wall `pkg-config --cflags mono-2` -o temp_hold.o mono_embed_libs.o mono_embed_host.c
ld -r temp_hold.o mono_embed_libs.o -o embeded_libs.o
# the C++ bits
clang++ -dynamiclib -o libVariantCaller.so -D PROGRAM -I../include -framework CoreFoundation -lobjc -liconv -Wall `pkg-config --cflags mono-2`  `pkg-config --libs-only-L mono-2` `pkg-config --variable=libdir mono-2`/libmono-2.0.a `pkg-config --libs-only-l mono-2 | sed -e "s/\-lmono-2.0 //"` embeded_libs.o ../VariantCaller.cpp
# NOW EXECUTABLE
clang++ -g -o variantcaller -D PROGRAM -I../include -I../ -framework CoreFoundation -lobjc -liconv -Wall `pkg-config --cflags mono-2`  `pkg-config --libs-only-L mono-2` `pkg-config --variable=libdir mono-2`/libmono-2.0.a `pkg-config --libs-only-l mono-2 | sed -e "s/\-lmono-2.0 //"` embeded_libs.o ../VariantCaller.cpp 
# Clean and test
rm *.dll
rm *.exe
rm Bio.*
rm *.mdb
rm *.o
rm *.s
./variantcaller
cd ../ 



