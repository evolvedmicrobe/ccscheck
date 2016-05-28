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
export PKG_CONFIG_PATH=$MONO_ROOT/lib/pkgconfig/
clang++ -o variantcaller -framework CoreFoundation -lobjc -liconv -Wall `pkg-config --cflags mono-2`  `pkg-config --libs-only-L mono-2` `pkg-config --variable=libdir mono-2`/libmono-2.0.a `pkg-config --libs-only-l mono-2 | sed -e "s/\-lmono-2.0 //"` VariantCaller.cpp
./VariantCaller
