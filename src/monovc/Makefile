.DEFAULT_GOAL := static_lib
# Two variables which can be overwritten by the user
MONO_ROOT = $(shell dirname `which mono`)/../
BOOST_LOCATION = /usr/local/include


# Set these environmental variables
export CPLUS_INCLUDE_PATH=$(BOOST_LOCATION):$CPLUS_INCLUDE_PATH
export PKG_CONFIG_PATH=$(MONO_ROOT)/lib/pkgconfig/

static_lib: .csharp
	mkdir -p build
	cp ../EmbeddedCCSCheck/bin/Release/* build/
	cp ../../lib/libbwacsharp.so build/
	# Make C wrapper and mono_embed_libs.o - contains all the dll files in the data section of an object file
	cd build && mkbundle -c -o mono_embed_host.c -oo mono_embed_libs.o --static --keeptemp --deps --nomain EmbeddedCCSCheck.exe ccscheck.exe Bio.BWA.dll Bio.Core.dll Bio.Desktop.dll Bio.Platform.Helpers.dll 
	# Now compile with the wrapper 
	cd build && $(CC) -g -c -Wall `pkg-config --cflags mono-2` -o temp_hold.o mono_embed_libs.o mono_embed_host.c
	# Compile the Variant Caller
	cd build && $(CXX) --std=c++11 -g -c -I../include -I/Users/nigel/git/ConsensusCore2/include -Wall `pkg-config --cflags mono-2` ../VariantCaller.cpp -o VariantCaller.o
	# Now link them to create our static library
	cd build && ld -r `pkg-config --variable=libdir mono-2`/libmono-2.0.a temp_hold.o mono_embed_libs.o VariantCaller.o -o VariantCaller.a

.csharp:
	xbuild /p:Configuration=Release ../ccscheck.sln

clean:
	rm ../EmbeddedCCSCheck/bin/Release/*
	rm -rf build


