UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(OS),Windows_NT)
	OSFLAGS = -shared -fPIC
	GCC = x86_64-w64-mingw32-gcc.exe
	TEST_OUT = src/build/honesttest_windows.plugin
else
	ifeq ($(UNAME_S),Linux)
		GCC = gcc
		OSFLAGS = -shared -fPIC -DSYSTEM=OPUNIX
		TEST_OUT = honesttest_unix.plugin
	endif
	ifeq ($(UNAME_S),Darwin)
		GCC = clang
		OSFLAGS = -bundle -DSYSTEM=APPLEMAC
		ifeq ($(UNAME_M),x86_64)
			TEST_OUT = honesttest_macosx86_64.plugin
		endif
		ifeq ($(UNAME_M),arm64)
			TEST_OUT = honesttest_macosxarm64.plugin
		endif
	endif
endif
CFLAGS = -O3 -Wall $(OSFLAGS)

all: clean honesttest

honesttest: honesttest.c stplugin.c
	$(GCC) $(CFLAGS) $^ -o $(TEST_OUT)

.PHONY: clean
clean:
	rm -f $(TEST_OUT)
