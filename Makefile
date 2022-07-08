# ---------------------------------------------------------------------
# OS parsing

# README.md
# honestdid.pkg
# stata.toc
# src/ado/honestdid.ado
# src/ado/honestdid.sthlp
# src/plugin/honestosqp.c
# src/plugin/honestosqp.h
# src/plugin/honestecos.c
# src/plugin/honestecos.h

ifeq ($(OS),Windows_NT)
	OSFLAGS = -shared -fPIC
	GCC = x86_64-w64-mingw32-gcc.exe
	OSQP_OUT = src/build/honestosqp_windows.plugin
	ECOS_OUT = src/build/honestecos_windows.plugin
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OSFLAGS = -shared -fPIC -DSYSTEM=OPUNIX
		OSQP_OUT = src/build/honestosqp_unix.plugin
		ECOS_OUT = src/build/honestecos_unix.plugin
	endif
	ifeq ($(UNAME_S),Darwin)
		OSFLAGS = -bundle -DSYSTEM=APPLEMAC
		OSQP_OUT = src/build/honestosqp_macosx.plugin
		ECOS_OUT = src/build/honestecos_macosx.plugin
	endif
	GCC = gcc
endif

ifeq ($(EXECUTION),windows)
	OSFLAGS = -shared
	GCC = x86_64-w64-mingw32-gcc
	OSQP_OUT = src/build/honestosqp_windows.plugin
	ECOS_OUT = src/build/honestecos_windows.plugin
endif

OSQP_H = /home/mauricio/bulk/lib/osqp/include
OSQP_A = /home/mauricio/bulk/lib/osqp/build/out/libosqp.a
ECOS_H = /home/mauricio/bulk/lib/ecos/include
ECOS_A = /home/mauricio/bulk/lib/ecos/libecos.a /home/mauricio/bulk/lib/ecos/libecos_bb.a
CFLAGS = -Wall -O3 $(OSFLAGS)

# ---------------------------------------------------------------------
# Main

## Compile directory
all: clean honestosqp honestecos

# ---------------------------------------------------------------------
# Rules

## Compile OSQP plugin
honestosqp: src/plugin/honestosqp.c src/plugin/stplugin.c
	$(GCC) $(CFLAGS) -I$(OSQP_H) -o $(OSQP_OUT) $^ $(OSQP_A)

## Compile ECOS plugin
honestecos: src/plugin/honestecos.c src/plugin/stplugin.c
	$(GCC) $(CFLAGS) -DLDL_LONG -DDLONG -I$(ECOS_H) -o $(ECOS_OUT) $^ $(ECOS_A)

.PHONY: clean
clean:
	rm -f $(OSQP_OUT) $(ECOS_OUT)

#######################################################################
#                                                                     #
#                    Self-Documenting Foo (Ignore)                    #
#                                                                     #
#######################################################################

.DEFAULT_GOAL := show-help

.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')

