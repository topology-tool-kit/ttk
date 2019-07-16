#! /bin/bash

echo "### build Spectra ###"

build_pkgs \
	curl

# Spectra is header only -> fetch and extract includes directly to final location
fetch_url https://codeload.github.com/yixuan/spectra/tar.gz/v0.8.1 \
    | tar zx -C usr --strip-components 1 --include "*/include/*"

