#! /bin/bash

SPECTRA_VERSION=0.8.1

set -e 

echo "### build Spectra ###"

build_pkgs \
	curl \
	ca-certificates

# Spectra is header only -> fetch and extract includes directly to final location
fetch_url https://codeload.github.com/yixuan/spectra/tar.gz/v0.8.1 \
    | tar zvx -C usr --strip-components 1 "spectra-${SPECTRA_VERSION}/include"

