#! /bin/bash

build_pkgs \
	curl

runtime_pkgs \
	libtbb2

# install OSPRay

curl -kL https://github.com/ospray/OSPRay/releases/download/v1.8.5/ospray-1.8.5.x86_64.linux.tar.gz | \
  tar zvx -C /usr --strip-components 1 \
  	  --exclude="scripts" \
	  --exclude="libtbb*"
