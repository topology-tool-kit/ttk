#!/bin/bash

aptget="apt-get -qq -o=Dpkg::Use-Pty=0 --no-install-recommends"

function build_pkgs() {
    # compute which packages must be newly installed
    # if we simply install existing packages, they might get
    # erroneously removed later
    dpkg-query -f '${db:Status-Abbrev} ${Package}\n' -W \
        | awk '/^ii/ {print $2}' \
        > /tmp/pkgs1

    for pkg in "$@"; 
        do echo ${pkg} >> /tmp/pkgs2; 
    done

    pkgs=$(sort /tmp/pkgs2 | uniq | sort - /tmp/pkgs1 /tmp/pkgs1 | uniq -u)

    echo "installing build packages: ${pkgs}"
    ${aptget} install ${pkgs}
    _build_pkgs="${_build_pkgs} ${pkgs}"

}

function runtime_pkgs() {
    echo "installing runtime packages"
    ${aptget} install $@
}

function fetch_url() {
    curl -qL "$1"
}

configure_args=""

function conf_args() {
    configure_args+=" $@"
}

# build basename
what=$(basename $1 .sh)
echo "---- begin $what ----"

# update/initialize the apt cache 
${aptget} update

# call build script
_build_dir="/tmp/${what}"

mkdir -p ${_build_dir}
pushd ${_build_dir}
source $1
popd

# remove build dir and build script
echo "removing build dir"
rm -rf ${_build_dir} $1

# remove build-only-packages
echo "removing build packages"
${aptget} purge ${_build_pkgs}

# also remove suggested / recommended packages that were installed as dependencies od build-only packages
${aptget} autoremove --purge -o APT::Autoremove::RecommendsImportant=0 -o APT::Autoremove::SuggestsImportant=0

# clear / remove apt cache
${aptget} clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# remove script

echo "---- end $what ----"
