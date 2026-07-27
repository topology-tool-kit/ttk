OIDN_VERSION=1.4.2

# only build on x86_64
if [ $(arch) == "x86_64" ]; then

    require-pkgs \
        python3-minimal

    fetch-src https://github.com/OpenImageDenoise/oidn/releases/download/v${OIDN_VERSION}/oidn-${OIDN_VERSION}.src.tar.gz

    cmake-default \
        -DOIDN_APPS=OFF

fi
