ISPC_VERSION=1.16.1


require-pkgs \
    clang        \
    libclang-dev \
    llvm-dev     \
    $([ "$(arch)" == "x86_64"  ] && echo "libc6-dev-i386-cross" ) \
    $([ "$(arch)" == "aarch64" ] && echo "libc6-dev-armhf-cross") \
    m4              \
    bison           \
    flex            \
    libz-dev


fetch-src https://github.com/ispc/ispc/archive/refs/tags/v${ISPC_VERSION}.tar.gz

conf_args \
    -DCMAKE_BUILD_TYPE=Release  \
    -DX86_ENABLED=$([ "$(arch)" == "x86_64"  ] && echo "ON" || echo "OFF") \
    -DARM_ENABLED=$([ "$(arch)" == "aarch64" ] && echo "ON" || echo "OFF") \
    -DWASM_ENABLED=OFF          \
    -DISPC_NO_DUMPS=ON          \
    -DISPC_STATIC_LINK=ON       \
    -DISPC_INCLUDE_EXAMPLES=OFF \
    -DISPC_INCLUDE_TESTS=OFF    \
    -DISPC_INCLUDE_UTILS=OFF

cmake-default