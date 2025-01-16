set +x

module unload nvidia-sdk
module load nvidia-sdk

cmake --version

cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -B $SRC_DIR/../build -S $RECIPE_DIR
cmake --build $SRC_DIR/../build -j 16

cp -r ../build/lib $PREFIX/lib
mkdir -p $PREFIX/bin
cp ../build/examples/sgpu.exe $PREFIX/bin

# touch WORKING
# while [ -f "WORKING" ]; do
#   sleep 1
# done
# echo "File 'WORKING' no longer exists! Continuing execution..."
