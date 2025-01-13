set +x

module unload nvidia-sdk
module load nvidia-sdk

cmake --version

# bash -i
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -B $SRC_DIR/../build -S "/home/joshkamm/ElectronicStructure/SlaterGPU/"
cmake --build $SRC_DIR/../build -j 16


sleep 90
