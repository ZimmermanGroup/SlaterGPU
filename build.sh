set +x

module unload nvidia-sdk
module load nvidia-sdk

cmake --version

# bash -i
# sleep 10
cmake -B $SRC_DIR/../build -S "/home/joshkamm/ElectronicStructure/SlaterGPU/"
cmake --build $SRC_DIR/../build -j 16

