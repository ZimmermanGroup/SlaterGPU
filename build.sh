set +x

cmake --version

cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -B $SRC_DIR/../build -S $RECIPE_DIR
cmake --build $SRC_DIR/../build -j 16
cmake --install $SRC_DIR/../build

# touch WORKING
# while [ -f "WORKING" ]; do
#   sleep 1
# done
# echo "File 'WORKING' no longer exists! Continuing execution..."
