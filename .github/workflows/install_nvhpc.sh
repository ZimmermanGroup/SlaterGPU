pixi global install cmake
pixi global install environment-modules
pixi global install helix

# Set environment variables
export NVHPC_SILENT="true"
export NVHPC_INSTALL_DIR="/home/runner/work"
export NVHPC_INSTALL_TYPE="auto"

# Check if we have a cache hit
if [ "$CACHE_HIT" == "true" ]; then
  echo "Using cached NVHPC installation"
  # Ensure modulefiles directory exists
  mkdir -p ~/work/modulefiles
else
  echo "No cache found. Downloading and installing NVHPC..."
  
  # Stream download directly into tar extraction
  wget -O - --progress=bar:force:noscroll https://developer.download.nvidia.com/hpc-sdk/25.1/nvhpc_2025_251_Linux_x86_64_cuda_12.6.tar.gz | tar xpz

  # Install NVHPC
  nvhpc_2025_251_Linux_x86_64_cuda_12.6/install
  
  # Clean up installation files to reduce cache size
  rm -rf nvhpc_2025_251_Linux_x86_64_cuda_12.6
  
  # Trim cache size by removing unnecessary files
  find $NVHPC_INSTALL_DIR/Linux_x86_64 -name "*.a" -delete || true
  find $NVHPC_INSTALL_DIR/Linux_x86_64 -name "examples" -type d -exec rm -rf {} + 2>/dev/null || true
  find $NVHPC_INSTALL_DIR/Linux_x86_64 -name "doc" -type d -exec rm -rf {} + 2>/dev/null || true
  find $NVHPC_INSTALL_DIR/Linux_x86_64 -name "samples" -type d -exec rm -rf {} + 2>/dev/null || true
fi
