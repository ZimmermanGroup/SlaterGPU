pixi global install cmake
pixi global install environment-modules
pixi global install helix

# Stream download directly into tar extraction
echo "Downloading and extracting NVHPC simultaneously..."
wget -O - --progress=bar:force:noscroll https://developer.download.nvidia.com/hpc-sdk/25.1/nvhpc_2025_251_Linux_x86_64_cuda_12.6.tar.gz | tar xpz

export NVHPC_SILENT="true"
export NVHPC_INSTALL_DIR="/home/runner/work"
export NVHPC_INSTALL_TYPE="auto"

# curl -fsSL https://pixi.sh/install.sh | bash
# source ~/.bashrc

nvhpc_2025_251_Linux_x86_64_cuda_12.6/install
