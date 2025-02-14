pixi global install cmake
pixi global install environment-modules
pixi global install helix
source ~/.pixi/envs/environment-modules/init/bash

wget --progress=bar:force:noscroll https://developer.download.nvidia.com/hpc-sdk/25.1/nvhpc_2025_251_Linux_x86_64_cuda_12.6.tar.gz
tar xpzf nvhpc_2025_251_Linux_x86_64_cuda_12.6.tar.gz
export NVHPC_SILENT="true"
export NVHPC_INSTALL_DIR="/home/runner/work"
export NVHPC_INSTALL_TYPE="auto"

# curl -fsSL https://pixi.sh/install.sh | bash
# source ~/.bashrc

nvhpc_2025_251_Linux_x86_64_cuda_12.6/install

module load ~/work/modulefiles/nvhpc/25.1
