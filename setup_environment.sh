#echo 'create environment'
#conda env create -n hotspots -f environment.yml

echo 'managing environment'
conda activate hotspots

mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
touch $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
touch $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

echo 'export PREFIX=/local/pcurran
export CSDHOME=$PREFIX/CCDC/CSD_2021
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$CONDA_PREFIX/lib/python3.7/site-packages/ccdc/_lib:$LD_LIBRARY_PATH
export GHECOM_EXE=$PREFIX/ghecom/src/ghecom
' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

echo 'unset PREFIX
unset CSDHOME
unset LD_LIBRARY_PATH
unset GHECOM_EXE' > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

cat $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
cat $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
echo 'complete'


