
echo 'managing environment'

mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
touch $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
touch $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

echo 'export PREFIX=/Applications
export CSDHOME=$PREFIX/CCDC/CSD_2020
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$CONDA_PREFIX/lib/python3.7/site-packages/ccdc/_lib:$LD_LIBRARY_PATH
export GHECOM_EXE=$PREFIX/ghecom_latest/ghecom
' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

echo 'unset PREFIX
unset CSDHOME
unset LD_LIBRARY_PATH
unset GHECOM_EXE' > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

cat $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
cat $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
echo 'complete'
