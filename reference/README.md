#### Installation

```
conda create -n sfl_env python=2.7 pip -y
source activate sfl_env

SNAP_VERSION="snap-4.1.0-4.1-centos6.5-x64-py2.6"
wget https://snap.stanford.edu/snappy/release/$SNAP_VERSION.tar.gz
tar -xzvf $SNAP_VERSION.tar.gz
rm $SNAP_VERSION.tar.gz
mv $SNAP_VERSION snap
cd snap
python setup.py install
cd ..

conda install -c cvxgrp cvxpy libgcc -y
pip install -r requirements.txt
```

#### Usage

See `./run.sh` for usage.

