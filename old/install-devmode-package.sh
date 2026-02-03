#/bin/bash
# https://packaging.python.org/guides/distributing-packages-using-setuptools/#working-in-development-mode
pushd .
cd ./package-project/src/
pip3 install -e ./
popd

echo; echo "Showing package info:"
#echo 'export RINEX_PATH="$HOME/Documents/spaceweather/tec/CSN_TEC/"' >> $HOME/.bashrc 
#source ~/.bashrc

pip3 show pytec

python3 -mpip install python-decouple
python3 ./package-consumer-project/consumepackage.py
source ~/.bashrc
