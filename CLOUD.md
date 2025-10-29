CLOUD 
# ssh key
ssh -i ~/.ssh/id_rsa cponti@104.198.208.95

sudo apt-get update
sudo apt-get install -y google-cloud-sdk

gsutil -m cp -r gs://my-flow-cytometry-bucket/my_cloud_workflow ~/
gsutil -m cp -r gs://my-flow-cytometry-bucket/my_project_Data ~/

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

python3.11 -m venv ~/my_cloud_workflow/venv
source ~/my_cloud_workflow/venv/bin/activate
pip install --upgrade pip setuptools wheel
pip install --no-build-isolation --only-binary=:all: pandas anndata scanpy fcsparser flowsom pytometry  matplotlib pooch leidenalg


python ~/my_cloud_workflow/main.py

scp -i ~/.ssh/id_rsa data_paths.py cponti@136.112.133.20:/home/cponti/test/

# mounting
sudo mkdir -p /mnt/gcs
sudo chown $USER:$USER /mnt/gcs
ls -ld /mnt/gcs
gcsfuse --implicit-dirs my-flow-cytometry-bucket /mnt/gcs
ls /mnt/gcs

python3 /mnt/gcs/my_cloud_workflow/main.py


# increase computation power 32 - e2-medium (2 vCPUs, 4 GB Memory)





# install gcsfuse 
# Remove the bad repo file first
sudo rm /etc/apt/sources.list.d/gcsfuse.list

# Add the last supported repo (noble)
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt gcsfuse-noble main" | sudo tee /etc/apt/sources.list.d/gcsfuse.list

# Add Google’s signing key
curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg

# Update and install gcsfuse
sudo apt update
sudo apt install gcsfuse

# GUI 
sudo apt update && sudo apt upgrade -y
sudo apt install -y wget curl nano
sudo apt install -y xfce4 xfce4-goodies
sudo apt install -y xrdp
sudo systemctl enable xrdp
sudo systemctl start xrdp
echo xfce4-session > ~/.xsession
sudo sed -i.bak '/fi/a startxfce4' /etc/xrdp/startwm.sh
sudo systemctl restart xrdp
sudo adduser XXXX # where XXXX is the username you wish to use, best to input a password here.
sudo usermod -aG sudo XXXX
sudo -i -u XXXX bash
echo xfce4-session > ~/.xsession
exit

# install pyenv 
curl -fsSL https://pyenv.run | bash
sudo apt install vim
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init - bash)"' >> ~/.bashrc
exec "$SHELL"
sudo apt install sqlite3
sudo apt install lzma-dev
sudo apt install liblzma-dev
pyenv install 3.12 
# creates /home/cponti/.pyenv/shims/python3.12 
which python3.12  
ls /usr/bin/ | grep python
/home/cponti/.pyenv/shims/python3.12 -m venv ~/my_cloud_workflow/venv2
pyenv venv 3.12.12 ~/my_cloud_workflow/venv2
pyenv global 3.12.12
python -m venv ~/my_cloud_workflow/venv2
pip install --no-build-isolation --only-binary=:all: pandas anndata scanpy fcsparser flowsom pytometry matplotlib pooch