# Cloud Workflow Setup Guide

This describes how to set up a cloud-based Linux environment for running the my_cloud_workflow Python analysis pipeline, including GUI access, Python environment management, and Google Cloud Storage mounting.

## 1. System Preparation & GUI Installation
Update and install essential packages, XFCE desktop, and XRDP for remote desktop access:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y wget curl nano
sudo apt install -y xfce4 xfce4-goodies
sudo apt install -y xrdp
sudo systemctl enable xrdp
sudo systemctl start xrdp
echo xfce4-session > ~/.xsession
sudo sed -i.bak '/fi/a startxfce4' /etc/xrdp/startwm.sh
sudo systemctl restart xrdp
```

Create a new user (replace `XXXX` with your username):
```bash
sudo adduser XXXX
sudo usermod -aG sudo XXXX
sudo -i -u XXXX bash
echo xfce4-session > ~/.xsession
exit
```

## 2. Install pyenv and Python
Install pyenv for managing Python versions:
```bash
curl -fsSL https://pyenv.run | bash
sudo apt install vim sqlite3 lzma-dev liblzma-dev
# Add pyenv to bashrc
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init - bash)"' >> ~/.bashrc
exec "$SHELL"
```

Install Python 3.12 and create a virtual environment:
```bash
pyenv install 3.12.12
pyenv global 3.12.12
python -m venv ~/my_cloud_workflow/venv2
```

## 3. Install Python Packages
Activate your environment and install required packages:
```bash
source ~/my_cloud_workflow/venv2/bin/activate
pip install --upgrade pip setuptools wheel
pip install --no-build-isolation --only-binary=:all: pandas anndata scanpy fcsparser flowsom pytometry matplotlib pooch leidenalg
```

## 4. (Optional) Install Miniconda
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

## 5. Mount Google Cloud Storage Bucket
```bash
sudo mkdir -p /mnt/gcs
sudo chown $USER:$USER /mnt/gcs
gcsfuse --implicit-dirs my-flow-cytometry-bucket /mnt/gcs
ls /mnt/gcs
```

## 6. Running the Workflow
Run the main analysis script (update the path as needed):
```bash
python3 /mnt/gcs/my_cloud_workflow/main.py
```

## Notes
- Replace `XXXX` with your desired username.
- For SSH access: `ssh -i ~/.ssh/id_rsa cponti@<your-server-ip>`
- Ensure all data and scripts are available in the correct directories (local or mounted GCS).
- For troubleshooting, check permissions and environment activation.

---
This is intended for Ubuntu/Debian-based cloud VMs. 