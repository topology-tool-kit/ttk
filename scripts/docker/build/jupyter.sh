#! /bin/bash
set -e

build_pkgs \
	git				\
	nodejs			\
	npm

runtime_pkgs \
	python3-pip				  	\
	python3-numpy				\
	python3-pillow		   	  	\
	python3-widgetsnbextension 	\
	jupyter-notebook

pip3 install git+https://github.com/Kitware/ipyparaview.git
jupyter nbextension enable --py --sys-prefix ipyparaview

cat <<EOF > /etc/jupyter/jupyter_notebook_config.py
c.NotebookApp.ip = '0.0.0.0'
c.NotebookApp.allow_remote_access = True
c.NotebookApp.open_browser = False
EOF
