#!/bin/bash
# Initialise python virtual env

# venv
python3 -m venv venv 
source ./venv/bin/activate

# upgrade pip etc
pip install --upgrade pip wheel pytest

# https://setuptools.pypa.io/en/latest/userguide/quickstart.html
pip install --upgrade build

pip install altair altair_saver pandas 
