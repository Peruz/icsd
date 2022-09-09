#!/bin/sh
# this file needs to be run from inside the docsrc/ directory

# make sure to activate the pyenv environement with the necessary packages before running the notebooks
#source ../../pyenv/bin/activate

# go to docsrc/ and run make html
make html

# then move the build html pages to the docs/ folder
rm -r ../docs
mkdir ../docs
touch ../docs/.nojekyll # create empty file

cp -r _build/html/* ../docs/

# then push the output of docs/* to the remote branch


