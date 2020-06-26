#!/bin/bash -x

# Make errors fatal
#set -e

# Only build the paper with Julia 1.0
if [ $TRAVIS_JULIA_VERSION == "1.0.0" ]
then

    # Display some info
    which latex
    which pdftex
    which tlmgr
    
    # Generate figure placeholders, or use the static versions?
    cd $TRAVIS_BUILD_DIR/tex/figures
    if [[ $TRAVIS_COMMIT_MESSAGE == *"--use-static-figures"* ]]; then
        cp static/* .
    else
        python figure_placeholders.py
    fi

    # Generate the Julia figures
    # https://github.com/JuliaPy/PyPlot.jl/issues/317#issuecomment-337348563
    echo "Generating julia figures..."
    cd $TRAVIS_BUILD_DIR/tex/figures/julia
    for f in *.jl; do
        echo "Running $f..."
        travis_wait 30 julia "$f" || echo "ERROR: failed to run $f."
    done

    # Generate the Python figures
    echo "Generating python figures..."
    cd $TRAVIS_BUILD_DIR/tex/figures/python
    for f in *.py; do
        echo "Running $f..."
        travis_wait 30 python "$f" || echo "ERROR: failed to run $f."
    done

    # Build the paper with tectonic
    cd $TRAVIS_BUILD_DIR/tex
	tectonic refining_trappist1_ttvs.tex

    # Get the branch name
    if [[ $TRAVIS_COMMIT_MESSAGE == *"--branch="* ]]; then
        PDFBRANCH="$(echo $TRAVIS_COMMIT_MESSAGE | sed -e 's/.*--branch=\([[:alnum:]_-]*\).*/\1/')"
    else
        PDFBRANCH=$TRAVIS_BRANCH-pdf
    fi

    # Force push the paper to GitHub
    cd $HOME
    mkdir tmp && cd tmp
    git init
    git checkout --orphan $PDFBRANCH
    mkdir tex
    cp $TRAVIS_BUILD_DIR/tex/refining_trappist1_ttvs.pdf tex/
    git add -f tex/refining_trappist1_ttvs.pdf

    # Include figures in the commit?
    if [[ $TRAVIS_COMMIT_MESSAGE == *"--keep-figures"* ]]; then
        mkdir tex/figures
        mkdir tex/figures/julia
        cp $TRAVIS_BUILD_DIR/tex/figures/julia/*.pdf tex/figures/julia/
        git add -f tex/figures/julia/*.pdf
        mkdir tex/figures/python
        cp $TRAVIS_BUILD_DIR/tex/figures/python/*.pdf tex/figures/python/
        git add -f tex/figures/python/*.pdf
    fi

    git -c user.name='travis' -c user.email='travis' commit -m "building the paper"
    git push -q -f https://$GITHUB_USER:$GITHUB_API_KEY@github.com/$TRAVIS_REPO_SLUG $PDFBRANCH >/dev/null 2>&1

fi
