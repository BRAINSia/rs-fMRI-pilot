# rs-fMRI-pilot
Code for resting-state fMRI pilot study

## Installation

    export CC=gcc-4.2
    export CXX=g++-4.2
    export FFLAGS=-ff2c
    pip install numpy
    pip install scipy

    export CC=clang
    export CXX=clang+
    unset FFLAGS
    pip install -d virtualenvs/rs-fMRI-pilot/tmp pygraphviz  # This will fail
    cd build
    tar -xvf ../tmp/pygraphviz-1.1.tar.gz
    vim pygraphviz-1.1/setup.py

Edit pygraphviz/setup.py like so:

--
    ...
    library_path='/usr/local/Cellar/graphviz/2.38.0/lib'
    include_path='/usr/local/Cellar/graphviz/2.38.0/include'
    ...
--

Finally:

    cd pygraphviz-1.1/
    python setup.py install
    unset CC CXX
    pip install traits networkx nibabel nose sphinx pygments jinja2 docutils docopt pydicom ipython


## Usage

To load environment:

    $ source /Shared/sinapse/progopt/setup_virtualenv
    $ workon rs-fMRI-pilot

To debug crash files:

    $ nipype_display_crash <CRASHFILE>

To unload environment:

    $ deactivate
