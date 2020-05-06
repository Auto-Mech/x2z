# run with: . /path/to/fake-install.sh
export THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export BUILD_DIR=$THIS_DIR/debug
export PATH=$BUILD_DIR/bin:$PATH
export PYTHONPATH=$BUILD_DIR/lib:$PYTHONPATH
