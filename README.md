Dependencies are described on my note book

Add to ~/.bashrc :
   # This will make the comp_chem_py modules available everywhere on your computer
   export PYTHONPATH=/${comp_chem_py}/external:$PYTHONPATH
   export PYTHONPATH=/${comp_chem_py}/src:$PYTHONPATH

   # This will make the comp_chem_py scripts modules available everywhere on your computer
   # TODO: either make a bin directory with simlink or add several directories to PATH:
   export PATH=$PATH:/${comp_chem_py}/bin

TODOs
   - Make sure all scripts are protected by an if __name__ == '__main__': command
   - document everything with sphinx
   - link to read the docs
   - add tests (maybe...)
   - consider making comp_chem_utils an external submodule
