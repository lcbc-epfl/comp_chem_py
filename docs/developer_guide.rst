Developer guide
===============

For now this is just a collection of things to remember when implementing 
something new in the library.

When adding a new script
------------------------

* Put everything that is not a function or pure data under the following
  condition::
    
        if __name__=="__main__":

  this is to avoid problems when compiling the documentation.
    
* A soft link shoud be added in the bin directory.

* Add the script path ``bin/your_script`` to the ``scripts`` list in the
  ``setup.py`` file in the root directory.


