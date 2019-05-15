#!/usr/bin/env python
"""
Load orbital ``.cube`` files in VMD and exports them to ``.png`` format.

Usage example::

    python vmd_plot_cube.py WAVEFUNCTION.1.pdb -c WAVEFUNCTION.*.cube


Requirements:

    * The VMD package should be available by just writting vmd in the terminal.
    * The conversion from .tga to .png requires the Imagemagick convert utility.
    * ...

This script has been written by learning and stilling from other scripts:

    * by Felix Plasser (http://www.chemical-quantum-images.blogspot.de)
    * by Jan-Michael Mewes (http://workadayqc.blogspot.de)

"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import os
import argparse

# ----------------------------------------------------------------
# THIS IS THE TEMPLATE FOR THE VMD LOAD SCRIPT
vmd_setup_graphic = """
# ==========================================
# SETUP GRAPHIC ENVIRONMENT
#
# General settings
# ------------------------------------------
axes location Off
display projection Orthographic
display rendermode GLSL
display depthcue off
color Display Background white
#
# PBC
#pbc set {13 13 13} -all
#pbc box
#pbc wrap -compound residue
#
menu graphics on
mol modstyle 0 0 CPK
#
# Orbital settings
mol addrep 0
mol addrep 0
# ------------------------------------------
# display positive side in blue
mol modmaterial 1 0 Translucent
mol modstyle 1 0 Isosurface  0.01 0 0 0 1 1
mol modcolor 1 0 ColorID 0
# ------------------------------------------
# display negative side in red
mol modmaterial 2 0 Translucent
mol modstyle 2 0 Isosurface -0.01 0 0 0 1 1
mol modcolor 2 0 ColorID 1
# ------------------------------------------
"""

load_script_name = "load_orbs_cubes.vmd"
prt_script_name = "print_orbs.vmd"
conv_script_name = "convert_tga_png.sh"

def create_scripts(cube_files):

    # ----------------------------------------------------------------
    # 1) CREATE VMD SCRIPT THAT WILL LOAD THE CUBE FILES IN VMD
    with open(load_script_name,"w") as script:

        script.write( vmd_setup_graphic )

        # loop over all cube files
        for cube in cube_files:

            script.write( "mol addfile {}\n".format(cube) )

    # ----------------------------------------------------------------
    # 2) CREATE VMD SCRIPT THAT WILL EXPORT PICTURES OF THE ORBITALS
    with open(prt_script_name,"w") as script:

        # loop over all cube files
        for idx,cube in enumerate(cube_files):

            script.write( "mol modstyle 1 0 Isosurface  0.01 {} 0 0 1 1\n".format(idx))
            script.write( "mol modstyle 2 0 Isosurface -0.01 {} 0 0 1 1\n".format(idx) )
            script.write( "render TachyonInternal {}.tga \n".format( cube[:-5] ) )

        script.write( "quit \n" )

    # ----------------------------------------------------------------
    # 3) CREATE CONVERSION SCRIPT 
    with open(conv_script_name,"w") as script:

        script.write( "#!/bin/bash \n")
        # loop over all cube files
        for idx,cube in enumerate(cube_files):
            script.write( "convert {0}.tga {0}.png \n".format( cube[:-5] ) )
            script.write( "rm {}.tga \n".format( cube[:-5] ) )


if __name__=="__main__":

    # ----------------------------------------------------------------
    # GET THE RELEVANT FILENAMES (STRUCTURE AND .CUBE)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xyz", type=str, help="File containing the structure of the system (.pdb, .xyz, ...)")
    parser.add_argument('-c','--cube', nargs='+', help='List of cube files to render', required=True)
    args = parser.parse_args()

    create_scripts(args.cube)

    # ----------------------------------------------------------------
    # ASK THE USER WHAT TO DO
    job = '' #raw_input("""
#    GENERATE SCRIPT AND SHOW COMMANDS : 1 [default]
#
#    DISPLAY CUBE FILES IN VMD DIRECTLY: 2
#    """)

    if job=='2':

        # ----------------------------------------------------------------
        # PRINT SOME INFORMATION TO THE USER
        raw_input("""
        We will now launch VMD with the orbital cube files.
        In the graphic window, make your modifications until you are satisfied.
        To export images of the orbitals run the following command in the vmd terminal:
         
        vmd > play {}
         
        Type <ENTER> when you are ready.
        """.format( prt_script_name ) )
         
         
        # ----------------------------------------------------------------
        # LAUNCH VMD WITH LOAD SCRIPT
        os.system( "vmd {} -e {}".format(args.xyz, load_script_name) )
         
         
        # ----------------------------------------------------------------
        # CONVERT THE .TGA FILES TO PNG FORMAT
        print("Conversion of images from .tga to .png format...")
        os.system( "bash {}".format(conv_script_name) )
        #for cube in args.cube:
        #    os.system( "convert {0}.tga {0}.png".format( cube[:-5] ) )
        #    os.system( "rm {}.tga".format( cube[:-5] ) )

    else:
        print("""
        COMMANDS TO RUN:
            # 1) Loads cube files and display nicely using vmd:
            vmd {0} -e {1}

            # 2) Generates .tga images of cube files from vmd prompt:
            vmd > play {2}

            # 3) Quit vmd and convert .tga to .png format:
            bash {3}
        """.format(args.xyz, load_script_name, prt_script_name, conv_script_name))

