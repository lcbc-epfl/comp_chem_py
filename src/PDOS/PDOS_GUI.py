#!/usr/bin/env python


__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

try:
    # python 2
    from Tkinter import *
except ImportError:
    # python 3
    from tkinter import *

from PDOS.modules.PDOS_module import calculate_and_plot_pdos


def try_plot_from_gui():
    """Check if inputs from GUI can be used to plot the PDOS, if yes then just do it..."""
    try:
        calculate_and_plot_pdos(gaussian_file.get(), group_file.get(), npts.get(), 
                fwhm.get(), xmin.get(), xmax.get(), plot_dos.get(), mulliken.get())
    except:
        error = Toplevel()
        error.title("*-error-*")
        Label(error, text = "An Input Parameter Is Not Suitable !!", foreground="red").pack()       
        Label(error, text = "Check Your Input Files !!", foreground="red").pack()    
        Label(error, text = "Be Carefull: The two first data have to be the corect names of the files,", foreground="red").pack()
        Label(error, text = "And Npts must receive an integer, FWMH must receive a  float !!", foreground="red").pack()  
        Label(error, text = "If you want, you can just enter the names of the input files and PLOT.").pack()


# for compatibility with autodoc in sphinx
if __name__ == '__main__':

    # create a window for the display
    window = Tk()
    window.title("*-PDOS-*")
    window.geometry("500x310")
    
    
    # Create the entry for the Gaussian file
    Label(window, text = "Name of the Gaussian Output File:").place(x=30,y=15)
    gaussian_file = StringVar()
    Entry(window, textvariable = gaussian_file, width=20).place(x=300,y=15)
    
    
    # Create the entry for the File with the Groupe of AOs
    Label(window, text = "Name of the file with the group of AOs:").place(x=30,y=40)
    group_file= StringVar()
    Entry(window, textvariable = group_file, width=20).place(x=300,y=40)
    
    
    # Create the entry for the npts
    Label(window, text = "Npts:").place(x=30,y=80)
    npts = IntVar(value = 0)
    Entry(window, textvariable = npts, width=10).place(x=100,y=80)
    
    
    # Create the entry for the FWHM
    Label(window, text = "FWHM:").place(x=210,y=80)
    fwhm = DoubleVar(value = 1.0)
    Entry(window, textvariable = fwhm, width=10).place(x=270,y=80)
    
    
    # Ask if we are gonna plot DOS
    Label(window, text = "Plot DOS ?").place(x=50,y=110)
    plot_dos = BooleanVar(value = True)
    Radiobutton(window, text = "Yes", variable = plot_dos, value = True).place(x=130,y=110)
    Radiobutton(window, text = "No", variable = plot_dos, value = False).place(x=130,y=130)
    
    
    # Create the entry for the interval of energy
    ypos = 160
    Label(window, text = "Energy from : ").place(x=30,y=ypos)
    xmin = DoubleVar(value = 0.0)
    Entry(window, textvariable = xmin, width=10).place(x=115,y=ypos)
    
    Label(window, text = "To : ").place(x=220,y=ypos)
    xmax = DoubleVar(value = 0.0)
    Entry(window, textvariable = xmax, width=10).place(x=270,y=ypos)
    
    
    # Choose between Mulliken and Lowdin partial charges
    Label(window, text = "Partial charges:").place(x=50,y=200)
    mulliken = BooleanVar(value = False)
    Radiobutton(window, text = "Lowdin", variable = mulliken, value = False).place(x=160,y=200)
    Radiobutton(window, text = "Mulliken", variable = mulliken, value = True).place(x=160,y=220)
    
    # Run main calculationg
    Button(window, text = "PLOT", command = try_plot_from_gui, foreground="blue").place(x=220,y=260)
    
    window.mainloop()
