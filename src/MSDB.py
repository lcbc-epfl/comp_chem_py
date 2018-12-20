#!/usr/bin/env python
"""Graphical User Interface for the molecular structure database"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

try:
    # python 2
    import Tkinter as tk
    import ttk
except ImportError:
    # python 3
    import tkinter as tk
    from tkinter import ttk

import tkFont
import tkFileDialog as tkfd
import os
import MySQLdb as db

from comp_chem_utils.molecule_data import mol_data, mol_file, xyz_file
from comp_chem_utils.mysql_tables import mol_info_table, mol_xyz_table

#======================================================================================
#class atom_type_frame(object):
#    """ scrollable Tkinter object containing atom type sections"""
#
#    def __init__(self, root):
#        self.outer_frame = tk.Frame(root)
#        self.canvas = tk.Canvas(self.outer_frame)
#        self.inner_frame = tk.Frame(self.canvas)
#        self.scroll = tk.Scrollbar(self.outer_frame, orient="vertical",command=self.canvas.yview)
#        self.ypos = 100
#
#    def setup(self, xpos, ypos):
#        self.outer_frame.place(x=xpos+90, y=ypos)
#
#        self.canvas.configure(yscrollcommand=self.scroll.set)
#        self.canvas.pack(side="left") 
#        self.canvas.create_window((0,0),window=self.inner_frame,anchor='nw')
#
#        self.scroll.pack(side="right",fill="y")
#        self.inner_frame.bind("<Configure>", self.bind)
#
#    def bind(self, event):
#        self.canvas.configure(scrollregion=self.canvas.bbox("all")) #,width=200,height=200)
#
#    def add_section(self, atype):
#        """ add new atom type section to window """
#
#        # get new atom type information from mymol
#        atype_sec = atom_type_section_tk()
#        atype_sec.charge = atype.charge
#        atype_sec.natoms = atype.natoms
#        atype_sec.symb = atype.symb 
#        xvals = atype.xvals
#        yvals = atype.yvals
#        zvals = atype.zvals
#
#        dwi = 12
#        xinc = 2*dwi + 90
#        xpos = [0,100,200,300]
#        # first display header
#        atype_sec.lab_type = tk.Label(self.inner_frame, text = "Type:    "+atype_sec.symb)
#        atype_sec.lab_type.place()
#        #atype_sec.lab_type.place(x=xpos[0],y=self.ypos)
#        atype_sec.lab_x = tk.Label(self.inner_frame, text = "X")
#        #atype_sec.lab_x.place(x=self.xpos+xinc,y=self.ypos)
#        atype_sec.lab_x.place(x=xpos[1],y=self.ypos)
#        atype_sec.lab_y = tk.Label(self.inner_frame, text = "Y")
#        #atype_sec.lab_y.place(x=self.xpos+2*xinc,y=self.ypos)
#        atype_sec.lab_y.place(x=xpos[2],y=self.ypos)
#        atype_sec.lab_z = tk.Label(self.inner_frame, text = "Z")
#        #atype_sec.lab_z.place(x=self.xpos+3*xinc,y=self.ypos)
#        atype_sec.lab_z.place(x=xpos[3],y=self.ypos)
#        #self.ypos += self.inc
#
#        ## setup variables based on atype
#        #for x, y, z in zip(xvals, yvals, zvals):
#        #    atype_sec.xvals.append(tk.DoubleVar(value = x))
#        #    atype_sec.yvals.append(tk.DoubleVar(value = y))
#        #    atype_sec.zvals.append(tk.DoubleVar(value = z))
#          
#        ## display entries
#        #for x, y, z in zip(atype_sec.xvals, atype_sec.yvals, atype_sec.zvals):
#        #    atype_sec.entry_x.append( tk.Entry(self.inner_frame, textvariable = x, width=dwi) )
#        #    atype_sec.entry_y.append( tk.Entry(self.inner_frame, textvariable = y, width=dwi) )
#        #    atype_sec.entry_z.append( tk.Entry(self.inner_frame, textvariable = z, width=dwi) )
#          
#        #for x, y, z in zip(atype_sec.entry_x, atype_sec.entry_y, atype_sec.entry_z):
#        #    x.pack()
#        #    y.pack()
#        #    z.pack()
#        #    #x.place(x=self.xpos+xinc,y=self.ypos)
#        #    #y.place(x=self.xpos+2*xinc,y=self.ypos)
#        #    #z.place(x=self.xpos+3*xinc,y=self.ypos)
#        #    #self.ypos += self.inc
#
#        return atype_sec
#
#
#======================================================================================
class mol_tk(object):
    """ structure handling the main Tkinter window """

    def __init__(self):
        self.window = tk.Tk()

        self.idd = 0
        self.idd_lab = tk.StringVar()
        self.set_idd_lab()
        self.nat = 0
        self.nat_lab = tk.StringVar()
        self.set_nat_lab()

        self.name = tk.StringVar()
        self.symb = tk.StringVar()
        self.note = tk.Text()
        self.charge = tk.IntVar()
        self.atom_types_secs = [] # list of atom type Tk sections

        self.xpos = 30
        self.ypos = 15
        self.inc = 30

        self.setup_menubar()
        self.setup_window()
        #self.atom_type_frame = atom_type_frame(self.window)
        #self.atom_type_frame.setup( self.xpos, self.ypos)
        # display window
        self.window.mainloop()


    def setup_menubar(self):
        """ setup the top menu bar with the following options
        Open        (Database, .xyz file, .mol file (DALTON format))
        Save to...  (Database, .xyz file, .mol file (DALTON format)) 
        View with   (Avogadro, UCSF Chimera)
        Add atom type
        Clean up
        """
        menubar = tk.Menu(self.window)
     
        # create a pulldown menu, and add it to the menu bar
        entry_from = tk.Menu(menubar, tearoff=0)
        entry_from.add_command(label="Database", command=self.open_db_window)
        entry_from.add_command(label=".xyz file", command=self.open_xyz_file)
        entry_from.add_command(label=".mol file (DALTON format)", command=self.open_mol_file)
        menubar.add_cascade(label="Open", menu=entry_from)
        
        save_menu = tk.Menu(menubar, tearoff=0)
        save_menu.add_command(label="Database", command=self.save_to_db_window)
        save_menu.add_command(label=".xyz file", command=self.save_to_xyz)
        save_menu.add_command(label=".mol file (DALTON format)", command=self.save_to_mol)
        menubar.add_cascade(label="Save to...", menu=save_menu)
        
        view_menu = tk.Menu(menubar, tearoff=0)
        view_menu.add_command(label="Avogadro", command=self.open_avogadro)
        view_menu.add_command(label="UCSF Chimera", command=self.open_chimera)
        menubar.add_cascade(label="View with...", menu=view_menu)

        menubar.add_command(label="Add atom type", command=self.add_atom_type_window)
        menubar.add_command(label="Clean up", command=self.clean_up)
        
        # display the menu
        self.window.config(menu=menubar)


    def setup_window(self):
        """ setup main window with empty data """

        self.window.title("MSDB")
        self.window.geometry("550x800")
     
        xpos2 = self.xpos + 100
        # Create entry for the charge
        tk.Label(self.window, text = "Charge:").place(x=self.xpos,y=self.ypos)
        tk.Entry(self.window, textvariable = self.charge, width=5).place(x=xpos2,y=self.ypos)
     
        # Create entry for the database id
        off = xpos2 + 50
        tk.Label(self.window, textvariable = self.idd_lab).place(x=self.xpos+off,y=self.ypos)
     
        # Create entry for the number of atoms
        off = xpos2 + 180
        tk.Label(self.window, textvariable = self.nat_lab).place(x=self.xpos+off,y=self.ypos)
        self.ypos += self.inc
     
        # Create entry for the entry name
        tk.Label(self.window, text = "Entry name:").place(x=self.xpos,y=self.ypos)
        tk.Entry(self.window, textvariable = self.name, width=44).place(x=xpos2,y=self.ypos)
        self.ypos += self.inc
     
        # Create entry for note
        note = """Use this field to specify extra information about the molecular structure: 
     
        - source (website, article...), 
        - experimental technique,
        - computational model (specify the program used the basis set...)"""
        tk.Label(self.window, text = "Note:").place(x=self.xpos,y=self.ypos)
        self.note = tk.Text(self.window, width=50, height=15)
        self.note.insert(tk.INSERT,note)
        self.note.place(x=xpos2,y=self.ypos)
        self.ypos += 260
     

    # functions called under the OPEN menu in the main window
    # -------------------------------------------------------
    def open_db_window(self):
        """ display database in new top level window and let the user do his stuff"""
        top = tk.Toplevel()
        top.title("DATABASE")
        top.geometry("550x550")

        print "TODO: ask user for username and password"""
        connect = db.connect( host=hostname, user=username, passwd=password, db=database)
        cursor = connect.cursor()

        def return_db_id():
            """ read ID chosen by user and return it to main program"""
            # collect information from GUI
            cur_items = tree.selection()
            if len(cur_items)>1:
                self.error_window("Add only one entry at a time!")
                return
            elif len(cur_items)<=0:
                self.error_window("Please select one entry from the database.")
                return

            idd = tree.item(cur_items[0])['values'][0]
         
            # quit GUI
            top.destroy()
            top.quit()

            # get mol_data from database
            mymol = mol_data()
            err = mymol.init_from_db(cursor, idd)
            if err:
                self.error_window(err)
            else:
                # clean up and update main window:
                self.clean_up()
                self.set_values(mymol)

         
        def delete_db_entry_tk():
            # collect information from GUI
            cur_items = tree.selection()
            if len(cur_items)<=0:
                self.error_window("Please select at least one entry from the database.")
                return
            for item in cur_items:
                idd = tree.item(item)['values'][0]
         
                # delete record corresponding to idd in database
                delete_db_entries( cursor, [idd] )
                # delete record corresponding to idd in GUI table
                tree.delete(item)


        header = ['ID', 'Formula', 'Name']
        ask_user = ttk.Frame(top)
        ask_user.pack(fill='both', expand=True)
         
        # create a treeview with dual scrollbars
        tree = ttk.Treeview(top, columns=header, show="headings")
        vsb = ttk.Scrollbar(top, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(top, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(column=0, row=0, sticky='nsew', in_=ask_user)
        vsb.grid(column=1, row=0, sticky='ns', in_=ask_user)
        hsb.grid(column=0, row=1, sticky='ew', in_=ask_user)
        ask_user.grid_columnconfigure(0, weight=1)
        ask_user.grid_rowconfigure(0, weight=1)

        # fill table
        for col in header:
            tree.heading(col, text=col.title())
            # adjust the column's width to the header string
            tree.column(col, width=tkFont.Font().measure(col.title()))

        # list of all available molecules:
        cursor.execute( mol_info_table().get_col( ['id','chem_name','name'] ) )
         
        for item in cursor:
            tree.insert('', 'end', values=item)
            # adjust column's width if necessary to fit each value
            for ix, val in enumerate(item):
                col_w = tkFont.Font().measure(val)
                if tree.column(header[ix],width=None)<col_w:
                    tree.column(header[ix], width=col_w)

        # Open/delete buttons
        choice = tk.Frame(top)
        choice.pack(side="top")
        open_b = tk.Button(choice, text = "Open", command = return_db_id)
        open_b.pack(side="left", ipadx=3, ipady=3, padx=20,pady=10)
        delete_b = tk.Button(choice, text = "Delete", command = delete_db_entry_tk)
        delete_b.pack(side="left", ipadx=3, ipady=3, padx=20,pady=10)

        # display
        top.mainloop()

        # commit changes to database and close connection
        connect.commit()
        connect.close()


    def open_mol_file(self):
        """ open mol file and display in window """
        mol = mol_file()
        filename = tkfd.askopenfilename(title = "Select .mol file",
            filetypes = (("mol files","*.mol"),("all files","*.*")))
        # read mol file
        err1 = mol.read_mol(filename)
        if err1:
            self.error_window(err1)
        else:
            # init mol_dat from mol_file
            mymol=mol_data()
            err2 = mymol.init_from_mol(mol)
            if err2:
                self.error_window(err2)
            else:
                # clean up and update window:
                self.clean_up()
                self.set_values(mymol)

    def open_xyz_file(self):
        """ open xyz file and display in window """
        xyz = xyz_file()
        filename = tkfd.askopenfilename(title = "Select .xyz file",
            filetypes = (("xyz files","*.xyz"),("all files","*.*")))
        # read xyz file
        err1 = xyz.read_xyz(filename)
        if err1:
            self.error_window(err1)
        else:
            # init mol_dat from xyz_file
            mymol=mol_data()
            err2 = mymol.init_from_xyz(xyz)
            if err2:
                self.error_window(err2)
            else: 
                # clean up and update window:
                self.clean_up()
                self.set_values(mymol)


    # functions called under the SAVE menu in the main window
    # -------------------------------------------------------
    def save_to_db_window(self):
        """ ask confirmation from user before saving to database """

        # create window
        top = tk.Toplevel()
        top.title("Saving to Database")
        name = self.name.get()
        if not name:
            name = "unknown"

        # YES function from top level window
        def save_to_db():
            """ Save mol_data to database"""
            # destroy user confirmation window
            top.destroy()
            top.quit()
             
            # open database
            print "TODO: ask user for username and password"""
            connect = db.connect( host=hostname, user=username, passwd=password, db=database)
            cursor = connect.cursor()
            # get mol_data from window:
            mymol = self.get_values()
            # check for duplicates
            ddl = mol_info_table().find_duplicates( mymol.chem_name )
            cursor.execute(ddl)
            idds = [idd for (idd, name) in cursor]

            what_to_do = 0
            if idds:
                what_to_do = self.duplicate_window( cursor )

            if what_to_do==0:
                # default: save entry to db
                err = mymol.out_to_db(cursor)
                if err:
                    self.error_window(err)
            elif what_to_do==1:
                # overwite duplicated entry 
                # VERY DANGEROUS!! I guess...
                # first we delete all entries in idds
                delete_db_entries( cursor, idds )

                # then we save the current entry
                err = mymol.out_to_db(cursor)
                if err:
                    self.error_window(err)
            else:
                # USER canceled...
                pass

            # save eventual changes to the database and close the connection
            connect.commit()
            connect.close()

        # NO functions from top level window
        def cancel_save_to_db():
            top.destroy()
            top.quit()

        # main question
        text = "Do you want to save the current entry ({0}) to the database? ".format(name)
        tk.Label(top, text = text).pack(padx=20,pady=20)

        # yes/no buttons
        yesno = tk.Frame(top)
        yesno.pack(side="top")
        yes = tk.Button(yesno, text = "Yes", command = save_to_db)
        yes.pack(side="left", ipadx=3, ipady=3, padx=10,pady=(0, 20))
        no = tk.Button(yesno, text = "No", command = cancel_save_to_db)
        no.pack(side="left", ipadx=3, ipady=3, padx=10,pady=(0, 20))

        # display
        top.mainloop()

    def save_to_xyz(self, temp=False):
        """ save GUI info to xyz file """
        # get mol_data from window:
        mymol = self.get_values()

        # save to disk
        xyz = mymol.out_to_xyz()
        if temp:
            xyz.fname = ".temp.xyz"
        else:
            xyz.fname = tkfd.asksaveasfilename(title = "Select .xyz file",
                filetypes = (("xyz files","*.xyz"),("all files","*.*")),
                initialfile = xyz.fname)

        if xyz.fname:
            err = xyz.out_to_file()
            if err:
                self.error_window(err)

        return xyz.fname 

    def save_to_mol(self):
        """ save GUI info to mol file """
        # get mol_data from window:
        mymol = self.get_values()

        # save to disk
        mol = mymol.out_to_mol()
        mol.fname = tkfd.asksaveasfilename(title = "Select .mol file",
                filetypes = (("mol files","*.mol"),("all files","*.*")),
                initialfile = mol.fname)
        if mol.fname:
            err = mol.out_to_file()
            if err:
                self.error_window(err)


    # functions called under the VIEW menu in the main window
    # -------------------------------------------------------
    def open_avogadro(self):
        """ open current molecule with avogadro"""
        xyzfname = self.save_to_xyz(True)
        os.system("avogadro "+xyzfname)
        os.remove(xyzfname)

    def open_chimera(self):
        """ open current molecule with USCF chimera"""
        xyzfname = self.save_to_xyz(True)
        os.system("chimera "+xyzfname)
        os.remove(xyzfname)


    # remaining functions called from the menu in the main window
    # -----------------------------------------------------------
    def add_atom_type_window(self):
        """ display top window with atoms to be chosen by user """
        top = tk.Toplevel()
        top.title("New atom type")
        top.geometry("340x400")

        def return_info():
            """ read info from top window and update main window"""
            # create empty atom type from GUI
            atype = atom_type()
            atype.natoms = natoms_in_type.get()
            if atype.natoms<=0:
                self.error_window("Enter positive number of atoms!")
                return

            cur_item = tree.selection()
            if len(cur_item) > 1:
                self.error_window("Select only one atom type at a time!")
                return

            atype.symb = tree.item(cur_item[0])['values'][0]
            atype.charge = tree.item(cur_item[0])['values'][1]

            # put zeros in coordinates
            atype.xvals = [0.0] * atype.natoms
            atype.yvals = [0.0] * atype.natoms
            atype.zvals = [0.0] * atype.natoms
             
            # add new atom type to main window:
            self.nat += atype.natoms
            self.set_nat_lab()
            self.add_atom_type(atype)

        header = ['Symbol', 'Charge']
        ask_user = ttk.Frame(top)
        ask_user.pack(fill='both', expand=True)
         
        # create a treeview with dual scrollbars
        tree = ttk.Treeview(top, columns=header, show="headings")
        vsb = ttk.Scrollbar(top, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(top, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(column=0, row=0, sticky='nsew', in_=ask_user)
        vsb.grid(column=1, row=0, sticky='ns', in_=ask_user)
        hsb.grid(column=0, row=1, sticky='ew', in_=ask_user)
        ask_user.grid_columnconfigure(0, weight=1)
        ask_user.grid_rowconfigure(0, weight=1)

        # fill table
        for col in header:
            tree.heading(col, text=col.title())
            # adjust the column's width to the header string
            tree.column(col, width=tkFont.Font().measure(col.title()))

        for i , symb in enumerate(atoms):
            item = (symb, i + 1)
            tree.insert('', 'end', values=item)
            # adjust column's width if necessary to fit each value
            for ix, val in enumerate(item):
                col_w = tkFont.Font().measure(val)
                if tree.column(header[ix],width=None)<col_w:
                    tree.column(header[ix], width=col_w)

        # add number of atoms section
        natoms_in_type = tk.IntVar(value = 1)
        lab = tk.Label(top, text = "Number of atoms:")
        lab.pack(side="left", ipadx=3, ipady=3, padx=10, pady=10 )
        ent = tk.Entry(top, textvariable = natoms_in_type, width=5)
        ent.pack(side="left", ipadx=3, ipady=3, padx=10, pady=10 )
        but = tk.Button(top, text = "Add", command = return_info)
        but.pack(side="left", ipadx=3, ipady=3, padx=40,pady=10 )

        top.mainloop()


    def clean_up(self):
        """ update GUI with empty mol_data"""
        # first we remove atom type sections
        for atype_sec in self.atom_types_secs:
            atype_sec.delete(self)
        self.atom_types_secs = [] # re-init
        mymol=mol_data()
        self.set_values(mymol)


    # functions to set/get the information from the main window to the mol_data structure:
    # ------------------------------------------------------------------------------------
    def set_values(self, mymol):
        """ set values in GUI based on mymol structure """
        self.idd = mymol.idd
        self.set_idd_lab()
        self.name.set(mymol.name)
        self.note.delete(1.0, tk.END)
        self.note.insert(tk.INSERT, mymol.note)
        self.charge.set(mymol.charge)

        # create atom type sections:
        nat = 0
        for atype in mymol.atom_types:
            self.add_atom_type(atype)
            nat += atype.natoms

        if nat!=mymol.natoms:
            print "ERROR: Inconsistency in mol_data (natoms)",nat,mymol.natoms
            self.error_window("Inconsistency in mol_data (natoms)")
        else:
            # set natoms
            self.nat = mymol.natoms
            self.set_nat_lab()

    def get_values(self):
        """ read values from GUI and output new mymol"""
        mymol = mol_data()
        mymol.name = self.name.get()
        mymol.note = self.note.get(1.0, tk.END).strip()
        mymol.charge = self.charge.get()
        mymol.natoms = self.nat

        # create atom type sections:
        for atype_sec in self.atom_types_secs:
            atype = atom_type()
            atype.set_from_GUI(atype_sec)
            mymol.atom_types.append(atype)

        mymol.natom_types = len(mymol.atom_types)
        mymol.get_chem_name()
        return mymol


    # other tool functions
    # --------------------
    def add_atom_type(self, atype):
        """ add new atom type section to window """

        #atype_sec = self.atom_type_frame.add_section( atype )

        # get new atom type information from mymol
        atype_sec = atom_type_section_tk()
        atype_sec.charge = atype.charge
        atype_sec.natoms = atype.natoms
        atype_sec.symb = atype.symb 
        xvals = atype.xvals
        yvals = atype.yvals
        zvals = atype.zvals
         
        dwi = 12
        xinc = 2*dwi + 90
        # first display header
        atype_sec.lab_type = tk.Label(self.window, text = "Type:    "+atype_sec.symb)
        atype_sec.lab_type.place(x=self.xpos,y=self.ypos)
        atype_sec.lab_x = tk.Label(self.window, text = "X")
        atype_sec.lab_x.place(x=self.xpos+xinc,y=self.ypos)
        atype_sec.lab_y = tk.Label(self.window, text = "Y")
        atype_sec.lab_y.place(x=self.xpos+2*xinc,y=self.ypos)
        atype_sec.lab_z = tk.Label(self.window, text = "Z")
        atype_sec.lab_z.place(x=self.xpos+3*xinc,y=self.ypos)
        self.ypos += self.inc
         
        # setup variables based on atype
        for x, y, z in zip(xvals, yvals, zvals):
            atype_sec.xvals.append(tk.DoubleVar(value = x))
            atype_sec.yvals.append(tk.DoubleVar(value = y))
            atype_sec.zvals.append(tk.DoubleVar(value = z))
         
        # display entries
        for x, y, z in zip(atype_sec.xvals, atype_sec.yvals, atype_sec.zvals):
            atype_sec.entry_x.append( tk.Entry(self.window, textvariable = x, width=dwi) )
            atype_sec.entry_y.append( tk.Entry(self.window, textvariable = y, width=dwi) )
            atype_sec.entry_z.append( tk.Entry(self.window, textvariable = z, width=dwi) )
         
        for x, y, z in zip(atype_sec.entry_x, atype_sec.entry_y, atype_sec.entry_z):
            x.place(x=self.xpos+xinc,y=self.ypos)
            y.place(x=self.xpos+2*xinc,y=self.ypos)
            z.place(x=self.xpos+3*xinc,y=self.ypos)
            self.ypos += self.inc

        # add new atom type section to self
        self.atom_types_secs.append(atype_sec)

    def error_window(self, err):
        """ display error message on top level window"""
        print "ERROR:",err
        top = tk.Toplevel()
        top.title("ERROR!")
        tk.Label(top, text = err).pack(padx=20,pady=20)
        top.mainloop()

    def duplicate_window(self, cursor):
        """ tell user we found duplicated entries and ask what to do about it"""
        top = tk.Toplevel()
        top.title("Duplicated entries in database")

        what_to_do = [-1] # need to declare as list for scope reasons
        # functions to return answer from user
        def replace():
            # quit top level window and return 1
            what_to_do[0] = 1
            top.destroy()
            top.quit()
        def saveas():
            # quit top level window and return 0
            what_to_do[0] = 0
            top.destroy()
            top.quit()
        def cancel():
            # quit top level window and return -1
            what_to_do[0] = -1
            top.destroy()
            top.quit()

        # show label
        label = "one or several entries with the same chemical \nformula already exist in the database, ID(s): "
        tk.Label(top, text = label).pack(padx=20,pady=20)
        for (idd, name) in cursor:
            label = "{:4d} {}".format(idd, name)
            tk.Label(top, text = label).pack(padx=20,pady=5)

        # add answer button
        choice = tk.Frame(top)
        choice.pack(side="top")

        repl = tk.Button(choice, text = "Replace all!", command = replace)
        repl.pack(side="left", ipadx=3, ipady=3, padx=20,pady=10)
        save = tk.Button(choice, text = "Save as new", command = saveas)
        save.pack(side="left", ipadx=3, ipady=3, padx=20,pady=10)
        canc = tk.Button(choice, text = "Cancel", command = cancel)
        canc.pack(side="left", ipadx=3, ipady=3, padx=20,pady=10)

        top.mainloop()
        return what_to_do[0]

    def set_nat_lab(self):
        """ set GUI value for natoms """
        self.nat_lab.set("Number of atoms: {0.nat:6d}".format(self))

    def set_idd_lab(self):
        """ set GUI value for database ID """
        self.idd_lab.set("Database ID: {0.idd:6d}".format(self))

#======================================================================================
class atom_type_section_tk(object):
    """ structure handling atom type sections in Tkinter window """

    def __init__(self):
        self.charge = 0
        self.natoms = 0
        self.symb = ''

        # header labels
        self.lab_type = None
        self.lab_x = None
        self.lab_y = None
        self.lab_z = None

        # coordinates info
        self.xvals = []
        self.yvals = []
        self.zvals = []
        self.entry_x = []
        self.entry_y = []
        self.entry_z = []

    def delete(self, mymol_tk):
        """ delete atom type section from GUI"""
        self.lab_type.destroy()
        self.lab_x.destroy()
        self.lab_y.destroy()
        self.lab_z.destroy()
        mymol_tk.ypos -= mymol_tk.inc # update ypos

        for x, y, z in  zip(self.entry_x, self.entry_y, self.entry_z):
            x.destroy()
            y.destroy()
            z.destroy()
            mymol_tk.ypos -= mymol_tk.inc # update ypos


#======================================================================================
def delete_db_entries(cursor, idds):
    """ delete set of entries from the data base """

    def delete():
        for idd in idds:
            cursor.execute( mol_info_table().delete_row( idd ) )
   
            # delete coordinate table
            cursor.execute( mol_xyz_table( idd ).delete_table() ) 
            print "Entry with ID = ", idd, " deleted from database"

        # quit GUI
        top.destroy()
        top.quit()
        return

    def cancel():
        # quit GUI
        top.destroy()
        top.quit()
        return

    # first we ask confirmation from user
    top = tk.Toplevel()
    top.title("Really?!")
    label = "Doyou really want to delete the following entries from the database?"
    tk.Label(top, text = label).pack(padx=20,pady=20)
    for idd in idds:
        tk.Label(top, text = str(idd) ).pack(padx=20,pady=5)

    # yes/no button
    choice = tk.Frame(top)
    choice.pack(side="top")
    yes_b = tk.Button(choice, text = "Yes", command = delete)
    yes_b.pack(side="left", ipadx=3, ipady=3, padx=20,pady=10)
    no_b = tk.Button(choice, text = "No", command = cancel)
    no_b.pack(side="left", ipadx=3, ipady=3, padx=20,pady=10)

    top.mainloop()


if __name__ == '__main__':
    # set up and display GUI
    mol_tk()

