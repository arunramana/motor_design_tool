
# import tkinter module
from tkinter import *
from tkinter.ttk import *
from PIL import Image, ImageTk

from dxf import *
from run_femm import *

import matplotlib.pyplot as plt
import ezdxf
from ezdxf.addons.drawing import RenderContext, Frontend
from ezdxf.addons.drawing.matplotlib import MatplotlibBackend
# import wx
import glob
import re

from functools import partial


#simple demo using tkinter

class DXF2IMG(object):

    default_img_format = '.png'
    default_img_res = 300
    def convert_dxf2img(self, names, img_format=default_img_format, img_res=default_img_res):
        for name in names:
            doc = ezdxf.readfile(name)
            msp = doc.modelspace()
            # Recommended: audit & repair DXF document before rendering
            auditor = doc.audit()
            # The auditor.errors attribute stores severe errors,
            # which *may* raise exceptions when rendering.
            if len(auditor.errors) != 0:
                raise exception("The DXF document is damaged and can't be converted!")
            else :
                fig = plt.figure()
                ax = fig.add_axes([0, 0, 1, 1])
                ctx = RenderContext(doc)
                ctx.set_current_layout(msp)
                ctx.current_layout.set_colors(bg='#FFFFFF')
                out = MatplotlibBackend(ax)
                Frontend(ctx, out).draw_layout(msp, finalize=True)

                img_name = re.findall("(\S+)\.",name)  # select the image name that is the same as the dxf file name
                first_param = ''.join(img_name) + img_format  #concatenate list and string
                fig.savefig(first_param, dpi=img_res)


def generate_dxf(master,labelname,no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distance1 , duct_radius1 , duct_angle1 ,duct_distance2 , duct_radius2 , duct_angle2, duct_distance3 , duct_radius3 , duct_angle3, notch_depth , notch_angle , no_turns):

    no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distance1 , duct_radius1 , duct_angle1 ,duct_distance2 , duct_radius2 , duct_angle2, duct_distance3 , duct_radius3 , duct_angle3, notch_depth , notch_angle , no_turns = no_slots.get(), no_poles.get() , stator_od.get() , stator_id.get() , rotor_od.get() , shaft_dia.get() , slot_opening.get() , tooth_tip_angle.get() , tooth_width.get() , yoke_thickness.get() , pole_length.get() , pole_width.get() , pole_arc.get() , bridge_thickness.get() , duct_distance1.get() , duct_radius1.get() , duct_angle1.get() ,duct_distance2.get() , duct_radius2.get() , duct_angle2.get(), duct_distance3.get() , duct_radius3.get() , duct_angle3.get(), notch_depth.get() , notch_angle.get() , no_turns.get()

    fname = "test4.dxf"

    width= master.winfo_screenwidth()*0.99
    height= master.winfo_screenheight()*0.99


    # print("dims",no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distance1 , duct_radius1 , duct_angle1 ,duct_distance2 , duct_radius2 , duct_angle2, duct_distance3 , duct_radius3 , duct_angle3, notch_depth , notch_angle , no_turns)

    no_slots = 12 if no_slots=="" or no_slots==" " else int(no_slots)

    no_poles = 10 if no_poles=="" or no_poles==" " else int(no_poles)


    stator_od = 125 if stator_od=="" or stator_od==" " else int(stator_od)
    stator_id = 84 if stator_id=="" or stator_id==" " else int(stator_id)
    rotor_od = 83 if rotor_od=="" or rotor_od==" " else int(rotor_od)
    shaft_dia = 20 if shaft_dia=="" or shaft_dia==" " else int(shaft_dia)

    slot_opening = 3 if slot_opening=="" or slot_opening==" " else int(slot_opening)

    tooth_tip_angle = 20 if tooth_tip_angle=="" or tooth_tip_angle==" " else int(tooth_tip_angle)

    tooth_width = 10 if tooth_width=="" or tooth_width==" " else int(tooth_width)

    yoke_thickness = 5 if yoke_thickness=="" or yoke_thickness==" " else int(yoke_thickness)


    pole_length = 18 if pole_length=="" or pole_length==" " else int(pole_length)
    pole_width = 5 if pole_width=="" or pole_width==" " else int(pole_width)
    pole_arc = 150 if pole_arc=="" or pole_arc==" " else int(pole_arc)

    bridge_thickness = 1 if bridge_thickness=="" or bridge_thickness==" " else float(bridge_thickness)


    #the three arrays below should be of same dimensions

    duct_distances = [duct_distance1,duct_distance2,duct_distance3]

    duct_angles = [duct_angle1,duct_angle2,duct_angle3]

    duct_radii = [duct_radius1,duct_radius2,duct_radius3]

    if "" in duct_distances or " " in duct_distances:
        duct_distances = [30,25,30]

    else:
        duct_distances = [int(duct_distance1),int(duct_distance2),int(duct_distance3)]



    if "" in duct_radii or " " in duct_radii:
        duct_radii = [2,2,2]
    else:
        duct_radii = [int(duct_radius1),int(duct_radius2),int(duct_radius3)]



    if "" in duct_angles or " " in duct_angles:
        duct_angles = [-18,0,18]

    else:
        duct_angles = duct_angles = [int(duct_angle1),(duct_angle2),(duct_angle3)]

    notch_depth = 4 if notch_depth=="" or notch_depth==" " else int(notch_depth)
    notch_angle = 7 if notch_angle=="" or notch_angle==" " else int(notch_angle)


    no_turns = 9 if no_turns=="" or no_turns==" " else int(no_turns)

    gdxf = GenerateDXF(fname, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns)

    gdxf.draw_motor()

    first =  DXF2IMG()
    b = first.convert_dxf2img(['test4.dxf'],img_format='.jpg')

    img_old = Image.open("test4.jpg")
    img_resized=img_old.resize((int(width*0.55),int(height*0.85)))
    photo = ImageTk.PhotoImage(img_resized)

    labelname.configure(image=photo)

    labelname.image = photo

def generate_femm(no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distance1 , duct_radius1 , duct_angle1 ,duct_distance2 , duct_radius2 , duct_angle2, duct_distance3 , duct_radius3 , duct_angle3, notch_depth , notch_angle , no_turns):

    no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distance1 , duct_radius1 , duct_angle1 ,duct_distance2 , duct_radius2 , duct_angle2, duct_distance3 , duct_radius3 , duct_angle3, notch_depth , notch_angle , no_turns = no_slots.get(), no_poles.get() , stator_od.get() , stator_id.get() , rotor_od.get() , shaft_dia.get() , slot_opening.get() , tooth_tip_angle.get() , tooth_width.get() , yoke_thickness.get() , pole_length.get() , pole_width.get() , pole_arc.get() , bridge_thickness.get() , duct_distance1.get() , duct_radius1.get() , duct_angle1.get() ,duct_distance2.get() , duct_radius2.get() , duct_angle2.get(), duct_distance3.get() , duct_radius3.get() , duct_angle3.get(), notch_depth.get() , notch_angle.get() , no_turns.get()

    fname = "test4.dxf"

    width= master.winfo_screenwidth()*0.99
    height= master.winfo_screenheight()*0.99


    # print("dims",no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distance1 , duct_radius1 , duct_angle1 ,duct_distance2 , duct_radius2 , duct_angle2, duct_distance3 , duct_radius3 , duct_angle3, notch_depth , notch_angle , no_turns)

    no_slots = 12 if no_slots=="" or no_slots==" " else int(no_slots)

    no_poles = 10 if no_poles=="" or no_poles==" " else int(no_poles)


    stator_od = 125 if stator_od=="" or stator_od==" " else int(stator_od)
    stator_id = 84 if stator_id=="" or stator_id==" " else int(stator_id)
    rotor_od = 83 if rotor_od=="" or rotor_od==" " else int(rotor_od)
    shaft_dia = 20 if shaft_dia=="" or shaft_dia==" " else int(shaft_dia)

    slot_opening = 3 if slot_opening=="" or slot_opening==" " else int(slot_opening)

    tooth_tip_angle = 20 if tooth_tip_angle=="" or tooth_tip_angle==" " else int(tooth_tip_angle)

    tooth_width = 10 if tooth_width=="" or tooth_width==" " else int(tooth_width)

    yoke_thickness = 5 if yoke_thickness=="" or yoke_thickness==" " else int(yoke_thickness)


    pole_length = 18 if pole_length=="" or pole_length==" " else int(pole_length)
    pole_width = 5 if pole_width=="" or pole_width==" " else int(pole_width)
    pole_arc = 150 if pole_arc=="" or pole_arc==" " else int(pole_arc)

    bridge_thickness = 1 if bridge_thickness=="" or bridge_thickness==" " else float(bridge_thickness)


    #the three arrays below should be of same dimensions

    duct_distances = [duct_distance1,duct_distance2,duct_distance3]

    duct_angles = [duct_angle1,duct_angle2,duct_angle3]

    duct_radii = [duct_radius1,duct_radius2,duct_radius3]

    if "" in duct_distances or " " in duct_distances:
        duct_distances = [30,25,30]

    else:
        duct_distances = [int(duct_distance1),int(duct_distance2),int(duct_distance3)]



    if "" in duct_radii or " " in duct_radii:
        duct_radii = [2,2,2]
    else:
        duct_radii = [int(duct_radius1),int(duct_radius2),int(duct_radius3)]



    if "" in duct_angles or " " in duct_angles:
        duct_angles = [-18,0,18]

    else:
        duct_angles = duct_angles = [int(duct_angle1),(duct_angle2),(duct_angle3)]

    notch_depth = 4 if notch_depth=="" or notch_depth==" " else int(notch_depth)
    notch_angle = 7 if notch_angle=="" or notch_angle==" " else int(notch_angle)


    no_turns = 9 if no_turns=="" or no_turns==" " else int(no_turns)

    u,v,w = 0,25,-25

    stack_depth = 36

    rf = GenerateFEMM(fname, no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distances , duct_radii , duct_angles , notch_depth , notch_angle , no_turns, stack_depth,u, v, w)

    rf.run_femm()


# creating main tkinter window/toplevel

def create_window():
    master = Toplevel()

    width= master.winfo_screenwidth()*0.99
    height= master.winfo_screenheight()*0.99
    master.geometry("%dx%d" % (width, height))

    master.title("Motor Design")

    return master

# this will create a label widget

def GUI(master):

    width= master.winfo_screenwidth()*0.99
    height= master.winfo_screenheight()*0.99

    no_slots = StringVar(master)
    no_poles = StringVar(master)
    stator_od = StringVar(master)
    stator_id = StringVar(master)
    rotor_od = StringVar(master)
    shaft_dia = StringVar(master)
    slot_opening = StringVar(master)
    tooth_tip_angle = StringVar(master)
    tooth_width = StringVar(master)
    yoke_thickness = StringVar(master)
    pole_length = StringVar(master)
    pole_width = StringVar(master)
    pole_arc = StringVar(master)
    bridge_thickness = StringVar(master)

    duct_distance1 = StringVar()
    duct_radius1 = StringVar()
    duct_angle1 = StringVar()

    duct_distance2 = StringVar()
    duct_radius2 = StringVar()
    duct_angle2 = StringVar()

    duct_distance3 = StringVar()
    duct_radius3 = StringVar()
    duct_angle3 = StringVar()

    notch_depth = StringVar()
    notch_angle = StringVar()
    no_turns = StringVar()


    l1 = Label(master, text = "Stator")
    l2 = Label(master, text = "Stator OD (mm)")
    l3 = Label(master, text = "Stator ID (mm)")
    l4 = Label(master, text = "No. Slots (No)")
    l5 = Label(master, text = "Slot Opening (mm)")
    l6 = Label(master, text = "Tooth Tip angle (deg)")
    l7 = Label(master, text = "Tooth Width (mm)")
    l8 = Label(master, text = "Yoke Thickness (mm)")

    l9 = Label(master, text = "Rotor")
    l10 = Label(master, text = "Rotor OD (mm)")
    l11 = Label(master, text = "Shaft Dia (mm)")
    l12 = Label(master, text = "No. Poles (No)")
    l13 = Label(master, text = "Pole Length (mm)")
    l14 = Label(master, text = "Pole Width (mm)")
    l15 = Label(master, text = "Pole Arc (deg)")
    l16 = Label(master, text = "Bridge Thickness (mm)")

    l17 = Label(master, text = "Duct Distance 1 (mm)")
    l18 = Label(master, text = "Duct Radius 1 (mm)")
    l19 = Label(master, text = "Duct Angle 1 (Deg)")

    l20 = Label(master, text = "Duct Distance 2 (mm)")
    l21 = Label(master, text = "Duct Radius 2 (mm)")
    l22 = Label(master, text = "Duct Angle 2 (Deg)")

    l23 = Label(master, text = "Duct Distance 3 (mm)")
    l24 = Label(master, text = "Duct Radius 3 (mm)")
    l25 = Label(master, text = "Duct Angle 3 (Deg)")

    l26 = Label(master, text = "Notch Depth (mm)")
    l27 = Label(master, text = "Notch Angle (Deg)")

    l28 = Label(master, text = "Circuit")
    l29 = Label(master, text = "No. Turns (No)")
    l30 = Label(master, text = "U Current (A)")
    l31 = Label(master, text = "V Current (A)")
    l32 = Label(master, text = "W Current (A)")


    l33 = Label(master, text = "Materials")
    l34 = Label(master, text = "Copper")
    l35 = Label(master, text = "Magnet")
    l36 = Label(master, text = "Steel")


    # grid method to arrange labels in respective
    # rows and columns as specified

    #stator

    l1.grid(row = 0, column = 1, sticky = W, pady = 2)

    l2.grid(row = 1, column = 0, sticky = W, pady = 2)
    l3.grid(row = 2, column = 0, sticky = W, pady = 2)
    l4.grid(row = 3, column = 0, sticky = W, pady = 2)
    l5.grid(row = 4, column = 0, sticky = W, pady = 2)
    l6.grid(row = 5, column = 0, sticky = W, pady = 2)
    l7.grid(row = 6, column = 0, sticky = W, pady = 2)
    l8.grid(row = 7, column = 0, sticky = W, pady = 2)

    #rotor
    l9.grid(row = 8, column = 1, sticky = W, pady = 2)

    l10.grid(row = 9, column = 0, sticky = W, pady = 2)
    l11.grid(row = 10, column = 0, sticky = W, pady = 2)
    l12.grid(row = 11, column = 0, sticky = W, pady = 2)
    l13.grid(row = 12, column = 0, sticky = W, pady = 2)
    l14.grid(row = 13, column = 0, sticky = W, pady = 2)
    l15.grid(row = 14, column = 0, sticky = W, pady = 2)
    l16.grid(row = 15, column = 0, sticky = W, pady = 2)

    #ducts
    l17.grid(row = 17, column = 0, sticky = W, pady = 2)
    l18.grid(row = 18, column = 0, sticky = W, pady = 2)
    l19.grid(row = 19, column = 0, sticky = W, pady = 2)

    l20.grid(row = 21, column = 0, sticky = W, pady = 2)
    l21.grid(row = 22, column = 0, sticky = W, pady = 2)
    l22.grid(row = 23, column = 0, sticky = W, pady = 2)

    l23.grid(row = 25, column = 0, sticky = W, pady = 2)
    l24.grid(row = 26, column = 0, sticky = W, pady = 2)
    l25.grid(row = 27, column = 0, sticky = W, pady = 2)

    #notch
    l26.grid(row = 29, column = 0, sticky = W, pady = 2)
    l27.grid(row = 30, column = 0, sticky = W, pady = 2)

    #Circuits
    l28.grid(row = 31, column = 1, sticky = W, pady = 2)
    l29.grid(row = 32, column = 0, sticky = W, pady = 2)
    l30.grid(row = 33, column = 0, sticky = W, pady = 2)
    l31.grid(row = 34, column = 0, sticky = W, pady = 2)
    l32.grid(row = 35, column = 0, sticky = W, pady = 2)

    '''
    #Materials
    l33.grid(row = 36, column = 1, sticky = W, pady = 2)
    l34.grid(row = 37, column = 0, sticky = W, pady = 2)
    l35.grid(row = 38, column = 0, sticky = W, pady = 2)
    l36.grid(row = 39, column = 0, sticky = W, pady = 2)
    '''

    # entry widgets, used to take entry from user
    #e1 = Entry(master)
    e2 = Entry(master,textvariable=stator_od)
    e3 = Entry(master,textvariable=stator_id)
    e4 = Entry(master,textvariable=no_slots)
    e5 = Entry(master,textvariable=slot_opening)
    e6 = Entry(master,textvariable=tooth_tip_angle)
    e7 = Entry(master,textvariable=tooth_width)
    e8 = Entry(master,textvariable=yoke_thickness)
    #e9 = Entry(master)
    e10 = Entry(master,textvariable=rotor_od)
    e11 = Entry(master,textvariable=shaft_dia)
    e12 = Entry(master,textvariable=no_poles)
    e13 = Entry(master,textvariable=pole_length)
    e14 = Entry(master,textvariable=pole_width)
    e15 = Entry(master,textvariable=pole_arc)
    e16 = Entry(master,textvariable=bridge_thickness)
    e17 = Entry(master,textvariable=duct_distance1)
    e18 = Entry(master,textvariable=duct_radius1)
    e19 = Entry(master,textvariable=duct_angle1)
    e20 = Entry(master,textvariable=duct_distance2)
    e21 = Entry(master,textvariable=duct_radius2)
    e22 = Entry(master,textvariable=duct_angle2)
    e23 = Entry(master,textvariable=duct_distance3)
    e24 = Entry(master,textvariable=duct_radius3)
    e25 = Entry(master,textvariable=duct_angle3)
    e26 = Entry(master,textvariable=notch_depth)
    e27 = Entry(master,textvariable=notch_angle)
    #e28 = Entry(master)
    e29 = Entry(master,textvariable=no_turns)
    e30 = Entry(master)
    e31 = Entry(master)
    e32 = Entry(master)

    '''
    #e33 = Entry(master)
    e34 = Entry(master)
    e35 = Entry(master)
    e36 = Entry(master)
    '''

    # this will arrange entry widgets
    #e1.grid(row = 0, column = 1, pady = 2)
    e2.grid(row = 1, column = 1, pady = 2)
    e3.grid(row = 2, column = 1, pady = 2)
    e4.grid(row = 3, column = 1, pady = 2)
    e5.grid(row = 4, column = 1, pady = 2)
    e6.grid(row = 5, column = 1, pady = 2)
    e7.grid(row = 6, column = 1, pady = 2)
    e8.grid(row = 7, column = 1, pady = 2)
    #e9.grid(row = 8, column = 1, pady = 2)
    e10.grid(row = 9, column = 1, pady = 2)
    e11.grid(row = 10, column = 1, pady = 2)
    e12.grid(row = 11, column = 1, pady = 2)
    e13.grid(row = 12, column = 1, pady = 2)
    e14.grid(row = 13, column = 1, pady = 2)
    e15.grid(row = 14, column = 1, pady = 2)
    e16.grid(row = 15, column = 1, pady = 2)
    e17.grid(row = 17, column = 1, pady = 2)
    e18.grid(row = 18, column = 1, pady = 2)
    e19.grid(row = 19, column = 1, pady = 2)
    e20.grid(row = 21, column = 1, pady = 2)
    e21.grid(row = 22, column = 1, pady = 2)
    e22.grid(row = 23, column = 1, pady = 2)
    e23.grid(row = 25, column = 1, pady = 2)
    e24.grid(row = 26, column = 1, pady = 2)
    e25.grid(row = 27, column = 1, pady = 2)
    e26.grid(row = 29, column = 1, pady = 2)
    e27.grid(row = 30, column = 1, pady = 2)
    #e28.grid(row = 27, column = 1, pady = 2)
    e29.grid(row = 32, column = 1, pady = 2)
    e30.grid(row = 33, column = 1, pady = 2)
    e31.grid(row = 34, column = 1, pady = 2)
    e32.grid(row = 35, column = 1, pady = 2)

    '''
    #e33.grid(row = 32, column = 1, pady = 2)
    e34.grid(row = 37, column = 1, pady = 2)
    e35.grid(row = 38, column = 1, pady = 2)
    e36.grid(row = 39, column = 1, pady = 2)
    '''

    img_old = Image.open("test4.jpg")
    img_resized=img_old.resize((int(width*0.55),int(height*0.85)))
    photo = ImageTk.PhotoImage(img_resized)

    imgLBL = Label(master, image=photo)
    imgLBL.image = photo


    call_result1 = partial(generate_dxf,master,imgLBL,no_slots, no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distance1 , duct_radius1 , duct_angle1 ,duct_distance2 , duct_radius2 , duct_angle2, duct_distance3 , duct_radius3 , duct_angle3, notch_depth , notch_angle , no_turns)

    call_result2 = partial(generate_femm,no_slots , no_poles , stator_od , stator_id , rotor_od , shaft_dia , slot_opening , tooth_tip_angle , tooth_width , yoke_thickness , pole_length , pole_width , pole_arc , bridge_thickness , duct_distance1 , duct_radius1 , duct_angle1 ,duct_distance2 , duct_radius2 , duct_angle2, duct_distance3 , duct_radius3 , duct_angle3, notch_depth , notch_angle , no_turns)

    b1 = Button(master, text = "Generate DXF", command = call_result1)

    b2 = Button(master, text = "Run FEMM", command = call_result2)

    # arranging button widgets
    b1.grid(row = 36, column = 0, sticky = E)
    b2.grid(row = 36, column = 1, sticky = E)


    # setting image with the help of label

    imgLBL.grid(row = 1, column = 4, columnspan = 100, rowspan = 100, padx = 5, pady = 5,)



master = create_window()

GUI(master)
mainloop()
