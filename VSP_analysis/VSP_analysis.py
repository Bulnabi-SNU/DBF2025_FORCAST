import openvsp as vsp
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import configparser
import pandas as pd
import contextlib
import io

vsp.VSPCheckSetup()


'''
코드 실행을 위해서 해야하는 것
1. airfoil data path 설정
2. custom_dir를 vsp3 저장 위치로 설정
'''


# airfoil data path
s9027_path = r"C:\Users\ksjsms\2025\AIAA\OpenVSP_ws\s9027.dat"
naca0008_path = r"C:\Users\ksjsms\2025\AIAA\OpenVSP_ws\naca0008.dat"
naca0009_path = r"C:\Users\ksjsms\2025\AIAA\OpenVSP_ws\naca0009.dat"

# Create necessary directories
custom_dir = r"C:\Users\ksjsms\2025\AIAA\OpenVSP_ws"
vsp3_dir = os.path.join(custom_dir, "vsp3")
analysis_dir = os.path.join(custom_dir, "analysis_results")
ect_dir = os.path.join(custom_dir, "ect")

# analysis parameter
vsp_file = os.path.join(vsp3_dir, "Mothership.vsp3")
alpha_start = -2.0   # Starting angle of attack (degrees)
alpha_end = 13      # Ending angle of attack (degrees)
alpha_step = 0.5     # Step size (degrees)
Re = 380000          # Reynolds number
Mach = 0             # Mach number (subsonic)

if not os.path.exists(vsp3_dir):
    os.makedirs(vsp3_dir)

if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)

if not os.path.exists(ect_dir):
    os.makedirs(ect_dir)



"""Mother ship Parameters"""
# Left/Right Boom Variable Parameters (mm)
boom_diameter = 20      #Change Boom Diameter
boom_length = 900.0    #Change Boom Length
boom_xoffset = 15      #Boom Offset along the x-axis
boom_yoffset = 285        #Boom Offset along the y-axis
boom_zoffset = 25        #Boom Offset along the z-axis   

# Center Fuselage Parameters (mm)
c_fuse_length = 750
c_fuse_width = 130
c_fuse_height = 130
c_fuse_xoffset = -400
c_fuse_yoffset = 0
c_fuse_zoffset = 20

# L/R Fuselage Parameters (mm)
lr_fuse_length = 360
lr_fuse_width = 90
lr_fuse_height = 90
lr_fuse_xoffset = -50
lr_fuse_yoffset = -285
lr_fuse_zoffset = -20

# Main Wing Parameters (mm)
mainwing_span = 1700
mainwing_c_root = 330
mainwing_c_tip = 90
mainwing_sweep = 0
mainwing_xoffset = 0
mainwing_yoffset = 0
mainwing_zoffset = 0
mainwing_option_tip = 3 #Edge option
mainwing_length_tip = 10
mainwing_offset_tip = 0

# Horizonal Tail Wing Parameters (mm)
tailwing_span = 590
tailwing_c_root = 140
tailwing_c_tip = 140
tailwing_sweep = 0
tailwing_xoffset = 810
tailwing_yoffset = 0
tailwing_zoffset = 25
tailwing_option_tip = 3
tailwing_length_tip = 10
tailwing_offset_tip = 0

# Vertical Wing Parameters (mm)
verwing_span = 70
verwing_c_root = 140
verwing_c_tip = 140
verwing_sweep = 0
verwing_xoffset = 810
verwing_yoffset = 285
verwing_zoffset = 25
verwing_xRotate = 90
verwing_option_tip = 3
verwing_length_tip = 10
verwing_offset_tip = 0

# Flaperon Parameters
aileron_start = 0.4
aileron_end = 0.9
aileron_c_ratio = 0.35

# Rudder Parameters
rudder_start = 0.1
rudder_end = 1
rudder_c_ratio = 0.35

# Elevator Parameters
elevator_start = 0.2
elevator_end = 0.9
elevator_c_ratio = 0.35




"""Create Wing"""
# Clear the current model
vsp.ClearVSPModel()

### Create the main wing ###
wing_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(wing_id,"Main Wing")

# Airfoil Selection
vsp.ChangeXSecShape(vsp.GetXSecSurf(wing_id,0),0,vsp.XS_FILE_AIRFOIL)
xsec_1 = vsp.GetXSec(vsp.GetXSecSurf(wing_id,0),0)
vsp.ReadFileAirfoil(xsec_1,s9027_path)

# Set the dimensions for the main wing
vsp.SetParmVal(wing_id, "Span", "XSec_1", mainwing_span / 2)  # Span of the each wing (Half of span)
vsp.SetParmVal(wing_id, "Root_Chord", "XSec_1", mainwing_c_root)  # Root chord length
vsp.SetParmVal(wing_id, "Tip_Chord", "XSec_1", mainwing_c_tip)  # Tip chord length
vsp.SetParmVal(wing_id, "Sweep", "XSec_1", mainwing_sweep) #Sweep Angle
vsp.SetParmVal(wing_id, "Sweep_Location", "XSec_1", 0.5) #Sweep location
vsp.SetParmVal(wing_id, "X_Rel_Location", "XForm", mainwing_xoffset)  # Position along X-axis
vsp.SetParmVal(wing_id, "Y_Rel_Location", "XForm", mainwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(wing_id, "Z_Rel_Location", "XForm", mainwing_zoffset)  # Position vertically
vsp.SetParmVal(wing_id, "CapUMaxOption", "EndCap" , mainwing_option_tip)
vsp.SetParmVal(wing_id, "CapUMaxLength", "EndCap" , mainwing_length_tip)
vsp.SetParmVal(wing_id, "CapUMaxOffset", "EndCap" , mainwing_offset_tip)

### Create Aileron ###
aileron_id = vsp.AddSubSurf(wing_id,vsp.SS_CONTROL)

vsp.SetParmVal(wing_id,"EtaFlag","SS_Control_1",1)
vsp.SetParmVal(wing_id,"EtaStart","SS_Control_1",aileron_start)
vsp.SetParmVal(wing_id,"EtaEnd","SS_Control_1",aileron_end)
vsp.SetParmVal(wing_id,"Length_C_Start","SS_Control_1",aileron_c_ratio)
vsp.Update()

### Create the horizonal tail wing ###
tailwing_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(tailwing_id,"Tail Wing")

# Airfoil Selection
vsp.ChangeXSecShape(vsp.GetXSecSurf(tailwing_id,0),0,vsp.XS_FILE_AIRFOIL)
xsec_2 = vsp.GetXSec(vsp.GetXSecSurf(tailwing_id,0),0)
vsp.ReadFileAirfoil(xsec_1,naca0008_path)

# Set the dimensions for the tail wing
vsp.SetParmVal(tailwing_id, "Span", "XSec_1", tailwing_span / 2)  # Span of the each wing (Half of span)
vsp.SetParmVal(tailwing_id, "Root_Chord", "XSec_1", tailwing_c_root)  # Root chord length
vsp.SetParmVal(tailwing_id, "Tip_Chord", "XSec_1", tailwing_c_tip)  # Tip chord length
vsp.SetParmVal(tailwing_id, "Sweep", "XSec_1", tailwing_sweep) #Sweep Angle
vsp.SetParmVal(tailwing_id, "X_Rel_Location", "XForm", tailwing_xoffset)  # Position along X-axis
vsp.SetParmVal(tailwing_id, "Y_Rel_Location", "XForm", tailwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(tailwing_id, "Z_Rel_Location", "XForm", tailwing_zoffset)  # Position vertically
vsp.SetParmVal(tailwing_id, "CapUMaxOption", "EndCap" , tailwing_option_tip)
vsp.SetParmVal(tailwing_id, "CapUMaxLength", "EndCap" , tailwing_length_tip)
vsp.SetParmVal(tailwing_id, "CapUMaxOffset", "EndCap" , tailwing_offset_tip)

# Create Elevator
elevator_id = vsp.AddSubSurf(tailwing_id,vsp.SS_CONTROL)

vsp.SetParmVal(tailwing_id,"EtaFlag","SS_Control_1",1)
vsp.SetParmVal(tailwing_id,"EtaStart","SS_Control_1",elevator_start)
vsp.SetParmVal(tailwing_id,"EtaEnd","SS_Control_1",elevator_end)
vsp.SetParmVal(tailwing_id,"Length_C_Start","SS_Control_1",elevator_c_ratio)
vsp.Update()

### Create the vertical wing (Right) ###
verwing_right_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(verwing_right_id,"Vertical Wing Right")

# Airfoil Selection
vsp.ChangeXSecShape(vsp.GetXSecSurf(verwing_right_id,0),0,vsp.XS_FILE_AIRFOIL)
xsec_2 = vsp.GetXSec(vsp.GetXSecSurf(verwing_right_id,0),0)
vsp.ReadFileAirfoil(xsec_1,naca0008_path)

# Set the dimensions for the vertial tail wing (Right)
vsp.SetParmVal(verwing_right_id, "Span", "XSec_1", verwing_span)  # Span of the wing
vsp.SetParmVal(verwing_right_id, "Root_Chord", "XSec_1", verwing_c_root)  # Root chord length
vsp.SetParmVal(verwing_right_id, "Tip_Chord", "XSec_1", verwing_c_tip)  # Tip chord length
vsp.SetParmVal(verwing_right_id, "Sweep", "XSec_1", verwing_sweep) #Sweep Angle
vsp.SetParmVal(verwing_right_id, "X_Rel_Location", "XForm", verwing_xoffset)  # Position along X-axis
vsp.SetParmVal(verwing_right_id, "Y_Rel_Location", "XForm", verwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(verwing_right_id, "Z_Rel_Location", "XForm", verwing_zoffset)  # Position vertically
vsp.SetParmVal(verwing_right_id, "X_Rel_Rotation", "XForm", verwing_xRotate)  # X-axis Rotation
vsp.SetParmVal(verwing_right_id, "CapUMaxOption", "EndCap" , verwing_option_tip)
vsp.SetParmVal(verwing_right_id, "CapUMaxLength", "EndCap" , verwing_length_tip)
vsp.SetParmVal(verwing_right_id, "CapUMaxOffset", "EndCap" , verwing_offset_tip)
vsp.SetParmVal(verwing_right_id, "Sym_Planar_Flag","Sym", 0)

# Create Rudder (Right)
rudder_right_id = vsp.AddSubSurf(verwing_right_id,vsp.SS_CONTROL)

vsp.SetParmVal(verwing_right_id,"EtaFlag","SS_Control_1",1)
vsp.SetParmVal(verwing_right_id,"EtaStart","SS_Control_1",rudder_start)
vsp.SetParmVal(verwing_right_id,"EtaEnd","SS_Control_1",rudder_end)
vsp.SetParmVal(verwing_right_id,"Length_C_Start","SS_Control_1",rudder_c_ratio)
vsp.Update()


### Create the vertical wing (Left) ###
verwing_left_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(verwing_left_id,"Vertical Wing Left")

# Airfoil Selection
vsp.ChangeXSecShape(vsp.GetXSecSurf(verwing_left_id,0),0,vsp.XS_FILE_AIRFOIL)
xsec_3 = vsp.GetXSec(vsp.GetXSecSurf(verwing_left_id,0),0)
vsp.ReadFileAirfoil(xsec_1,naca0008_path)

# Set the dimensions for the vertial tail wing (Left)
vsp.SetParmVal(verwing_left_id, "Span", "XSec_1", verwing_span)  # Span of the wing
vsp.SetParmVal(verwing_left_id, "Root_Chord", "XSec_1", verwing_c_root)  # Root chord length
vsp.SetParmVal(verwing_left_id, "Tip_Chord", "XSec_1", verwing_c_tip)  # Tip chord length
vsp.SetParmVal(verwing_left_id, "Sweep", "XSec_1", verwing_sweep) #Sweep Angle
vsp.SetParmVal(verwing_left_id, "X_Rel_Location", "XForm", verwing_xoffset)  # Position along X-axis
vsp.SetParmVal(verwing_left_id, "Y_Rel_Location", "XForm", -1 * verwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(verwing_left_id, "Z_Rel_Location", "XForm", verwing_zoffset)  # Position vertically
vsp.SetParmVal(verwing_left_id, "X_Rel_Rotation", "XForm", verwing_xRotate)  # X-axis Rotation
vsp.SetParmVal(verwing_left_id, "CapUMaxOption", "EndCap" , verwing_option_tip)
vsp.SetParmVal(verwing_left_id, "CapUMaxLength", "EndCap" , verwing_length_tip)
vsp.SetParmVal(verwing_left_id, "CapUMaxOffset", "EndCap" , verwing_offset_tip)
vsp.SetParmVal(verwing_left_id, "Sym_Planar_Flag","Sym", 0)

# Create Rudder (Left)
rudder_left_id = vsp.AddSubSurf(verwing_left_id,vsp.SS_CONTROL)

vsp.SetParmVal(verwing_left_id,"EtaFlag","SS_Control_1",1)
vsp.SetParmVal(verwing_left_id,"EtaStart","SS_Control_1",rudder_start)
vsp.SetParmVal(verwing_left_id,"EtaEnd","SS_Control_1",rudder_end)
vsp.SetParmVal(verwing_left_id,"Length_C_Start","SS_Control_1",rudder_c_ratio)
vsp.Update()

# Save the generated VSP files in the "vsp3" directory
vsp.WriteVSPFile(os.path.join(vsp3_dir, "wing.vsp3"))
print("Wing Created successfully!")


"""Create Fuselage"""

# Delete Wing
vsp.VSPRenew()

### Create left boom ###
left_boom_id = vsp.AddGeom("FUSELAGE", "")
vsp.SetGeomName(left_boom_id,"Left Boom")

# Set the dimensions for the left boom
vsp.SetParmVal(left_boom_id, "Length", "Design", boom_length)  # Length of the boom
vsp.SetParmVal(left_boom_id, "X_Rel_Location", "XForm", boom_xoffset)  # Position along X-axis
vsp.SetParmVal(left_boom_id, "Y_Rel_Location", "XForm", boom_yoffset)  # Offset to the Y-axis
vsp.SetParmVal(left_boom_id, "Z_Rel_Location", "XForm", boom_zoffset)  # Offset to the Z-axis
vsp.Update()

# Get the number of cross-section surfaces for the left boom
left_num_xsec_surfs = vsp.GetNumXSecSurfs(left_boom_id)

# Loop through all cross-section surfaces
for i in range(left_num_xsec_surfs):
    left_xsec_surf = vsp.GetXSecSurf(left_boom_id, i)  # Get the XSecSurf ID for this surface
    left_num_xsecs = vsp.GetNumXSec(left_xsec_surf)  # Get the number of cross-sections in this surface

    # Loop through all cross-sections in this surface
    for j in range(left_num_xsecs):
        vsp.ChangeXSecShape(left_xsec_surf, j, vsp.XS_CIRCLE)  # Change the shape to a circle
        left_xsec = vsp.GetXSec(left_xsec_surf, j)  # Get the specific cross-section
        left_dia = vsp.GetXSecParm(left_xsec, "Circle_Diameter")  # Get the parameter for the diameter
        vsp.SetParmVal(left_dia, boom_diameter)  # Set the diameter value to 12.0
        vsp.ResetXSecSkinParms(left_xsec) #Reset all skinning options to 0 and OFF

### Create right boom ###
right_boom_id = vsp.AddGeom("FUSELAGE", "")
vsp.SetGeomName(right_boom_id,"Right Boom")

# Set the dimensions for the right boom
vsp.SetParmVal(right_boom_id, "Length", "Design", boom_length)  # Length of the boom
vsp.SetParmVal(right_boom_id, "X_Rel_Location", "XForm", boom_xoffset)  # Position along X-axis
vsp.SetParmVal(right_boom_id, "Y_Rel_Location", "XForm", -1 * boom_yoffset)  # Offset to the Y-axis
vsp.SetParmVal(right_boom_id, "Z_Rel_Location", "XForm", boom_zoffset)  # Offset to the Z-axis
vsp.Update()

# Get the number of cross-section surfaces for the right boom
right_num_xsec_surfs = vsp.GetNumXSecSurfs(right_boom_id)
vsp.SetGeomName(right_boom_id,"Right Boom")

# Loop through all cross-section surfaces
for i in range(right_num_xsec_surfs):
    right_xsec_surf = vsp.GetXSecSurf(right_boom_id, i)  # Get the XSecSurf ID for this surface
    right_num_xsecs = vsp.GetNumXSec(right_xsec_surf)  # Get the number of cross-sections in this surface

    # Loop through all cross-sections in this surface
    for j in range(right_num_xsecs):
        vsp.ChangeXSecShape(right_xsec_surf, j, vsp.XS_CIRCLE)  # Change the shape to a circle
        right_xsec = vsp.GetXSec(right_xsec_surf, j)  # Get the specific cross-section
        right_dia = vsp.GetXSecParm(right_xsec, "Circle_Diameter")  # Get the parameter for the diameter
        vsp.SetParmVal(right_dia, boom_diameter)  # Set the diameter value to 12.0
        vsp.ResetXSecSkinParms(right_xsec) #Reset all skinning options to 0 and OFF

### Create the left motor mount ####
left_fuse = vsp.AddGeom("FUSELAGE", "")
vsp.SetGeomName(left_fuse,"Left Fuselage")

# Set the dimensions for the fuselage
vsp.SetParmVal(left_fuse, "Length", "Design", lr_fuse_length)  # Length of the boom
vsp.SetParmVal(left_fuse, "X_Rel_Location", "XForm", lr_fuse_xoffset)  # Position along X-axis
vsp.SetParmVal(left_fuse, "Y_Rel_Location", "XForm", lr_fuse_yoffset)  # Offset to the left along Y-axis
vsp.SetParmVal(left_fuse, "Z_Rel_Location", "XForm", lr_fuse_zoffset)  # Offset to the left along Y-axis
vsp.Update()

# Get the number of cross-section surfaces for the left fuselage
leftfuse_num_xsec_surfs = vsp.GetNumXSecSurfs(left_fuse)

# Loop through all cross-section surfaces
for i in range(leftfuse_num_xsec_surfs):
    leftfuse_xsec_surf = vsp.GetXSecSurf(left_fuse, i)  # Get the XSecSurf ID for this surface
    leftfuse_num_xsecs = vsp.GetNumXSec(leftfuse_xsec_surf)  # Get the number of cross-sections in this surface

    # Loop through all cross-sections in this surface
    for j in range(1,leftfuse_num_xsecs-1):
        vsp.ChangeXSecShape(leftfuse_xsec_surf, j, vsp.XS_EDIT_CURVE)  # Change the shape to a rounded rectangle
        leftfuse_xsec = vsp.GetXSec(leftfuse_xsec_surf, j)  # Get the specific cross-section
        leftfuse_wid=vsp.GetXSecParm(leftfuse_xsec, "Width") # Get the parameter for the width
        leftfuse_hei=vsp.GetXSecParm(leftfuse_xsec, "Height") # Get the parameter for the height
        vsp.SetParmVal(leftfuse_wid, lr_fuse_width) 
        vsp.SetParmVal(leftfuse_hei, lr_fuse_height)
       
### Create the right motor mount ###
right_fuse = vsp.AddGeom("FUSELAGE", "")
vsp.SetGeomName(right_fuse,"Right Fuselage")

# Set the dimensions for the fuselage
vsp.SetParmVal(right_fuse, "Length", "Design", lr_fuse_length)  # Length of the boom
vsp.SetParmVal(right_fuse, "X_Rel_Location", "XForm", lr_fuse_xoffset)  # Position along X-axis
vsp.SetParmVal(right_fuse, "Y_Rel_Location", "XForm", -1 * lr_fuse_yoffset)  # Offset to the right along Y-axis
vsp.SetParmVal(right_fuse, "Z_Rel_Location", "XForm", lr_fuse_zoffset)  # Offset to the right along Y-axis
vsp.Update()

# Get the number of cross-section surfaces for the right fuselage
rightfuse_num_xsec_surfs = vsp.GetNumXSecSurfs(right_fuse)

# Loop through all cross-section surfaces
for i in range(rightfuse_num_xsec_surfs):
    rightfuse_xsec_surf = vsp.GetXSecSurf(right_fuse, i)  # Get the XSecSurf ID for this surface
    rightfuse_num_xsecs = vsp.GetNumXSec(rightfuse_xsec_surf)  # Get the number of cross-sections in this surface

    # Loop through all cross-sections in this surface
    for j in range(1,rightfuse_num_xsecs-1):
        vsp.ChangeXSecShape(rightfuse_xsec_surf, j, vsp.XS_EDIT_CURVE)  # Change the shape to a rounded rectangle
        rightfuse_xsec = vsp.GetXSec(rightfuse_xsec_surf, j)  # Get the specific cross-section
        rightfuse_wid=vsp.GetXSecParm(rightfuse_xsec, "Width") # Get the parameter for the width
        rightfuse_hei=vsp.GetXSecParm(rightfuse_xsec, "Height") # Get the parameter for the height
        vsp.SetParmVal(rightfuse_wid, lr_fuse_width) 
        vsp.SetParmVal(rightfuse_hei, lr_fuse_height)

### Create the center fuselage ###
c_fuse = vsp.AddGeom("FUSELAGE", "")
vsp.SetGeomName(c_fuse,"Center Fuselage")

# Set the dimensions for the fuselage
vsp.SetParmVal(c_fuse, "Length", "Design", c_fuse_length)  # Length of the boom
vsp.SetParmVal(c_fuse, "X_Rel_Location", "XForm", c_fuse_xoffset)  # Position along X-axis
vsp.SetParmVal(c_fuse, "Y_Rel_Location", "XForm", c_fuse_yoffset)  # Offset to the c along Y-axis
vsp.SetParmVal(c_fuse, "Z_Rel_Location", "XForm", c_fuse_zoffset)  # Offset to the c along Y-axis
vsp.Update()

# Get the number of cross-section surfaces for the c fuselage
cfuse_num_xsec_surfs = vsp.GetNumXSecSurfs(c_fuse)

# Loop through all cross-section surfaces
for i in range(cfuse_num_xsec_surfs):
    cfuse_xsec_surf = vsp.GetXSecSurf(c_fuse, i)  # Get the XSecSurf ID for this surface
    cfuse_num_xsecs = vsp.GetNumXSec(cfuse_xsec_surf)  # Get the number of cross-sections in this surface

    # Loop through all cross-sections in this surface
    for j in range(1,cfuse_num_xsecs-1):
        vsp.ChangeXSecShape(cfuse_xsec_surf, j, vsp.XS_EDIT_CURVE)  # Change the shape to a rounded rectangle
        cfuse_xsec = vsp.GetXSec(cfuse_xsec_surf, j)  # Get the specific cross-section
        cfuse_wid=vsp.GetXSecParm(cfuse_xsec, "Width") # Get the parameter for the width
        cfuse_hei=vsp.GetXSecParm(cfuse_xsec, "Height") # Get the parameter for the height
        vsp.SetParmVal(cfuse_wid, c_fuse_width) 
        vsp.SetParmVal(cfuse_hei, c_fuse_height)

# Save fuselageassemble and Mothership files similarly
vsp.WriteVSPFile(os.path.join(vsp3_dir, "fuselage.vsp3"))
print("Fuselage Created successfully!")


"""Create wing again"""
### Create the main wing ###
wing_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(wing_id,"Main Wing")

# Airfoil Selection
vsp.ChangeXSecShape(vsp.GetXSecSurf(wing_id,0),0,vsp.XS_FILE_AIRFOIL)
xsec_1 = vsp.GetXSec(vsp.GetXSecSurf(wing_id,0),0)
vsp.ReadFileAirfoil(xsec_1,s9027_path)

# Set the dimensions for the main wing
vsp.SetParmVal(wing_id, "Span", "XSec_1", mainwing_span / 2)  # Span of the each wing (Half of span)
vsp.SetParmVal(wing_id, "Root_Chord", "XSec_1", mainwing_c_root)  # Root chord length
vsp.SetParmVal(wing_id, "Tip_Chord", "XSec_1", mainwing_c_tip)  # Tip chord length
vsp.SetParmVal(wing_id, "Sweep", "XSec_1", mainwing_sweep) #Sweep Angle
vsp.SetParmVal(wing_id, "Sweep_Location", "XSec_1", 0.5) #Sweep location
vsp.SetParmVal(wing_id, "X_Rel_Location", "XForm", mainwing_xoffset)  # Position along X-axis
vsp.SetParmVal(wing_id, "Y_Rel_Location", "XForm", mainwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(wing_id, "Z_Rel_Location", "XForm", mainwing_zoffset)  # Position vertically
vsp.SetParmVal(wing_id, "CapUMaxOption", "EndCap" , mainwing_option_tip)
vsp.SetParmVal(wing_id, "CapUMaxLength", "EndCap" , mainwing_length_tip)
vsp.SetParmVal(wing_id, "CapUMaxOffset", "EndCap" , mainwing_offset_tip)

# Create Aileron
aileron_id = vsp.AddSubSurf(wing_id,vsp.SS_CONTROL)

vsp.SetParmVal(wing_id,"EtaFlag","SS_Control_1",1)
vsp.SetParmVal(wing_id,"EtaStart","SS_Control_1",aileron_start)
vsp.SetParmVal(wing_id,"EtaEnd","SS_Control_1",aileron_end)
vsp.SetParmVal(wing_id,"Length_C_Start","SS_Control_1",aileron_c_ratio)
vsp.Update()

### Create the horizonal tail wing ###
tailwing_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(tailwing_id,"Tail Wing")

# Set the dimensions for the tail wing
vsp.SetParmVal(tailwing_id, "Span", "XSec_1", tailwing_span / 2)  # Span of the each wing (Half of span)
vsp.SetParmVal(tailwing_id, "Root_Chord", "XSec_1", tailwing_c_root)  # Root chord length
vsp.SetParmVal(tailwing_id, "Tip_Chord", "XSec_1", tailwing_c_tip)  # Tip chord length
vsp.SetParmVal(tailwing_id, "Sweep", "XSec_1", tailwing_sweep) #Sweep Angle
vsp.SetParmVal(tailwing_id, "X_Rel_Location", "XForm", tailwing_xoffset)  # Position along X-axis
vsp.SetParmVal(tailwing_id, "Y_Rel_Location", "XForm", tailwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(tailwing_id, "Z_Rel_Location", "XForm", tailwing_zoffset)  # Position vertically
vsp.SetParmVal(tailwing_id, "CapUMaxOption", "EndCap" , tailwing_option_tip)
vsp.SetParmVal(tailwing_id, "CapUMaxLength", "EndCap" , tailwing_length_tip)
vsp.SetParmVal(tailwing_id, "CapUMaxOffset", "EndCap" , tailwing_offset_tip)

# Create Elevator
elevator_id = vsp.AddSubSurf(tailwing_id,vsp.SS_CONTROL)

vsp.SetParmVal(tailwing_id,"EtaFlag","SS_Control_1",1)
vsp.SetParmVal(tailwing_id,"EtaStart","SS_Control_1",elevator_start)
vsp.SetParmVal(tailwing_id,"EtaEnd","SS_Control_1",elevator_end)
vsp.SetParmVal(tailwing_id,"Length_C_Start","SS_Control_1",elevator_c_ratio)
vsp.Update()

### Create the vertical wing (Right) ###
verwing_right_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(verwing_right_id,"Vertical Wing Right")

# Set the dimensions for the vertial tail wing (Right)
vsp.SetParmVal(verwing_right_id, "Span", "XSec_1", verwing_span)  # Span of the wing
vsp.SetParmVal(verwing_right_id, "Root_Chord", "XSec_1", verwing_c_root)  # Root chord length
vsp.SetParmVal(verwing_right_id, "Tip_Chord", "XSec_1", verwing_c_tip)  # Tip chord length
vsp.SetParmVal(verwing_right_id, "Sweep", "XSec_1", verwing_sweep) #Sweep Angle
vsp.SetParmVal(verwing_right_id, "X_Rel_Location", "XForm", verwing_xoffset)  # Position along X-axis
vsp.SetParmVal(verwing_right_id, "Y_Rel_Location", "XForm", verwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(verwing_right_id, "Z_Rel_Location", "XForm", verwing_zoffset)  # Position vertically
vsp.SetParmVal(verwing_right_id, "X_Rel_Rotation", "XForm", verwing_xRotate)  # X-axis Rotation
vsp.SetParmVal(verwing_right_id, "CapUMaxOption", "EndCap" , verwing_option_tip)
vsp.SetParmVal(verwing_right_id, "CapUMaxLength", "EndCap" , verwing_length_tip)
vsp.SetParmVal(verwing_right_id, "CapUMaxOffset", "EndCap" , verwing_offset_tip)
vsp.SetParmVal(verwing_right_id, "Sym_Planar_Flag","Sym", 0)

# Create Rudder (Right)
rudder_right_id = vsp.AddSubSurf(verwing_right_id,vsp.SS_CONTROL)

vsp.SetParmVal(verwing_right_id,"EtaFlag","SS_Control_1",1)
vsp.SetParmVal(verwing_right_id,"EtaStart","SS_Control_1",rudder_start)
vsp.SetParmVal(verwing_right_id,"EtaEnd","SS_Control_1",rudder_end)
vsp.SetParmVal(verwing_right_id,"Length_C_Start","SS_Control_1",rudder_c_ratio)
vsp.Update()

### Create the vertical wing (Left) ###
verwing_left_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(verwing_left_id,"Vertical Wing Left")

# Set the dimensions for the vertial tail wing (Right)
vsp.SetParmVal(verwing_left_id, "Span", "XSec_1", verwing_span)  # Span of the wing
vsp.SetParmVal(verwing_left_id, "Root_Chord", "XSec_1", verwing_c_root)  # Root chord length
vsp.SetParmVal(verwing_left_id, "Tip_Chord", "XSec_1", verwing_c_tip)  # Tip chord length
vsp.SetParmVal(verwing_left_id, "Sweep", "XSec_1", verwing_sweep) #Sweep Angle
vsp.SetParmVal(verwing_left_id, "X_Rel_Location", "XForm", verwing_xoffset)  # Position along X-axis
vsp.SetParmVal(verwing_left_id, "Y_Rel_Location", "XForm", -1 * verwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(verwing_left_id, "Z_Rel_Location", "XForm", verwing_zoffset)  # Position vertically
vsp.SetParmVal(verwing_left_id, "X_Rel_Rotation", "XForm", verwing_xRotate)  # X-axis Rotation
vsp.SetParmVal(verwing_left_id, "CapUMaxOption", "EndCap" , verwing_option_tip)
vsp.SetParmVal(verwing_left_id, "CapUMaxLength", "EndCap" , verwing_length_tip)
vsp.SetParmVal(verwing_left_id, "CapUMaxOffset", "EndCap" , verwing_offset_tip)
vsp.SetParmVal(verwing_left_id, "Sym_Planar_Flag","Sym", 0)

# Create Rudder (Left)
rudder_left_id = vsp.AddSubSurf(verwing_left_id,vsp.SS_CONTROL)

vsp.SetParmVal(verwing_left_id,"EtaFlag","SS_Control_1",1)
vsp.SetParmVal(verwing_left_id,"EtaStart","SS_Control_1",rudder_start)
vsp.SetParmVal(verwing_left_id,"EtaEnd","SS_Control_1",rudder_end)
vsp.SetParmVal(verwing_left_id,"Length_C_Start","SS_Control_1",rudder_c_ratio)
vsp.Update()

# Save the Full Assemble model
vsp.WriteVSPFile(os.path.join(vsp3_dir, "Mothership.vsp3"))
print("Mothership Created successfully!")

################################################ Aerodynamic analysis ################################################
def calculate_coefficient(vsp_file, alpha_start, alpha_end, alpha_step, Re, Mach):
    """
    Calculates the coefficients for a model using VSPAERO.

    Args:
        vsp_file (str): Path to the .vsp3 model file.
        alpha_start (float): Starting angle of attack (degrees).
        alpha_end (float): Ending angle of attack (degrees).
        alpha_step (float): Step size for angle of attack (degrees).
        Re (float): Reynolds number.
        Mach (float): Freestream Mach number.

    Returns:
        dict: A dictionary with angles of attack as keys and corresponding coefficients as values.
    """
    # Clear previous data
    vsp.ClearVSPModel()
    vsp.VSPRenew()
    vsp.VSPCheckSetup()

    # Load the VSP model
    if not os.path.exists(vsp_file):
        raise FileNotFoundError(f"Model file {vsp_file} not found.")
    vsp.ReadVSPFile(vsp_file)

    # Configure VSPAERO
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetVSPAERORefWingID(wing_id) ##### wing_id 확인할것
    sweep_analysis = "VSPAEROSweep"

    # Compute geometry
    vsp.SetAnalysisInputDefaults(geom_analysis)
    vsp.ExecAnalysis(geom_analysis)

    # Configure sweep analysis for coefficient
    vsp.SetAnalysisInputDefaults(sweep_analysis)
    vsp.SetIntAnalysisInput(sweep_analysis, "AnalysisMethod", [vsp.VORTEX_LATTICE])
    vsp.SetIntAnalysisInput(sweep_analysis, "GeomSet", [vsp.SET_ALL])

    # Set the Number of points for alpha
    point_number = round(int((alpha_end - alpha_start) / alpha_step) + 1)

    # **Set the reference geometry set**
    vsp.SetDoubleAnalysisInput(sweep_analysis, "MachStart", [Mach])
    vsp.SetDoubleAnalysisInput(sweep_analysis, "ReCref", [Re])
    vsp.SetDoubleAnalysisInput(sweep_analysis, "AlphaStart", [alpha_start])
    vsp.SetDoubleAnalysisInput(sweep_analysis, "AlphaEnd", [alpha_end])
    vsp.SetIntAnalysisInput(sweep_analysis, "AlphaNpts", [point_number])

    vsp.Update()

    #Calculate Wing Surface Area
    Sref = vsp.GetParmVal(vsp.GetParm(wing_id,"TotalArea","WingGeom"))
   
    # Execute sweep analysis
    results_id = vsp.ExecAnalysis(sweep_analysis)

    # Move all files except .vsp3 from "vsp3" folder to "ect" folder
    for file in os.listdir(vsp3_dir):
        if not file.endswith(".vsp3"):
            src = os.path.join(vsp3_dir, file)
            dst = os.path.join(ect_dir, file)
        
            # 파일 존재 여부 확인 후 처리
            if os.path.exists(dst):
                print(f"File already exists at destination: {dst}. Overwriting...")
                os.remove(dst)  # 기존 파일 삭제
                
            os.rename(src, dst) # 파일 이동
            print(f"Moved {src} to {dst}")

    # Extract lift coefficient data
    sweepResults = vsp.GetStringResults(results_id, "ResultsVec")

    Alpha_res = [0] * point_number
    CL_res = [0] * point_number
    CDi_res = [None] * point_number
    CDo_res = [None] * point_number
    CDtot_res = [None] * point_number

    for i in range (point_number):
        alpha_vec = vsp.GetDoubleResults(sweepResults[i], "Alpha")
        Alpha_res[i] = alpha_vec[len(alpha_vec) - 1]

        cl_vec = vsp.GetDoubleResults(sweepResults[i], "CL")
        CL_res[i] = cl_vec[len(cl_vec) - 1]

        cdi_vec = vsp.GetDoubleResults(sweepResults[i], "CDi")
        CDi_res[i] = cdi_vec[len(cdi_vec) - 1]
        
        cdo_vec = vsp.GetDoubleResults(sweepResults[i], "CDo")
        CDo_res[i] = cdo_vec[len(cdo_vec) - 1]
        
        cdtot_vec = vsp.GetDoubleResults(sweepResults[i], "CDtot")
        CDtot_res[i] = cdtot_vec[len(cdtot_vec) - 1]

    return Alpha_res,CL_res,Sref, CDi_res, CDo_res, CDtot_res

# Calculate lift coefficients
alpha,CL,Sref,CDi,CDo,CDtot = calculate_coefficient(vsp_file, alpha_start, alpha_end, alpha_step, Re, Mach)

# Create a DataFrame for better handling and plotting
data = pd.DataFrame({
    'Alpha (deg)': alpha,
    'C_L': CL,
    'C_Di': CDi,
    'C_Do': CDo,
    'C_Dtot': CDtot
})

# Create CSV file
span = mainwing_span + mainwing_length_tip * 10
filename = f"aero_result_span_{span}.csv"
save_path = os.path.join(analysis_dir, filename)
data.to_csv(save_path, index=False)

# Display the DataFrame (optional)
print("\nAerodynamic Coefficients:")
print(data.to_string(index=False))

# Plot Coefficients (C_L, C_Di, C_Do, C_Dtot)
plt.figure(figsize=(12, 8))
plt.plot(
    data['Alpha (deg)'],
    data['C_L'],
    marker='o',
    linestyle='-',
    color='b',
    label='C_L (Lift Coefficient)'
)
plt.plot(
    data['Alpha (deg)'],
    data['C_Di'],
    marker='s',
    linestyle='--',
    color='r',
    label='C_Di (Induced Drag)'
)
plt.plot(
    data['Alpha (deg)'],
    data['C_Do'],
    marker='^',
    linestyle='-.',
    color='g',
    label='C_Do (Zero-Lift Drag)'
)
plt.plot(
    data['Alpha (deg)'],
    data['C_Dtot'],
    marker='d',
    linestyle=':',
    color='k',
    label='C_Dtot (Total Drag)'
)

# Adding Titles and Labels
plt.title('Aerodynamic Coefficients vs Angle of Attack', fontsize=16)
plt.xlabel('Angle of Attack (deg)', fontsize=14)
plt.ylabel('Coefficient Value', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()

