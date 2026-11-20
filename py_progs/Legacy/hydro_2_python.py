#!/usr/bin/env python

'''
Parses a hydro input file into an astropy table.

Synopsis:
    This is a simple program which will
    parse a hydro input file into an astropy
    table that python can read in as a wind 
    

Description:  

    The code attempts to work out if the 
    input file is an hdf file or an ascii file.
    The hdf file is produced by default from
    zeus, and is much the easier option 
    however the code will have a try and work
    out how the r and theta files are named
    using the ndf filename. 
    In this initial incarnation, it only deals with
    r theta files, future iterations will deal
    with different dimensionalities, and coordinate
    types as and when it becomes necessary.

    There are two subroutines, which deal with 
    importing the data from either of the input 
    file formats, and returning an indentical
    dictionary that the main routine can use.
    

Arguments:  

    filename of the input file - if missed out code 
        will prompt for it

Returns:

    an astropy formatted ascii file, containing
    r, theta, the three components of velocity, 
    the temperature of each cell and the density.

Notes:
                                       
History:
15may    nsh    Coded


'''


import sys
import numpy as np
from pyhdf import SD
from astropy import constants as consts
from astropy.io import ascii
from astropy.table import Table


'''The following routine reads in hdf data from the file fname
   It takes a filename as input argument, and returns a dictionary
   The upper level of the dictionay contains
   Filename - the filename originally supplied
   Coord_sys - the coordinate system - NB, it only knows about spol at the moment
   Time - the time stamp of the simulation
   N_Data - the number of different data types
   Dims - the number of dimensions
   Data - This is in turn a dictionary that contains the different data
    This dictionary contains infividual dictionaries, each of which has
        data - an array of the data
        x1 - the first corrdinate (usually theta)
        x2 - the second coordinate (usually r)'''


def get_hdf_data(fname):


    hdf = SD.SD(fname)
    info= hdf.datasets()

    #Lets see what is inside

    hdf.info()
    data_sets=[]

    for name in sorted(info.keys()):
        if name[0:4]=="Data":
            sds=hdf.select(name)
            long_name=sds.attributes()["long_name"]
            for i in range(len(sds.attributes()["long_name"])):
                if long_name[i:i+2]=="AT":
                    junk=i-1
            short_name=long_name[:junk]
            data_sets.append([name,long_name,short_name])

#Get the time from the last long name 
    if long_name[junk+9:] != "********":
        time=float(long_name[junk+9:])
    else:
        time=0.0
        
#Get the dimensions from the last data set

    dims=len((sds.info()[2]))

#Get the coordinate system from the last data set

    coord_sys=sds.attributes()["coordsys"]
    if coord_sys=='spherical polar':
        coord_sys='spol'
    else:
        print(("get_hdf_data: I don't understand coordinate system ",coord_sys)) 
        exit()

#Now we know which of the datasets contain real data, we can extract all the data

    
    alldat={}
    alldat["Filename"]=fname    
    alldat["Coord_sys"]=coord_sys
    alldat["Time"]=time
    alldat["Data"]={}
    alldat["Data_names"]=np.array(data_sets)[:,2]
    alldat["N_data"]=len(data_sets)
    alldat["Dims"]=dims

#Loop over all the data sets in the hdf file - name each of the resulting dictionaries with the short name

    for i in range (len(data_sets)):
        print((data_sets[i][2]))
        sds=hdf.select(data_sets[i][0])
        data = sds.get()
        c1=info[data_sets[i][0]][0][0]
        c2=info[data_sets[i][0]][0][1]
        sds=hdf.select(c1)
        x2=sds.get()
        sds=hdf.select(c2)
        x1=sds.get()
        alldat[data_sets[i][2]]={}
        alldat[data_sets[i][2]]=data
    
    alldat["r_cent"]=x1
    alldat["theta_cent"]=x2

#HDF files only give us the centre of the grid, python needs the inner edges as well.
    
    r_edge=[]
    r_ratio=(x1[2]-x1[1])/(x1[1]-x1[0])
    dr=(x1[1]-x1[0])/(0.5*(1.0+r_ratio))
    r_edge.append(x1[0]-0.5*dr)
    print((r_edge[0],r_ratio))
    for i in range(len(x1)-1):
        r_edge.append(r_edge[-1]+dr)
        dr=dr*r_ratio

    theta_edge=[]
    theta_ratio=(x2[2]-x2[1])/(x2[1]-x2[0])
    dtheta=(x2[1]-x2[0])/(0.5*(1.0+theta_ratio))
    theta_min=x2[0]-0.5*dtheta
    if theta_min<0.0:
        theta_min=0.0
    theta_edge.append(theta_min)
    print((x2[0]))
    print((theta_edge[0],theta_ratio))
    for i in range(len(x2)-1):
        theta_edge.append(theta_edge[-1]+dtheta)
        dtheta=dtheta*theta_ratio
    if (theta_edge[-1]+(x2[-1]-theta_edge[-1])*2.0)>(np.pi/2.0):
        x2[-1]=(theta_edge[-1]+(np.pi/2.0))/2.0


    alldat["r_edge"]=r_edge
    alldat["theta_edge"]=theta_edge


    return(alldat)

'''The following routine reads in data from the file fname, radial
   grid data from r_file an theta grid data from theta_file
   It takes three filenames as input arguments, and returns a dictionary
   The upper level of the dictionay contains
   Filename - the filename originally supplied
   Coord_sys - the coordinate system - NB, it only knows about spol at the moment
   Time - the time stamp of the simulation
   N_Data - the number of different data types
   Dims - the number of dimensions
   Data - This is in turn a dictionary that contains the different data
    This dictionary contains infividual dictionaries, each of which has
        data - an array of the data
        x1 - the first corrdinate (usually theta)
        x2 - the second coordinate (usually r)'''



def get_ndf_data(fname,r_file,theta_file):




    inp=open(fname,"r")
    ir=[]
    itheta=[]
    data=[]

    raw_names=inp.readline()

    for line in inp.readlines():
        data_temp=line.split()
        ir.append(int(data_temp[0]))
        itheta.append(int(data_temp[1]))
        temp=[]
        for i in range(len(data_temp)-2):
            temp.append(float(data_temp[i+2]))
        data.append(temp)
    inp.close()
    
    r_cent=[]
    r_edge=[]
    theta_cent=[]
    theta_edge=[]
    
    inp=open(theta_file,"r")
    for line in inp.readlines():
        data_temp=line.split()
        try:
            if int(data_temp[0]) >=np.min(itheta) and int(data_temp[0]) <= np.max(itheta):
                theta_cent.append(float(data_temp[2]))
                theta_edge.append(float(data_temp[1]))
            if int(data_temp[0]) == np.max(itheta):
                break
        except:
            print(("Something wrong with theta data file ",theta_file))
    inp.close()

    inp=open(r_file,"r")
    for line in inp.readlines():
        data_temp=line.split()
        try:
            if int(data_temp[0]) >=np.min(ir) and int(data_temp[0]) <= np.max(ir):
                r_cent.append(float(data_temp[2]))
                r_edge.append(float(data_temp[1]))
            if int(data_temp[0]) == np.max(ir):
                break
        except:
            print(("Something wrong with r data file ",r_file))
    inp.close()

    data_sets=["DENSITY","1-VELOCITY","2-VELOCITY","3-VELOCITY","TOTAL ENERGY"]
    
    alldat={}
    alldat["Filename"]=fname    
    alldat["Coord_sys"]="spol"
    alldat["Time"]=0
    alldat["N_data"]=5
    alldat["Dims"]=2
    alldat["Data_names"]=data_sets
    

    alldat["theta_cent"]=np.array(theta_cent)
    alldat["theta_edge"]=np.array(theta_edge)
    alldat["r_cent"]=np.array(r_cent)
    alldat["r_edge"]=np.array(r_edge)


    for i in range (len(data_sets)):
        alldat[data_sets[i]]=np.reshape(np.array(data)[:,i],(len(theta_cent),len(r_cent)))
        
    
    return(alldat)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys



    #Start of the main code

    print ("Welcome to hydro_2_python") 
    MU=0.6

    if len(sys.argv)>1:
        fname=sys.argv[1]
    else:
        fname=input("Enter name of file to convert (main data file): ")



    if fname[0:3]=="ndf":
        print ("I think this is an ascii file")
        r_file="grid_r_"+fname[-3:]+".dat"
        theta_file="grid_theta_"+fname[-3:]+".dat"
        data=get_ndf_data(fname,r_file,theta_file)
    elif fname[0:3]=="hdf":
        print ("I think this is an hdf file")
        data=get_hdf_data(fname)

    #We should now have three components of velocity, internal energy, and density.
    #We need to compute the temperature in order to supply it to python.


    temp=(2.0/3.0)*data["TOTAL ENERGY"]
    data["TEMPERATURE"]=temp/((data["DENSITY"]/(consts.m_p.cgs*MU))*consts.k_B.cgs)

    # Open an output file 
    fname=data["Filename"]+".zeu"

    # Preamble

    out=open(fname,'w')
    out.write("# This is a file generated by hydro_to_python\n")
    out.write("# We can put any number of comments in behind # signs\n")
    out.write("# By default, the order of coordinates are \n")
    out.write("#                r, theta phi for spherical polars\n")
    out.write("#                         x,y,z        for carteisan\n")
    out.write("#                         or w, z, phi    for cylindrical\n")
    out.write("# Coordinate_system "+data["Coord_sys"]+"\n")
    out.write("# Dimensions "+str(data["Dims"])+"\n")



    ndata=len(data)
    data_names=list(data.keys())
    titles=[]
    titles=titles+["ir","r_cent","r_edge"]
    titles=titles+["itheta","theta_cent","theta_edge"]
    titles=titles+["v_r","v_theta","v_phi","density","temperature"]







    col0=np.array([])
    col1=np.array([])
    col2=np.array([])
    col3=np.array([])
    col4=np.array([])
    col5=np.array([])
    col6=np.array([])
    col7=np.array([])
    col8=np.array([])
    col9=np.array([])
    col10=np.array([])

    fmt='%013.6e'

    #This next line defines formats for the output variables. This is set in a dictionary
    fmts={    'ir':'%03i',    
    'r_cent':fmt,
    'r_edge':fmt,
    'itheta':'%i',    
    'theta_cent':fmt,
    'theta_edge':fmt,
    'iphi':'%03i',    
    'phi_cent':fmt,
    'phi_edge':fmt,
    'ix':'%03i',    
    'x_cent':fmt,
    'x_edge':fmt,
    'v_r':fmt,
    'v_theta':fmt,
    'v_phi':fmt,
    'density':fmt,
    'temperature':fmt}





    for j in range(len(data["theta_cent"])):
        col0=np.append(col0,np.arange(len(data["r_cent"])))
        col1=np.append(col1,data["r_cent"])
        col2=np.append(col2,data["r_edge"])
        col3=np.append(col3,np.ones(len(data["r_cent"]))*j)
        col4=np.append(col4,np.ones(len(data["r_cent"]))*data["theta_cent"][j])
        col5=np.append(col5,np.ones(len(data["r_cent"]))*data["theta_edge"][j])
        col6=np.append(col6,data["1-VELOCITY"][j])
        col7=np.append(col7,data["2-VELOCITY"][j])
        col8=np.append(col8,data["3-VELOCITY"][j])
        col9=np.append(col9,data["DENSITY"][j])
        col10=np.append(col10,data["TEMPERATURE"][j])

    out_dat=Table([col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10],names=titles)
    ascii.write(out_dat,out,formats=fmts)





    out.close()




    data=ascii.read(fname)

