import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm
from pyhdf.SD import SD,SDC
from scipy import interpolate
import math
import datetime
import sys

class nakajimaKingRetrieval:
# Created on March 4 2019
# @author: yiwang
# This class works to retrieve optical thickness and particle size, 
# which are two key viariables to control the climate effect of clouds.
# the algorithms are developed by Nakajima and King, J ATMOS SCI, 1990
# the variable used in this class are depicted in the following figure.
#       idxsv________idxsvr
#        /  |         /|  interpolation 2d
#       /   |        / |   2rd facet
#      /    |       /  |
#   idxs____|___idxsr  |
#     |     |      |   |
#     |  idxv______|__idxvr
#  SZA|    /z[1,0] |   /z[1,1] 
#     |  x/        |  /
#     |  /VZA      | /   interpolation 2d
#     | /          |/     1st facet
#   idxo____y____idxr
#  z[0,0]  RAA   z[0,1]      
    
    def getIdx(self, sza, vza, raa):
            # return the location in the Look-up table regarding to the geometries
            # sza: the look-up tables of solar zenith angle
            # vza: the look-up tables of viewing zenith angle
            # raa: the look-up tables of relative azimuthal angle
            # idx: the index in the look-up table
        if not (sza and vza and raa):
            print('ERROR: look-up tables of SZA/VZA/RAA is null')
            return []
        a1=np.where(sza00==sza)[0]
        a2=np.where(vza00==vza)[0]
        a3=np.where(raa00==raa)[0]
        a4=a2[np.nonzero(np.in1d(a2,a1))[0]]
        idx = a3[np.nonzero(np.in1d(a3,a4))[0]]
        return idx

    def interplot3d(self, Vsza, Vvza, Vraa, arr):
            # a function for 3D interplotion
            # Vsza: input value of sza
            # Vvza: input value of vza
            # Vraa: input value of raa
            # znew: the index of geometry in the 3D look-up table
        if not (Vsza and Vvza and Vraa):
            print('ERROR: geometry info is null')
            return []
        idxo = self.getIdx(int(Vsza/3)*3,int(Vvza/3)*3,int(Vraa/3)*3)
        idxs = self.getIdx(int(Vsza/3)*3+3,int(Vvza/3)*3,int(Vraa/3)*3)
        idxv = self.getIdx(int(Vsza/3)*3,int(Vvza/3)*3+3,int(Vraa/3)*3)
        idxr = self.getIdx(int(Vsza/3)*3,int(Vvza/3)*3,int(Vraa/3)*3+3)
        idxsv = self.getIdx(int(Vsza/3)*3+3,int(Vvza/3)*3+3,int(Vraa/3)*3)
        idxvr = self.getIdx(int(Vsza/3)*3,int(Vvza/3)*3+3,int(Vraa/3)*3+3)
        idxsr = self.getIdx(int(Vsza/3)*3+3,int(Vvza/3)*3,int(Vraa/3)*3+3)
        idxsvr = self.getIdx(int(Vsza/3)*3+3,int(Vvza/3)*3+3,int(Vraa/3)*3+3)
        
        x = [int(Vvza/3)*3,int(Vvza/3)*3+3]
        y = [int(Vraa/3)*3,int(Vraa/3)*3+3]
        
        if (idxo and idxs and idxv and idxr and idxsv and idxvr and idxsr and idxsvr):
            # 1st facet
            z1 = np.zeros((2,2))
            z1.fill(np.nan)
            z1[0,0] = arr[idxo]
            z1[0,1] = arr[idxr]
            z1[1,0] = arr[idxv]
            z1[1,1] = arr[idxvr]
            
            f1 = interpolate.interp2d(x, y, z1, kind='linear')
            znew1 = f1(Vvza,Vraa)
            
            # 2rd facet
            z2 = np.zeros((2,2))
            z2.fill(np.nan)
            z2[0,0] = arr[idxs]
            z2[0,1] = arr[idxsr]
            z2[1,0] = arr[idxsv]
            z2[1,1] = arr[idxsvr]
            f2 = interpolate.interp2d(x, y, z2, kind='linear')
            znew2 = f2(Vvza,Vraa)
            
            # interpolate in SZA between 1st facet and 2rd facet
            xp = [int(Vsza/3)*3,int(Vsza/3)*3+3]
            fp = [znew1[0],znew2[0]]
            znew = np.interp(Vsza, xp, fp)
        else:
            znew = np.nan
        return znew

    def loadLUT(self, dirlut):
        # loading U and Q component from the LUT
        # dirlut: look-up table directory
        # u: U component in the polarization
        # q: Q component in the polarization
        if not dirlut:
            print('ERROR: LUT directory is not correct')
            return []
        add=dirlut
        hflut = SD(add,SDC.READ)
        
        u_obj = hflut.select('UpStokesU')
        u=u_obj.get()
        
        q_obj = hflut.select('UpStokesQ')
        q=q_obj.get()
        return u, q

    def loadLUTrflc(self, dirlut):
        # loading I component from the LUT
        # dirlut: look-up table directory
        # i: I component in the polarization
        if not dirlut:
            print('ERROR: LUT directory is not correct')
            return []
        add=dirlut
        hflut = SD(add,SDC.READ)
        
        i_obj = hflut.select('UpStokesI')
        i=i_obj.get()
        return i
    
    def retrieve(self, rflcArr, imodel, rflc02, rflc07):
        # interpolate in tau from reflectivity
        # rflcArr: matrix contains reflectivity
        # the index of ice model
        # rflc02: reflectivity at band 02
        # rflc07: reflectivity at band 07
        # tauValue: optical thickness
        # reffnew: particle size
        if not rflcArr:
            print('ERROR: reflectivity matrix is null')
            return []
        if not imodel:
            print('ERROR: please define the ice model index')
            return []
        if rflc02 and rflc07:
            rflc07new = np.zeros(14)
            for ireff in range(14):
                xp = rflcArr[imodel, ireff, 0, :]
                taunew = np.interp(rflc02, xp, tauList)
                x07 = rflcArr[imodel, ireff, 1, :]
                rflc07new[ireff] = np.interp(taunew, tauList, x07)
            f = interpolate.interp1d(rflc07new, reffList)
            reffnew = f(rflc07)
            rflc02list = np.zeros(10)
            for itau in range(10):
                xp2 = rflcArr[0, :, 0, itau]
                rflc02list[itau] = np.interp(reffnew, reffList, xp2)
                
            ftau = interpolate.interp1d(rflc02list, tauList)
            tauValue = ftau(rflc02)                
            return tauValue, reffnew
        else:
            print('ERROR: reflectivity at band 02/07 is not correct')
            return []

# ================================================
# MAIN program

# -----------------------------------------------
# input parameters
# dirdefault: default directory of look-up tables
# addr0: directory of LUT folder                
# f_input: file of input: lat, lon, sza, vza, raa, rflc02, rflc07 
# f0: file name of output: lat, lon, model 1 tau, model 1 particle size,  
#    model 2 tau, model 2 particle size

if sys.argv[1]:
    dirdefault = sys.argv[1]
else:
    dirdefault = './tr01_r050_c8_gm080_11.6339_b02_r0.00000_c000.00000_adding_out.hdf'
    
if sys.argv[2]:
    addr0 = sys.argv[2]
else:
    addr0 = '../LUT/'

if sys.argv[3]:
    f_input = np.loadtxt(sys.argv[3]))
else:
    f_input = np.loadtxt('./out/069496.054_inv.dat')#,dtype='float')

if sys.argv[4]:
    f0 = open(sys.argv[4],'w')
else:
    f0 = open('f_retrievals.txt','w')

print(datetime.datetime.now(), 'loading LUTs')
hf00 = SD(dirdefault,SDC.READ)
datasets_dic = hf00.datasets()

# set up the global variables and assign the values from the default LUT
global sza00
global vza00
global raa00

i00_obj = hf00.select('UpStokesI')
i00 = i00_obj.get()

u00_obj = hf00.select('UpStokesU')
u00 = u00_obj.get()

q00_obj = hf00.select('UpStokesQ')
q00 = q00_obj.get()

sza_obj = hf00.select('SolarZenithAngle')
sza00 = sza_obj.get()

vza_obj = hf00.select('ViewingZenithAngle')
vza00 = vza_obj.get()

raa_obj = hf00.select('RelativeAzimuthAngle')
raa00 = raa_obj.get()

sphra_obj = hf00.select('UpSphericalAlbedo')
sphra00 = sphra_obj.get()

# -----------------------------------------------
# setting parameters of the LUTs, without angle info (angle info used from the above LUT)
tauFilenameList     = ['000.00000','000.10000','000.40000','001.00000','002.50000',\
                       '006.30000','010.00000','030.00000','070.00000','150.00000']
reffFilenameListTHM = ['0.57199','1.14818','1.75205','2.38922','3.05789','3.76296',\
                       '4.50792','5.30667','6.16103','7.08433','9.24419','12.0296',\
                           '16.1004','23.7867']
reffFilenameListMC6 = ['1.66200','3.32397','4.98597','6.64798','8.30994','9.97198',\
                       '11.6339','13.2959','14.9579','16.6200','19.9439','23.2679',\
                           '26.5919','29.9158']
tauList             = [0.0, 0.1, 0.4, 1.0, 2.5, 6.3, 10.0, 30.0, 70.0, 150.0]
reffList            = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, \
                       60.0, 70.0, 80.0, 90.0]
modelList           = ['Ar01/out/ar01_thm1_gt', 'Tr01/out/tr01_r050_c8']
bandList            = ['02', '07']

rflcArr = np.zeros(shape=(2, 14, 2, 10)) # dim info: imodel: model; i2:reff; i3:band; i4:tau

# loading the retrieval class
retrieval = nakajimaKingRetrieval()

# -----------------------------------------------
# iteratedly loading LUTs
for i in range(len(f_input)): # loop for pixel
    Vsza = float(f_input[i,7]) # load SZA
    Vvza = float(f_input[i,8]) # load VZA
    Vraa = float(f_input[i,9]) * -1 # load RAA
    rflc02 = float(f_input[i,22]) # reflectivity at band 02
    rflc07 = float(f_input[i,23]) # reflectivity at band 07

    imodel = 0 # loading LUTs for THM
    for i3 in range(2): # loop for channel
        for i2 in range(14): # loop for reff
            for i4 in range(10): # loop for tau
                imodel = 0
                fileName = modelList[imodel]+'_gm080_'+reffFilenameListTHM[i2]\
                    +'_b'+bandList[i3]+'_r0.00000_c'+tauFilenameList[i4]+'_adding_out.hdf'
                i00 = retrieval.loadLUTrflc(addr0+fileName)
                rflcArr[imodel, i2, i3, i4] = retrieval.interplot3d\
                    (Vsza,Vvza,Vraa,i00) # interpolate I

                imodel = 1 # loading LUTs for MC6
                fileName = modelList[imodel]+'_gm080_'+reffFilenameListMC6[i2]\
                    +'_b'+bandList[i3]+'_r0.00000_c'+tauFilenameList[i4]+'_adding_out.hdf'
                i00 = retrieval.loadLUTrflc(addr0+fileName)
                rflcArr[imodel, i2, i3, i4] = retrieval.interplot3d(Vsza,Vvza,Vraa,i00) # interpolate I

    print(datetime.datetime.now(), 'loading LUT finished, retrieval starts ...')

    # -------------------------------------------------------------
    # call retrieve function
    f0.write(f_input[i,0]+' '+f_input[i,1]+' ') # lat, lon
    imodel = 0 # for model THM
    tau, reff = retrieval.retrieve(rflcArr, imodel, rflc02, rflc07)
    f0.write(str(tau)+' '+str(reff)+' ') # tau, reff of model 1

    imodel = 1 # for model MC6
    tau, reff = retrieval.retrieve(rflcArr, imodel, rflc02, rflc07)
    f0.write(str(tau)+' '+str(reff)+'\n') # tau, reff of model 2

f0.close()        
print(datetime.datetime.now(), 'computation is done, plotting starts')

# -----------------------------------------------------------
# plotting

if sys.argv[4]:
    f_data = np.loadtxt(sys.argv[4])
else:
    f_data = np.loadtxt('f_retrievals.txt')

for i in range(len(f_data)):
    lat0  = float(f_input[i,0]) # load lat
    lon0  = float(f_input[i,1]) # load lon
    tau1  = float(f_input[i,2]) # load tau of model 1
    reff1 = float(f_input[i,3]) # load reff of model 1
    tau2  = float(f_input[i,4]) # load tau of model 2
    reff2 = float(f_input[i,5]) # load reff of model 2

map = Basemap(projection='cea',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')

x, y = map(lon0,lat0)
# draw coastlines, state and country boundaries, edge of map.
map.drawcoastlines(linewidth=0.5)

hb = map.scatter(x,y,c=tau2/tau1,marker='h',s=3,lw=0,cmap=plt.cm.Spectral_r)

cbar=plt.colorbar(hb, orientation='horizontal',shrink=0.4)
cbar.set_label('The optical thickness difference of two ice models')
plt.clim(0,90)

plt.show()

pngname = "./glbtaudiff"+".pdf"
print("save ", pngname)
plt.savefig(pngname, dpi=150, facecolor='w', edgecolor='w',
    orientation='portrait', papertype=None, format=None,
    transparent=False, bbox_inches='tight', pad_inches=0.1)

print(datetime.datetime.now(),'plotting is done')
