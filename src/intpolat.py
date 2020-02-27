import matplotlib as mpl
#mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pyhdf.SD import SD,SDC
from scipy import interpolate
import math
import datetime
import sys

class nakajimaKingRetrieval:
        def getIdx(self, sza, vza, raa):
# the values of  sza, vza, raa must be the multiple of 3.
            a1=np.where(sza00==sza)[0]
            a2=np.where(vza00==vza)[0]
            a3=np.where(raa00==raa)[0]
            a4=a2[np.nonzero(np.in1d(a2,a1))[0]]
            idx = a3[np.nonzero(np.in1d(a3,a4))[0]]
            return idx

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

        def interplot3d(self, Vsza, Vvza, Vraa, arr):
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
#       1st facet
                z1 = np.zeros((2,2))
                z1.fill(np.nan)
                z1[0,0] = arr[idxo]
#        print idxr,arr[idxr]
                z1[0,1] = arr[idxr]
                z1[1,0] = arr[idxv]
                z1[1,1] = arr[idxvr]

                f1 = interpolate.interp2d(x, y, z1, kind='linear')
                znew1 = f1(Vvza,Vraa)

#       2rd facet
                z2 = np.zeros((2,2))
                z2.fill(np.nan)
                z2[0,0] = arr[idxs]
                z2[0,1] = arr[idxsr]
                z2[1,0] = arr[idxsv]
                z2[1,1] = arr[idxsvr]
                f2 = interpolate.interp2d(x, y, z2, kind='linear')
                znew2 = f2(Vvza,Vraa)

#     interpolate in SZA between 1st facet and 2rd facet
                xp = [int(Vsza/3)*3,int(Vsza/3)*3+3]
                fp = [znew1[0],znew2[0]]
                znew = np.interp(Vsza, xp, fp)
            else:
                znew = np.nan

            return znew

        def loadLUT(self, dirlut):
            add=dirlut
            hflut = SD(add,SDC.READ)
            #print hf00.info()

            u_obj = hflut.select('UpStokesU')
            u=u_obj.get()

            q_obj = hflut.select('UpStokesQ')
            q=q_obj.get()
            return u, q

        def loadLUTrflc(self, dirlut):
            add=dirlut
            hflut = SD(add,SDC.READ)
#    print hf00.info()

            i_obj = hflut.select('UpStokesI')
            i=i_obj.get()
#    sfcb_obj = hflut.select('UpSphericalAlbedo')
#    sfcb=sfcb_obj.get()
            return i#, sfcb[0]
        def retrieve(self, rflcArr, imodel, rflc02, rflc07):
        #     interpolate in tau from reflectivity
                rflc07new = np.zeros(14)
                for ireff in range(14):
                        xp = rflcArr[imodel, ireff, 0, :]
                        taunew = np.interp(rflc02, xp, tauList)
#       print ireff, taunew

                        x07 = rflcArr[imodel, ireff, 1, :]
                        rflc07new[ireff] = np.interp(taunew, tauList, x07)
#        print ireff, rflc07new

#print rflc07, rflc07new, reffList
#reffnew = np.interp(rflc07, rflc07new, reffList)
#print reffnew

                f = interpolate.interp1d(rflc07new, reffList)
#print f(rflc07)
                reffnew = f(rflc07)
#f2 = interpolate.interp1d(rflc07new, reffList, kind='cubic')
#print f2(rflc07)

                rflc02list = np.zeros(10)
                for itau in range(10):
                        xp2 = rflcArr[0, :, 0, itau]
                        rflc02list[itau] = np.interp(reffnew, reffList, xp2)

                ftau = interpolate.interp1d(rflc02list, tauList)
#print ftau(rflc02)
                tauValue = ftau(rflc02)
                #print 'Tau:', tauValue, 'Reff:', reffnew
                return tauValue, reffnew

# -----------------------------------------------

# -----------------------------------------------
# MAIN program
#
# parameters:
# Vsza: solar zenith angle
# Vvza: viewing zenith angle
# Vraa: relative azimuthal angle. !! Check with 0 degree definition carefully!!
# rflc02: reflectivity on obsoptive band, like MODIS band 02
# rflc07: reflectivity on emissivity band, like MODIS band 07
# -----------------------------------------------

#Vsza = 55.5 # load SZA
#Vvza = 60.5 # load VZA
#Vraa = 30.5 # load RAA

#rflc02 = 0.45
#rflc07 = 0.15

# you can wrap all parameters into routine interface
# for example:
# Vsza = float(sys.argv[1])
# Vvza = float(sys.argv[2])
# running command:
# $ python intpolat.py 55.5, 60.5
# or you can load all parameters using file
# for example:
# f_input=np.loadtxt('./f_input.txt',dtype='float')

print datetime.datetime.now()

# loading in input file

add1='./tr01_r050_c8_gm080_11.6339_b02_r0.00000_c000.00000_adding_out.hdf'
hf00 = SD(add1,SDC.READ)
#print hf00.info()

global sza00
global vza00
global raa00

datasets_dic = hf00.datasets()
#for idx, sds in enumerate(datasets_dic.keys()):
#    print idx, sds

i00_obj = hf00.select('UpStokesI')
i00 = i00_obj.get()

u00_obj = hf00.select('UpStokesU')
u00 = u00_obj.get()

q00_obj = hf00.select('UpStokesQ')
q00 = q00_obj.get()

sza_obj = hf00.select('SolarZenithAngle')
sza00 = sza_obj.get()
#print np.shape(sza00)

vza_obj = hf00.select('ViewingZenithAngle')
vza00 = vza_obj.get()

raa_obj = hf00.select('RelativeAzimuthAngle')
raa00 = raa_obj.get()

sphra_obj = hf00.select('UpSphericalAlbedo')
sphra00 = sphra_obj.get()

# -----------------------------------------------
# load in other LUTs, without angle info (angle info used from the above LUT)

# ------------------------------------------------------------------------

tauFilenameList     = ['000.00000','000.10000','000.40000','001.00000','002.50000','006.30000','010.00000','030.00000','070.00000','150.00000']
reffFilenameListTHM = ['0.57199','1.14818','1.75205','2.38922','3.05789','3.76296','4.50792','5.30667','6.16103','7.08433','9.24419','12.0296','16.1004','23.7867']
reffFilenameListMC6 = ['1.66200','3.32397','4.98597','6.64798','8.30994','9.97198','11.6339','13.2959','14.9579','16.6200','19.9439','23.2679','26.5919','29.9158']
tauList             = [0.0, 0.1, 0.4, 1.0, 2.5, 6.3, 10.0, 30.0, 70.0, 150.0]
reffList            = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0]
modelList           = ['Ar01/out/ar01_thm1_gt', 'Tr01/out/tr01_r050_c8']
bandList            = ['02', '07']
addr1 = 'ar01_thm1_gt_gm080_0.57199_b02_r0.00000_c000.00000_adding_out.hdf'
addr2 = 'tr01_r050_c8_gm080_11.6339_b02_r0.00000_c000.00000_adding_out.hdf'
addr0 = '../LUT/'

rflcArr = np.zeros(shape=(2, 14, 2, 10))
# imodel: model; i2:reff; i3:band; i4:tau

retrieval = nakajimaKingRetrieval()


#Vsza = 55.5 # load SZA
#Vvza = 60.5 # load VZA
#Vraa = 30.5 # load RAA

#rflc02 = 0.45
#rflc07 = 0.15

f_input=np.loadtxt('./out/069496.054_inv.dat')#,dtype='float')

f0=open('f_retrievals.txt','w')
for i in range(len(f_input)):#len(f_input)):
        Vsza = float(f_input[i,7]) # load SZA
        Vvza = float(f_input[i,8]) # load VZA
        Vraa = float(f_input[i,9]) * -1 # load RAA

        rflc02 = float(f_input[i,22])
        rflc07 = float(f_input[i,23])
        print Vsza, Vvza, Vraa, rflc02, rflc07

        # loading LUTs for THM
        imodel = 0
        for i3 in range(2):
                for i2 in range(14):
                        for i4 in range(10):
                                fileName = modelList[imodel]+'_gm080_'+reffFilenameListTHM[i2]+'_b'+bandList[i3]+'_r0.00000_c'+tauFilenameList[i4]+'_adding_out.hdf'
                                #print addr0+fileName
                                i00 = retrieval.loadLUTrflc(addr0+fileName)
                                rflcArr[imodel, i2, i3, i4] = retrieval.interplot3d(Vsza,Vvza,Vraa,i00) # interpolate I

        # loading LUTs for MC6
        imodel = 1
        for i3 in range(2):
                for i2 in range(14):
                        for i4 in range(10):
                                fileName = modelList[imodel]+'_gm080_'+reffFilenameListMC6[i2]+'_b'+bandList[i3]+'_r0.00000_c'+tauFilenameList[i4]+'_adding_out.hdf'
                                #print addr0+fileName
                                i00 = retrieval.loadLUTrflc(addr0+fileName)
                                rflcArr[imodel, i2, i3, i4] = retrieval.interplot3d(Vsza,Vvza,Vraa,i00) # interpolate I

#print datetime.datetime.now(), 'loading LUT finished'


# -------------------------------------------------------------
# call retrieve function

        imodel = 0
        tau, reff = retrieval.retrieve(rflcArr, imodel, rflc02, rflc07)
        print 'THM: ', 'Tau: ', tau, 'Reff: ', reff
        f0.write(str(tau)+' '+str(reff)+' ')

        imodel = 1
        tau, reff = retrieval.retrieve(rflcArr, imodel, rflc02, rflc07)
        print 'MC6: ', 'Tau: ', tau, 'Reff: ', reff
        f0.write(str(tau)+' '+str(reff)+'\n')
        

#rflc07list = np.zeros(14)
#for ireff in range(14):
#        xp3 = rflcArr[0, ireff, 0, :]
#       #print tauValue, tauList, xp3
#        rflc07list[ireff] = np.interp(tauValue, tauList, xp3)
#
#print rflc07list, reffList
#freff = interpolate.interp1d(rflc07list, reffList)
#reffValue = freff(rflc07)
#print 'Tau:', tauValue, 'Reff:', reffValue

print datetime.datetime.now(), 'computation done'
