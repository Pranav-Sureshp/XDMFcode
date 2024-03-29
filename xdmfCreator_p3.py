"""
AUTHOR
Ryan Stefan

DESCRIPTION
Read a series of hdf files in the local directory and create an .xdmf file 
that can be used to viewed them with ParaView.
Processes all ALMA results2d/3d and plt2d/3d files automatically that are
present in the local directory.
Determine which hdf files should be processed by moving them in and out of the
directory.

USAGE
python slice.py

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import glob

# =============================================================================
# FUNCTIONS
# =============================================================================
class H5DatasetNames:
    def __init__(self):
        self.dsNames = []
        self.dsSize = []
        self.dsShape = []

    def __call__(self, name, h5obj):
        # Only h5py datasets have dtype attribute, so this will skip the names
        # that have only partial paths, like /data, /data/var3d, etc.
        # hasattr checks if the class has that member name
        if hasattr(h5obj,'dtype') and not name in self.dsNames:
            self.dsNames += [name]      # Full path to dataset
            self.dsSize += [h5obj.size] # Number of elements in dataset
            self.dsShape += [h5obj.shape] # Number of elements in dataset
         # No return will make visititem (and visit) function recursive

# =============================================================================
# Write vtu file header
def vtuHeader(fp):
  fp.write('<?xml version="1.0"?>\n')
  fp.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
  fp.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
  fp.write('  <Domain Name="MSI">\n')
  fp.write('    <Grid Name="CellTime" GridType="Collection" CollectionType="Temporal">\n')

# =============================================================================
# Write vtu file footer
def vtuFooter(fp):
  fp.write('    </Grid>\n')
  fp.write('  </Domain>\n')
  fp.write('</Xdmf>\n')

# =============================================================================
# Write vtu PointData header
def vtuGridStart(fp,time,nCells):
  fp.write('      <Grid Name="rectilinear_scalar" GridType="Uniform">\n')
  fp.write('        <Time Value="%15.8e" />\n' % (time))
  fp.write('        <Topology TopologyType="3DRECTMesh" NumberOfElements="%d %d %d"/>\n' % (nCells[2],nCells[1],nCells[0]))

# =============================================================================
# Write vtu PointData footer
def vtuGridEnd(fp):
  fp.write('      </Grid>\n')
  
# =============================================================================
# Write vtu PointData object
def vtuGridAddGeom(fp,nDim,nCells,fileName,xGrpNames):
  fp.write('        <Geometry GeometryType="VXVYVZ">\n')
  fp.write('          <DataItem Name="x" Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">\n' % nCells[0])
  #fp.write('            plt3d_000000000.h5:/coords/x\n')
  fp.write('            %s:%s\n' % (fileName,xGrpNames[0]))
  fp.write('          </DataItem>\n')
  if (nDim>=2):
    fp.write('          <DataItem Name="y" Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">\n' % nCells[1])
    #fp.write('            plt3d_000000000.h5:/coords/y\n')
    fp.write('            %s:%s\n' % (fileName,xGrpNames[1]))
    fp.write('          </DataItem>\n')
  if (nDim>=3):
    fp.write('          <DataItem Name="z" Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">\n' % nCells[2])
    #fp.write('            plt3d_000000000.h5:/coords/z\n')
    fp.write('            %s:%s\n' % (fileName,xGrpNames[2]))
    fp.write('          </DataItem>\n')
  fp.write('        </Geometry>\n')

# =============================================================================
# Write vtu PointData object
def vtuGridAddData(fp,nDim,nCells,attrName,fileName,fGrpName):
  fp.write('        <Attribute Name="%s" AttributeType="Scalar" Center="Node">\n' % (attrName))
  if (nDim==1):
    fp.write('          <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">\n' % (nCells[0]))
    #fp.write('            plt3d_000000000.h5:/data/var1d/rh/000001\n')
    fp.write('            %s:%s\n' % (fileName,fGrpName))
  if (nDim==2):
    fp.write('          <DataItem Dimensions="%d %d" NumberType="Float" Precision="4" Format="HDF">\n' % (nCells[1],nCells[0]))
    #fp.write('            plt3d_000000000.h5:/data/var2d/rh/000001\n')
    fp.write('            %s:%s\n' % (fileName,fGrpName))
  if (nDim==3):
    fp.write('          <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="4" Format="HDF">\n' % (nCells[2],nCells[1],nCells[0]))
    #fp.write('            plt3d_000000000.h5:/data/var3d/rh/000001\n')
    fp.write('            %s:%s\n' % (fileName,fGrpName))
  fp.write('          </DataItem>\n')
  fp.write('        </Attribute>\n')

# =============================================================================
# High level function to wrtie out all the xdmf content.
def writePltFile(nDim,wkdir,filePrefix,fileExt,attrPos,grpNameTime,grpNameStep,grpNameCoords):
  nDim3 = 3 # Structure of 2d xdmf file assumes a 3d setup with 2d hdf dataset inside.
    
  # Get list of all h5 files in this directory
  dirListAbs = glob.glob(wkdir+filePrefix+'*.'+fileExt)
  nFiles = len(dirListAbs)
  print( 'Number of files: ', nFiles)
  if (nFiles==0):
    return
    
  dirListRel = [''] * nFiles # Allocate array of strings
  for ifi in range(0,nFiles):
    dirListRel[ifi] = os.path.basename(dirListAbs[ifi]) # File name only, remove path

  fp = open(wkdir+filePrefix+'new.xdmf','w')
  vtuHeader(fp)

  # Loop on h5 files
  for ifi in range(0,nFiles):
    print( '  Processing file: ', dirListRel[ifi])
    fh5 = h5py.File(dirListAbs[ifi],'r')
    h5Ds = H5DatasetNames()
    fh5.visititems(h5Ds) # Add dataset items to a list object, h5Ds.dsNames

    # Get coords. Assume all x,y,z are present even for 2d, where one dimension will be 1
    nCellsTopoTmp = np.zeros(nDim3,dtype=np.int) # Always use 3 coords, x,y,z. In 2d assume one component will be 1 cell thick
    for id in range(0,len(nCellsTopoTmp)):
      tmp = fh5.get(grpNameCoords[id])
      nCellsTopoTmp[id] = tmp.shape[0]
    nCellsTopo = tuple(nCellsTopoTmp) # nx,ny,nz

    # Get dataset
    dsGrpNames = []
    dsGrpSize = []
    dsGrpShape = []
    dsGrpAttrName = []
    for ids in range(0,len(h5Ds.dsShape)): # Loop on all group objects to separate 1d, 2d, 3d objects.     
      if (len(h5Ds.dsShape[ids])==nDim): # Look for 2d or 3d datasets only
        dsGrpNames += ['/'+h5Ds.dsNames[ids]]
        dsGrpSize += [h5Ds.dsSize[ids]]
        if (nDim==2):
          # For 2d we need to know which coordinate is only 1 cell thick, this info is not present in
          # the dataset. Use the coords in this case.
          dsGrpShape += [nCellsTopo[::-1]] # reverse order for hdf file, nz,ny,nx
        else:
          # For 3d we can take dimensions directly from the dataset.
          dsGrpShape += [h5Ds.dsShape[ids]]
        dsGrpAttrName += [h5Ds.dsNames[ids].split('/')[attrPos]]
    
    # Get time values
    tval = fh5.get(grpNameTime)
    time = np.array(tval)[0]
    print( '  Time: ', time)
    
    # Get step values
    sval = fh5.get(grpNameStep)
    step = int(round(np.array(sval)[0]))
    print ('  Step: ', step )
        
    fh5.close()
    
    vtuGridStart(fp,time,nCellsTopo)
    vtuGridAddGeom(fp,nDim3,nCellsTopo,dirListRel[ifi],grpNameCoords)

    # Loop through 3d datasets and write to xdmf file
    nDsGrp = len(dsGrpNames)
    for ids in range(0,nDsGrp):
      nCellsRev = dsGrpShape[ids] # tuple object holding shape arrays
      nCells = nCellsRev[::-1] # hdf dimensions are reversed
      if (nCellsTopo!=nCells):
        print ('  ERROR: Number of cells in topology not same as dataset',nCellsTopo,nCells)
        return
      vtuGridAddData(fp,nDim3,nCells,dsGrpAttrName[ids],dirListRel[ifi],dsGrpNames[ids])
    vtuGridEnd(fp)
    
  vtuFooter(fp)
  fp.close()

# =============================================================================
# High level function to wrtie out all the xdmf content.
# High level function to wrtie out all the xdmf content.
def writeResFile(nDim,wkdir,filePrefix,fileExt,attrPos,grpNameTime,grpNameStep,grpNameCoords):
  nDim3 = 3 # Structure of 2d xdmf file assumes a 3d setup with 2d hdf dataset inside.
  
  # Get list of all h5 files in this directory
  dirListAbs = glob.glob(wkdir+filePrefix+'*.'+fileExt)
  nFiles = len(dirListAbs)
  print ('Number of files: ', nFiles)
  if (nFiles==0):
    return
    
  dirListRel = [''] * nFiles # Allocate array of strings
  for ifi in range(0,nFiles):
    dirListRel[ifi] = os.path.basename(dirListAbs[ifi]) # File name only, remove path

  fp = open(wkdir+filePrefix+'new.xdmf','w')
  vtuHeader(fp)

  # Loop on h5 files
  for ifi in range(0,nFiles):
    print ('  Processing file: ', dirListRel[ifi])
    fh5 = h5py.File(dirListAbs[ifi],'r')
    h5Ds = H5DatasetNames()
    fh5.visititems(h5Ds) # Add dataset items to a list object, h5Ds.dsNames

    # Get coords. Assume all x,y,z are present even for 2d, where one dimension will be 1
    nCellsTopoTmp = np.zeros(nDim3,dtype=np.int)
    for id in range(0,len(nCellsTopoTmp)):
      tmp = fh5.get(grpNameCoords[id])
      nCellsTopoTmp[id] = tmp.shape[0]
    nCellsTopo = tuple(nCellsTopoTmp) # nx,ny,nz

    # Get datasets
    dsGrpNames = []
    dsGrpSize = []
    dsGrpShape = []
    dsGrpAttrName = []
    # Loop on all group objects to separate 1d and 3d objects.
    for ids in range(0,len(h5Ds.dsShape)):        
      # Save 3d datasets
      if (len(h5Ds.dsShape[ids])==nDim):
        dsGrpNames += ['/'+h5Ds.dsNames[ids]]
        dsGrpSize += [h5Ds.dsSize[ids]]
        if (nDim==2):
          # For 2d we need to know which coordinate is only 1 cell thick, this info is not present in
          # the dataset. Use the coords in this case.
          dsGrpShape += [nCellsTopo[::-1]] # reverse order for hdf file, nz,ny,nx
        else:
          # For 3d we can take dimensions directly from the dataset.
          dsGrpShape += [h5Ds.dsShape[ids]]
        dsGrpAttrName += [h5Ds.dsNames[ids].split('/')[attrPos]]
    
    # Get time values
    tval = fh5.get(grpNameTime)
    timeRes = np.array(tval)
    print ('  Times: ', timeRes)
    
    nTimes = len(timeRes)
    if (len(dsGrpNames)%nTimes!=0):
      print ('  ERROR: Number of group names inconsistent with number of times in file ',dirListRel[ifi])
      return
    nDsGrp = len(dsGrpNames)/nTimes
    
    # Get step values
    sval = fh5.get(grpNameStep)
    stepRes = np.array(sval).astype(int)
    print ('  Steps: ',stepRes)
        
    fh5.close()
    
    for itf in range(0,nTimes): 
      vtuGridStart(fp,timeRes[itf],nCellsTopo)
      vtuGridAddGeom(fp,nDim3,nCellsTopo,dirListRel[ifi],grpNameCoords)

      # Loop through 3d datasets and write to xdmf file
      for ids in range(0,nDsGrp):
        ii = ids*nTimes+itf # Get index for the correct time. Data is like f1,f2,f3,g1,g2,g3,h1,h2,h3,etc.
        nCellsRev = dsGrpShape[ii] # tuple object holding shape arrays
        nCells = nCellsRev[::-1] # hdf dimensions are reversed
        if (nCellsTopo!=nCells):
          print ('  ERROR: Number of cells in topology not same as dataset',nCellsTopo,nCells)
          return
        vtuGridAddData(fp,nDim3,nCells,dsGrpAttrName[ii],dirListRel[ifi],dsGrpNames[ii])
      vtuGridEnd(fp)
    
  vtuFooter(fp)
  fp.close()
  
# =============================================================================
# High level function to wrtie out all the xdmf content.
def writeRestartFile(nDim,wkdir,filePrefix,fileExt,attrPos,grpNameTime,grpNameStep,grpNameCoords):
  nDim3 = 3 # Structure of 2d xdmf file assumes a 3d setup with 2d hdf dataset inside.
    
  # Get list of all h5 files in this directory
  dirListAbs = glob.glob(wkdir+filePrefix+'*.'+fileExt)
  nFiles = len(dirListAbs)
  print( 'Number of files: ', nFiles)
  if (nFiles==0):
    return
    
  dirListRel = [''] * nFiles # Allocate array of strings
  for ifi in range(0,nFiles):
    dirListRel[ifi] = os.path.basename(dirListAbs[ifi]) # File name only, remove path

  fp = open(wkdir+filePrefix+'new.xdmf','w')
  vtuHeader(fp)

  # Loop on h5 files
  for ifi in range(0,nFiles):
    print( '  Processing file: ', dirListRel[ifi])
    fh5 = h5py.File(dirListAbs[ifi],'r')
    h5Ds = H5DatasetNames()
    fh5.visititems(h5Ds) # Add dataset items to a list object, h5Ds.dsNames

    # Get coords. Assume all x,y,z are present even for 2d, where one dimension will be 1
    nCellsTopoTmp = np.zeros(nDim3,dtype=np.int) # Always use 3 coords, x,y,z. In 2d assume one component will be 1 cell thick
    for id in range(0,len(nCellsTopoTmp)):
      tmp = fh5.get(grpNameCoords[id])
      nCellsTopoTmp[id] = tmp.shape[0]
    nCellsTopo = tuple(nCellsTopoTmp) # nx,ny,nz

    # Get dataset
    dsGrpNames = []
    dsGrpSize = []
    dsGrpShape = []
    dsGrpAttrName = []
    for ids in range(0,len(h5Ds.dsShape)): # Loop on all group objects to separate 1d, 2d, 3d objects.     
      if (len(h5Ds.dsShape[ids])==nDim): # Look for 2d or 3d datasets only
        dsGrpNames += ['/'+h5Ds.dsNames[ids]]
        dsGrpSize += [h5Ds.dsSize[ids]]
        if (nDim==2):
          # For 2d we need to know which coordinate is only 1 cell thick, this info is not present in
          # the dataset. Use the coords in this case.
          dsGrpShape += [nCellsTopo[::-1]] # reverse order for hdf file, nz,ny,nx
        else:
          # For 3d we can take dimensions directly from the dataset.
          dsGrpShape += [h5Ds.dsShape[ids]]
        dsGrpAttrName += [h5Ds.dsNames[ids].split('/')[attrPos]]
    
    # Get time values
    tval = fh5.get(grpNameTime)
    time = np.array(tval)[0]
    print( '  Time: ', time)
    
    # Get step values
    sval = fh5.get(grpNameStep)
    step = int(round(np.array(sval)[0]))
    print ('  Step: ', step )
        
    fh5.close()
    
    vtuGridStart(fp,time,nCellsTopo)
    vtuGridAddGeom(fp,nDim3,nCellsTopo,dirListRel[ifi],grpNameCoords)

    # Loop through 3d datasets and write to xdmf file
    nDsGrp = len(dsGrpNames)
    for ids in range(0,nDsGrp):
      nCellsRev = dsGrpShape[ids] # tuple object holding shape arrays
      nCells = nCellsRev[::-1] # hdf dimensions are reversed
      if (nCellsTopo!=nCells):
        print ('  ERROR: Number of cells in topology not same as dataset',nCellsTopo,nCells)
        return
      vtuGridAddData(fp,nDim3,nCells,dsGrpAttrName[ids],dirListRel[ifi],dsGrpNames[ids])
    vtuGridEnd(fp)
    
  vtuFooter(fp)
  fp.close()

# =============================================================================
# MAIN
# =============================================================================

# Default/initial values
codeVersion = "0.6.0"
fname_in = 'xdmfCreator.in'
sep_in = '#'

# -----------------------------------------------------------------------------
# Read inputs
print ("xdmfCreator, version " + codeVersion)

# Get command line arguments. If no arguments given use defaults
#if len(sys.argv) == 1:
#  nDim = nDim0
#elif len(sys.argv) == 2:
#  nDim = sys.argv[1]
#else:
#  print "ERROR: Too many command line arguments, only one allowed: xdmfCreator.py [nDim]"
#  sys.exit(0)

# Read the input file
print ("Reading input file: %s" % fname_in)
try:
  with open(fname_in) as file:
    lines_in = file.readlines()
except:
  print ("ERROR: Cannot open file: ",fname_in)
  sys.exit(0)

for line in lines_in:
   wordsOnly = line.split(sep_in,1)[0] # Remove comments after separator
   words = wordsOnly.split()
   if len(words) == 2:
      if words[0] == 'filePrefixRes2d':
         filePrefixRes2d = words[1].strip('\"')
      if words[0] == 'fileExtRes2d':
         fileExtRes2d = words[1].strip('\"')
      if words[0] == 'grpNameTimeRes2d':
         grpNameTimeRes2d = words[1].strip('\"')
      if words[0] == 'grpNameStepRes2d':
         grpNameStepRes2d = words[1].strip('\"')
      if words[0] == 'attrPosRes2d':
         tmp = np.fromstring(words[1], dtype=int, sep=" ")
         attrPosRes2d = tmp[0]

      if words[0] == 'filePrefixRes3d':
         filePrefixRes3d = words[1].strip('\"')
      if words[0] == 'fileExtRes3d':
         fileExtRes3d = words[1].strip('\"')
      if words[0] == 'grpNameTimeRes3d':
         grpNameTimeRes3d = words[1].strip('\"')
      if words[0] == 'grpNameStepRes3d':
         grpNameStepRes3d = words[1].strip('\"')
      if words[0] == 'attrPosRes3d':
         tmp = np.fromstring(words[1], dtype=int, sep=" ")
         attrPosRes3d = tmp[0]

      if words[0] == 'grpNameCoordsResx':
         grpNameCoordsResx = words[1].strip('\"')
      if words[0] == 'grpNameCoordsResy':
         grpNameCoordsResy = words[1].strip('\"')
      if words[0] == 'grpNameCoordsResz':
         grpNameCoordsResz = words[1].strip('\"')

      if words[0] == 'filePrefixPlt2d':
         filePrefixPlt2d = words[1].strip('\"')
      if words[0] == 'fileExtPlt2d':
         fileExtPlt2d = words[1].strip('\"')
      if words[0] == 'grpNameTimePlt2d':
         grpNameTimePlt2d = words[1].strip('\"')
      if words[0] == 'grpNameStepPlt2d':
         grpNameStepPlt2d = words[1].strip('\"')
      if words[0] == 'attrPosPlt2d':
         tmp = np.fromstring(words[1], dtype=int, sep=" ")
         attrPosPlt2d = tmp[0]

      if words[0] == 'filePrefixPlt3d':
         filePrefixPlt3d = words[1].strip('\"')
      if words[0] == 'fileExtPlt3d':
         fileExtPlt3d = words[1].strip('\"')
      if words[0] == 'grpNameTimePlt3d':
         grpNameTimePlt3d = words[1].strip('\"')
      if words[0] == 'grpNameStepPlt3d':
         grpNameStepPlt3d = words[1].strip('\"')
      if words[0] == 'attrPosPlt3d':
         tmp = np.fromstring(words[1], dtype=int, sep=" ")
         attrPosPlt3d = tmp[0]
 
      if words[0] == 'grpNameCoordsPltx':
         grpNameCoordsPltx = words[1].strip('\"')
      if words[0] == 'grpNameCoordsPlty':
         grpNameCoordsPlty = words[1].strip('\"')
      if words[0] == 'grpNameCoordsPltz':
         grpNameCoordsPltz = words[1].strip('\"')
 
      if words[0] == 'filePrefixRestart':
         filePrefixRestart = words[1].strip('\"')
      if words[0] == 'fileExtRestart':
         fileExtRestart = words[1].strip('\"')
      if words[0] == 'grpNameTimeRestart':
         grpNameTimeRestart = words[1].strip('\"')
      if words[0] == 'grpNameStepRestart':
         grpNameStepRestart = words[1].strip('\"')
      if words[0] == 'attrPosRestart':
         tmp = np.fromstring(words[1], dtype=int, sep=" ")
         attrPosRestart = tmp[0]
 
      if words[0] == 'grpNameCoordsRestx':
         grpNameCoordsRestx = words[1].strip('\"')
      if words[0] == 'grpNameCoordsResty':
         grpNameCoordsResty = words[1].strip('\"')
      if words[0] == 'grpNameCoordsRestz':
         grpNameCoordsRestz = words[1].strip('\"')

wkdir = os.getcwd() + os.path.sep

# ----- PLT FILES -------------------------------------------------------------
grpNameCoords = [grpNameCoordsPltx,grpNameCoordsPlty,grpNameCoordsPltz]

print ("")
print ("Processing 2d plt files")
writePltFile(2,wkdir,filePrefixPlt2d,fileExtPlt2d,attrPosPlt2d,grpNameTimePlt2d,grpNameStepPlt2d,grpNameCoords)

print ("")
print ("Processing 3d plt files")
writePltFile(3,wkdir,filePrefixPlt3d,fileExtPlt3d,attrPosPlt3d,grpNameTimePlt3d,grpNameStepPlt3d,grpNameCoords)

# ----- RESULTS FILES ----------------------------------------------------------------
grpNameCoords = [grpNameCoordsResx,grpNameCoordsResy,grpNameCoordsResz]

print ("")
print ("Processing 2d results files")
writeResFile(2,wkdir,filePrefixRes2d,fileExtRes2d,attrPosRes2d,grpNameTimeRes2d,grpNameStepRes2d,grpNameCoords)

print ("")
print ("Processing 3d results files")
writeResFile(3,wkdir,filePrefixRes3d,fileExtRes3d,attrPosRes3d,grpNameTimeRes3d,grpNameStepRes3d,grpNameCoords)

# ----- Restart FILES -------------------------------------------------------------
grpNameCoords = [grpNameCoordsRestx,grpNameCoordsResty,grpNameCoordsRestz]

print ("")
print ("Processing 3d restart files")
writeRestartFile(3,wkdir,filePrefixRestart,fileExtRestart,attrPosRestart,grpNameTimeRestart,grpNameStepRestart,grpNameCoords)

print ("Done.")
#sys.exit(0)
