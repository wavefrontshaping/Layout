# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:41:53 2016

core.py
============================
The core module of SLMlayout


@author: S.M. Popoff, M.W. Matthes
"""
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
import time
import numpy as np
cimport numpy as np
import matplotlib.pyplot as plt
import ctypes as ct
import itertools 
from matplotlib.path import Path
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython cimport array
import array
import logging


def get_logger(name):
        loglevel = logging.DEBUG
        logger = logging.getLogger(name)
        if not getattr(logger, 'handler_set', None):
            logger.setLevel(logging.INFO)
            logFormatter = logging.Formatter("%(asctime)s - %(name)-18.18s [%(levelname)-7.7s]  %(message)s") #[%(threadName)-12.12s] 
            consoleHandler = logging.StreamHandler()
            consoleHandler.setFormatter(logFormatter)
            logger.addHandler(consoleHandler)
            logger.setLevel(loglevel)
            logger.handler_set = True
        return logger

# the Layout class displays messages in a logging fashion. Output should be directed into logs and stdout in the main program.
# most messages are at DEBUG level, thus you can prevent their display by setting the threshold higher. 
logger = get_logger(__name__)

def logical_xor(str1, str2):
    return bool(str1) ^ bool(str2)


def createPolygon(shape, vertices):
    x, y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x,y)).T
    mask = Path(vertices).contains_points(points)
    return mask.reshape(shape)

def one_polygon_vertices(x,y,radius,sides):
    angles = np.linspace(-np.pi/2, -2.5*np.pi,sides, endpoint=False)+np.pi/2
    # angles = np.arange(-np.pi/2,-2.5*np.pi,-np.pi/3)+np.pi/2
    # print(angles-angles2)
    return [ [x+radius*np.cos(a),y+radius*np.sin(a)] for a in angles]

def scale_coordinates(generator, image_width, image_height, side_length, center, radius):
    scaled_width = int(image_width / side_length) + 2
    scaled_height = int(image_height / side_length) + 2
    scaled_center = [c / side_length for c in center]
    scaled_radius = radius / side_length

    for pos,coords in generator(scaled_width, scaled_height, scaled_center, scaled_radius):
        yield (pos[0]* side_length, pos[1]* side_length),[(x * side_length, y * side_length) for (x, y) in coords]


def fromFile(file_path):
    '''
    Return a Layout object created from a saved file.
    
    Parameters
    ----------
    file_path : string
        Path of the saved file corresponding to the Layout object to load.
        
    Returns
    -------
    layout : Layout
        The loaded Layout object.
    
    '''
    layout = Layout()
    data = np.load(file_path,  allow_pickle = True)
    layout.__dict__.update(dict(data))
    return layout

class Layout:
    '''
    The Layout super class.
    Contains the standard methods to generate mask from a vector
    and to handle the segments.
    '''
    def __init__(self):
        self._parts = []
        self.nParts = 0
        self._res = [0,0]
        self._pos_vec = []
        self._surfaces = None
        
    def saveToFile(self, file_path):
        '''
        Save the current state of the Layout to a .npz file (it uses numpy.savez)

        Parameters
        ----------
        file_path : string
            Path of the file to save the Layout object to.
        '''
        np.savez(file_path, **self.__dict__)
        
    def getMaskFromImage(self, complex_pattern, leePeriod = 1, angle = 0):
        '''
        Generate a binary amplitude mask from a complex image of the sime size as the full Layout.
        WARNING: This function is not optimized in Cython.

        Parameters
        ----------
        vec : np.array
            An image of the same size as the Layout corresponding to the target optical field. Can be complex valued.

        leePeriod : int, optional
            The period of the Lee Hologram grating. The default value is 1 (no grating).

        angle : float, optional
            The angle of the Lee Hologram grating. The default angle is 0 (vertical grating).

        Returns
        -------
        array : np.array
            The image that has the same resolution as the layout encoded in 8-bit.
        '''
        assert not np.max(np.abs(complex_pattern)) > 1.
  
        phase_shift = np.floor(np.angle(complex_pattern)/(2*np.pi)*leePeriod)
        X,Y = np.meshgrid(np.arange(self._res[1]),np.arange(self._res[0]))
        mask = (np.mod((np.cos(angle)*X+np.sin(angle)*Y+phase_shift),leePeriod) \
                < leePeriod*0.5) 
        mask *= (np.mod((-np.sin(angle)*X+np.cos(angle)*Y),leePeriod) \
                < leePeriod*np.abs(complex_pattern))
        mask = mask.astype(int)
        mask *= (2**8-1)
        return mask.astype(np.uint8)

    def getModulatedMask(self):
        '''
        Compute an binary mask corresponding to the mouldated area on the modulator.
        Return an array of the same resolution as the modulator with 1. where it corresponds to 
        a pixel that is modulated (is inside a macropixel) and 0. otherwise.
        '''

        modulated_mask = np.zeros([self._res[0],self._res[1]]) 
        for part in self._parts:
            modulated_mask[part[:,0],part[:,1]] = 1
        return modulated_mask

    def addPhaseRamp( 
            self,
            img,     
            pi2_value,
            period, 
            angle = 0.,
            mask = False,
            dtype = np.float):
        '''
        Add a phase ramp to a real pattern. 
        FOR PHASE MODULATOR ONLY.

        img : np.array
            The input image. Must have the same resolution as the modulator.

        pi_value : int
            Pixel value corrsponding to a 2pi modulation. 

        period : float
            FOR PHASE MODULATORS ONLY.
            Add a phase ramp to the image with a period of phase_ramp_period.
            If phase_ramp_period = None, no phase ramp is added (default).

        period : float, optional
            FOR PHASE MODULATORS ONLY.
            Set the angle of the phase ramp is phase_ramp_period is not None, ignored otherwise.
            Defaults to 0.

        mask: boolean, optional
            If mask = True, all the pixels that are not in modulated macropixels 
            are set to 0.
            Defaults to False.

        dtype : dtype, optional
            dtype of the returned array. Default dtype is np.float.

        '''
        
        X,Y = np.meshgrid(np.arange(self._res[0]),np.arange(self._res[1]))
        phase_ramp = (X*np.cos(angle) + Y*np.sin(angle))*pi2_value/period
        # wrap between 0 and pi_value
        modified = np.mod(phase_ramp + img, pi2_value+1)
        if mask:
            modulated_mask = self.getModulatedMask()
            modified *= modulated_mask
        return modified.astype(dtype)
   
    def getImageFromVec(
            self,
            vec, 
            overlap = False, 
            dtype = complex
            ):
        '''
        Creates a image from a list of values corresponding to the field one want to display on each segment of the layout.
        The returned mask can then be used to send images to a DMD using the ALP4lib module.


        Parameters
        ----------
        vec : list or np.array
            A list or vector containing the values of the field one wants to assign to each segment of the layout.

        overlap : bool, optional
            If overlap = True, the values will be added for pixels belonging to two or more segments of the layout (can happen with no gap between the segments)
        
        dtype : dtype, optional
            dtype of the returned array. Default dtype is complex.

        Returns
        -------
        array : np.array
            The image that has the same resolution as the layout.
        '''
        if not len(vec) == self.nParts:
            raise ValueError('Vector should have the same number of elements as the number of parts in the layout.')
        
        img = np.zeros([self._res[0],self._res[1]],dtype = dtype) 
        
        for part,x in zip(self._parts,vec):
            if overlap:
                img[part[:,0],part[:,1]] += x
            else:
                img[part[:,0],part[:,1]] = x

        return img



    def _calculatePosVec(self):
        self._pos_vec = []
        for part in self._parts:
            self._pos_vec.append((np.array(part[:,1])//8).astype('uint')+(part[:,0]*self._res[1]//8).astype('uint'))
    

    def getBitPlaneFromVec(self, vec, int leePeriod = 1, angle = 0, int inversion = False, dataFormat = 'Python'):
        '''
        Creates a bitplane from a list of values corresponding to the field one want to display on each segment of the layout.
        The returned bitplane can then be used to send images to a DMD using the ALP4lib module.
        Each byte corresponds to the binary value of 8 successive pixels.
        
        If leePeriod > 1, add a periodic modulation to create +1/-1 orders.
        vec can have a phase that will be discretized in the N=leePeriod phase values.
        The amplitude value is encoded by removing lines in the direction orthogonal to the Lee grating


        Parameters
        ----------
        vec : list or np.array
            A list or vector containing the complex values of the field one wants to assign to each segment of the layout.

        leePeriod : int, optional
            The period of the Lee Hologram grating. The default value is 1 (no grating).

        angle : float, optional
            The angle of the Lee Hologram grating. The default angle is 0 (vertical grating).

        inversion : bool, optional
            If inversion is set to True, all the bits are switched (0/1 becomes 1/0). The default value is False.
        
        dataFormat : string, optional
            if dataFormat == 'Python', return a numpy array, if dataFormat == 'C' returns a pointer to a C array (that can be sent using the ALP4lib SeqPut() function).
            By default, returns a numpy array (dataFormat = 'Python').
        

            
        Returns
        -------
        array : np.array or C pointer
            A bitplane that can be sent to the DMD using the ALP4lib SeqPut() or SeqPutEx() function.
        '''
 

        cdef int dataSize = <int>(self._res[0]*self._res[1]//8)
        cdef int offset
        cdef int cperiod = <unsigned int>leePeriod
        cdef float tilt_y = <float>np.cos(angle)
        cdef float tilt_x = <float>np.sin(angle)
        cdef double [:] mvec_phi = np.angle(np.array(vec)).astype(np.double)
        cdef double [:] mvec_amp = np.abs(np.array(vec)).astype(np.double)
        cdef int i_vec, i_part


        if not self._pos_vec:
            self._calculatePosVec()


        # Because of overlaps, elements may not have the same number of pixels,
        # it is then not easy to build a n-dimension array in C to fit the data
        # -> use 1D array and use a counter to increment the index
        cdef int [:] mpos = np.concatenate(self._pos_vec).astype(np.intc) 
        cdef int [:] mpartx = np.concatenate([p[:,0] for p in self._parts]).astype(np.intc) 
        cdef int [:] mparty = np.concatenate([p[:,1] for p in self._parts]).astype(np.intc)
        # store the lengths of each part for incrementing the counter 
        cdef int [:] clengths = np.array([len(p) for p in self._pos_vec]).astype(np.intc) 
        cdef unsigned int count_pos = 0

        
        img = (ct.c_ubyte*dataSize)()
 

            
        #for i,x in enumerate(vec):
        for i_vec from 0<= i_vec < <int>len(vec):

            if mvec_amp[i_vec] != 0:
                # encode the phase of the field in the spatial phase of the grating
                offset =  <int>(np.floor(mvec_phi[i_vec]/(2.*np.pi)*leePeriod))

                for i_part from 0<= i_part < clengths[i_vec] by 1:
                    # first condition is the grating for encoding the phase
                    # second condition encodes the amplitude by reducing removing lines in the direction orthogonal to the first grating
                    if ((tilt_y*mparty[count_pos+i_part] + offset + tilt_x*mpartx[count_pos+i_part]) % cperiod)  < (0.5*cperiod) and \
                       ((tilt_x*mparty[count_pos+i_part] - tilt_y*mpartx[count_pos+i_part]) % cperiod)  < (mvec_amp[i_vec]*cperiod) :

                            img[mpos[count_pos+i_part]] ^=(1<<(~mparty[count_pos+i_part]&7))

            count_pos += clengths[i_vec]

        if inversion:
            for i_part from 0<= i_part < dataSize by 1:
                img[i_part] ^= 255
                  
        # PyMem_Free(c_part)

        if dataFormat == 'C':
            return img
        elif dataFormat == 'Python':
            return np.asarray(img)
    
  

    def getMaskFromBitPlane(self,bitPlane):
        '''
        Generate the spatial mask as a numpy array from a bit plane.
        
        Requires the resolution to be a mutiple of 8.

        Parameters
        ----------
        bitPlane : ndarray
            A bitplane generated by getBitPlaneFromVec()
            
        Returns
        -------
        array : list
            An numpy array representing the binary amplitude mask as displayed by the DMD.
        '''
        assert(np.mod(self._res[1],8) == 0)

        cdef unsigned int [:,:] img = np.empty(self._res,dtype=np.uintc)
        cdef unsigned int ix, iy, ibit
        cdef unsigned char [:] mbp = bitPlane

        for ix from 0<= ix < self._res[0] by 1:
            for iy from 0<= iy < self._res[1]//8 by 1:
                for ibit from 0<= ibit < 8:
                    img[ix,iy*8+ibit] = (mbp[ix*self._res[1]//8+iy]>>7-ibit)&1

        return np.asarray(img)
  


    def getProjections(self, array, method = 'sum'):
        '''
        Returns the projection of an array on the segments of the layout.
        
        Parameters
        ----------
        array : ndarray
            It has to have the same dimension as the layout.
        
        method : string, optional
            If method = 'sum', returns the sum of the values of the array on each segment of the layout,
            if method = 'average', returns the average value of the array on each segment of the layout,
            if method = 'complex', returns a complex value whose amplitude is the sum of the absolute values
            in the segments and the phase is the phase at the center of the segments,
            if method = 'center', returns the value at the center of each segment.
            The default value is 'average'.
            
        Returns
        -------
        out_vec : list
            List of coefficient corresponding to the projection on each segment of the layout.
        '''
        if not method in ['sum','center','complex','average']:
            raise(ValueError,'Unvalid method.')
        

        out_vec = np.empty(self.nParts, dtype = complex)

        for idx,part in enumerate(self._parts):
            if method == 'complex' or  method == 'center':
                #coords_center = np.mean(part,axis = 0).astype(np.int)
                coords_center = tuple(zip(*[self._grid[idx]]))
                if method == 'complex':
                    #angle_center = np.angle(array[coords_center[0],coords_center[1]])
                    angle_center = np.angle(array[coords_center])
                    out_vec[idx] = np.sum(np.abs(array[part[:,0],part[:,1]])) * np.exp(1j * angle_center)
                elif method == 'center':
                    out_vec[idx] = array[coords_center]
            elif method == 'sum':
                out_vec[idx] = np.sum(array[part[:,0],part[:,1]])
            elif method == 'average':
                out_vec[idx] = np.mean(array[part[:,0],part[:,1]])
        return out_vec

    def getCenters(self,vec = None):
        '''
        Get the center position of the segments indexed by the input vector.
        '''  
        if vec:
            return self._grid[vec]
        else:
            return self._grid
    def sortSegments(self,order='dist2center', rearrange = True):
        '''
        Sort the segment with respect to the chosen criteria.
        Returns the index list of the sorted elements.
        
        Parameters
        ----------
        order : string, optional
            Sorting criteria, can take following values:
            - 'dist2center' :  Sorts the elements by there distance to the center of the figure (default).
            - 'angles' : Sorts the elements by there angular position.
            - 'random' : Randomly shuffles the elements.

        rearrange : bool, optional
            If rearrange is set to True, the order of the elements will be changed in the layout object,
            if not, the new order is returned via the index list but the layout object remains unchanged. 
            
        Returns
        -------
        idx : list
            List index of the sorted segments of the layout.
        '''
        if order == 'dist2center':
            idx = np.array(self._rparts).argsort().astype(int)
        elif order == 'angles':
            idx = np.array(self._angles).argsort().astype(int)
        elif order == 'random':
            idx = np.random.permutation(np.arange(self.nParts))  
        elif order == 'square':      
            list1 =  sorted(self._grid , key=lambda k: [k[0], k[1]])
            idx = np.zeros(self.nParts, dtype =int)
            for i in range(self.nParts):
                for j in range(self.nParts):
                    if np.all(list1[j] == self._grid[i]):
                        idx[j] = i
        else:
            raise ValueError('Invalid sorting order argument.')
        if rearrange:
            self._parts = [self._parts[i] for i in idx]
            self._rparts = self._rparts[idx]
            self._grid = self._grid[idx]
            self._angles = self._angles[idx]
        return idx

    def showLayout(self):
        '''
        Display the layout, shows each cell in a different color.
        '''
        pattern = self.getImageFromVec(np.arange(10,10+self.nParts), dtype = float)
        plt.figure()
        plt.imshow(pattern)
        plt.show()
        plt.colorbar()

    def calculateSurfaces(self):
        self._surfaces = []
        for part in self._parts:
            self._surfaces.append(part.shape[0])      


    def getSurface(self):
        if not self._surfaces:
            self.calculateSurfaces()
        return self._surfaces

    def checkOverlaps(self,display=False):
        pattern = np.zeros(self._res,dtype=float)
        for i in range(self.nParts):
            vec = [0.]*self.nParts
            vec[i] = 1
            pattern += self.getImageFromVec(vec, overlap=  True, dtype = float)

        if display:
            plt.figure();plt.title('Check overlaps' )
            plt.imshow(pattern,interpolation = 'None')

        if np.max(pattern) == 1:
            return False
        else: 
            return True

    def removeOverlaps(self):
        pattern = np.zeros(self._res,dtype=float)
        self._surfaces = []
        for i in range(self.nParts):
            vec = [0.]*self.nParts
            vec[i] = 1
            pattern += self.getImageFromVec(vec, overlap = True, dtype = float)
            # remove part in commom with one of the previous hexagon            

            ind = np.argwhere(pattern == 2)
            temp = self._parts[i].tolist()


            for pt in ind:
                temp.remove(pt.tolist())
            self._parts[i] = np.array(temp)           

            # Remove the duplicates from the sum pattern 
            pattern[pattern == 2] = 1


        


        
