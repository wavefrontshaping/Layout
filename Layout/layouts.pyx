from .core import Layout
from .logger import get_logger
import numpy as np
import matplotlib.pyplot as plt

logger = get_logger(__name__)

def checkLine(p1, p2, shape):
    """
    Uses the line foing from p1 to p2 and check each pixel of base_array
    to know if they are on the right side of the line.

    Returns boolean array, with True on right side and False on left side of the line.
    """
    idxs = np.indices(shape) # Create 3D array of indices

    p1 = np.array(p1).astype(float)#.astype(float)
    p2 = np.array(p2).astype(float)

    if abs(p2[0] - p1[0])> 1e-4:
        # Calculate max column idx for each row idx based on interpolated line between two points
        max_col_idx = (idxs[0] - p1[0]) / (p2[0] - p1[0]) * (p2[1] - p1[1]) +  p1[1] 
        # Get the direction of the line
        sign = np.sign(p2[0] - p1[0])
        return idxs[1] * sign <= max_col_idx * sign
    else: # vertical line
        sign = np.sign(p2[1] - p1[1])
        return idxs[0]*sign >= p1[0]*sign
    
def createPolygon(shape, vertices):
    """
    Creates np.array with dimensions defined by shape
    Fills polygon defined by vertices with ones, all other values zero
    """
    fill = np.ones(shape).astype(bool)  # Initialize boolean array defining shape fill

    # Create check array for each edge segment, combine into fill array
    for k in range(len(vertices)):
        fill = np.all([fill, checkLine(vertices[k-1], vertices[k], shape)], axis=0)


    return fill

class Squares(Layout):
    def __init__(self,radius,cellSize, resolution, center = None, gap = 0,verbose = True):
        
        Layout.__init__(self)  
        self._cellSize = (cellSize)
        self._radius = radius
        self._res = resolution
        self._surfaces = []
        self._grid = []
        self._gap = gap
        
        if center == None:
            self._center = [float(self._res[0]-1.)/2,float(self._res[1]-1.)/2]
        else:
            self._center = center
        
        if verbose:
            logger.info('Creation of hexagonal layout.')
#        print('Setting up the grid.')
#        self.setupGrid()
#        print('Creation of the hexagons.' )
        self._parts = []
        if verbose:
            logger.info('Creation of the hexagons.' )
        self.getSquareCell()
        if verbose:
            logger.info('Setting up the grid.')
        self.getSquareSegments()
        if verbose:
            logger.info(('-> Number of segments = %g' % self.nParts))
#        self.getFirstHexagons()
#        self.drawHexagons()
#        self.getHexagons()
#        self.checkOverlaps()
        # If no gap, we need to check that the segments do not overlap
        # We recompile the segments so that there is no overlap
        # Small disparities in the segment surfaces arise
        if self._gap == 0:
            if verbose:
                logger.info('Removing overlaps.')
            self.removeOverlaps()
            
        self.calculateSurfaces()
        if self._gap == 0:
            if verbose:
                logger.info(('-> Maximum relative variation of segment surfaces = %0.3f' % (float(max(self._surfaces)-min(self._surfaces))/float(max(self._surfaces)))))

        if verbose:
            logger.info('Sorting segments.')
        self.sortSegments()
        
              

        
    def getSquareCell(self):
        # need to draw two hexagons manually 
        # they will be used to generate the whole pattern
        cellSize = self._cellSize
        cellSize += np.mod(cellSize,2)
        self._resCell = [cellSize]*2
        r = self._cellSize/2
        angles = np.arange(-3*np.pi/4,-11*np.pi/4,-np.pi/2)
        self._firstHex = []
        self._oneSqua = np.argwhere(np.ones((self._cellSize,self._cellSize))).astype(int)

    

    def getSquareSegments(self):
        self.nParts = 0
        ny_max = int(np.floor(float(self._radius)/self._cellSize))
        rsq = []
        self._parts = []
        self._grid = []
        self._angles = []
        x_shift = int(np.floor(self._cellSize+self._gap))
       # x_shift -= np.mod(x_shift+1,2) 
        y_shift = int(np.floor(self._gap+self._cellSize))
        #y_shift += np.mod(y_shift+1,2) 
        for y in np.arange(-ny_max-1,ny_max+1):
            ind = np.mod(y,2)
            # one out of two line is shifted
            add_shift = ind*int(0.5*(self._cellSize+self._gap))
            if (float(self._radius)/y_shift)**2-y**2 > 0:
                nx_max = int(np.floor(np.sqrt((float(self._radius)/y_shift)**2-y**2)))
            else:
                nx_max = 2
#            print '<<', nx_max, ny_max, '>>'
            for x in np.arange(-nx_max-1,nx_max+1):
                pos = [int(self._center[0]+x*x_shift+add_shift),int(self._center[1]+y*y_shift)]
                _rsq = (pos[0]-self._center[0])**2 +  (pos[1]-self._center[1])**2
                if (_rsq <= self._radius**2):
                    self._parts.append((self._oneSqua+np.array(pos)-[self._resCell[0]/2,self._resCell[0]/2]).astype(int))
                    self._grid.append(pos)
                    self.nParts += 1
                    rsq.append(_rsq)
                    self._angles.append(np.arctan2(1.*(pos[0]-self._center[0]),1.*(pos[1]-self._center[1])))

        self._parts = np.array(self._parts)
        self._grid = np.array(self._grid)
        self._angles = np.array(self._angles)
        self._rparts = np.sqrt(rsq)
        
     
  
    

         



class Hexagons(Layout):
    def __init__(self,radius, cellSize, resolution, center = None, gap = 0,verbose = True, checkOverlaps = True):
        
        Layout.__init__(self)  
        self._cellSize = cellSize
        self._radius = radius
        self._res = resolution
        self._grid = []
        self._gap = gap
        
        if center == None:
            self._center = [float(self._res[0]-1.)/2,float(self._res[1]-1.)/2]
        else:
            self._center = center

        logger.info('Creation of the hexagonal layout')

        self._parts = []

        logger.debug('Creation of the hexagons and removal of the center column')
        self.getHexagonCell()
        self.removeOneCol()

        logger.debug('Setting up the grid')
        self.getHexagonSegments()

        logger.info('-> Number of segments = %g' % self.nParts)

        # We need to check that the segments do not overlap
        if checkOverlaps and self.checkOverlaps():
        
            # We recompile the segments so that there is no overlap
            # Small disparities in the segment surfaces arise
            logger.debug('Removing overlaps')
            self.removeOverlaps()
                
        self.calculateSurfaces()
        if self._gap == 0:
            logger.debug('-> Maximum relative variation of segment surfaces = %0.3f' % (float(max(self._surfaces)-min(self._surfaces))/float(max(self._surfaces))))
        
        logger.debug('Sorting segments accoring to distance from center (default)')
        self.sortSegments()
        
              
    def getSurface(self):
        return self._surfaces
        
    def getHexagonCell(self):
        # need to draw two hexagons manually 
        # they will be used to generate the whole pattern
        cellSize = int(2./np.sqrt(3)*self._cellSize)
        cellSize += np.mod(cellSize,2)
        self._resCell = [cellSize]*2
        r = 1./np.sqrt(3)*self._cellSize
        angles = np.arange(-np.pi/2,-2.5*np.pi,-np.pi/3)
        self._firstHex = []
        # First one
        vertices = [ [cellSize//2+r*np.cos(a),cellSize//2+r*np.sin(a)] for a in angles]
#        self._firstHex.append(createPolygon(self._resCell,vertices))
        self._oneHex = np.argwhere(createPolygon(self._resCell,vertices)).astype(int)

        return self._oneHex

    def removeOneCol(self):
        '''
        Removes the center column of self.oneHex
        '''
        center = int(np.max(self._oneHex[:,1]) / 2)
        singleHex = [[duo[0],duo[1]-(np.sign(duo[1]-center)+1)/2] for duo in self._oneHex if duo[1] != center]
        self._oneHex = np.array(singleHex)
        
        
        
 
        

    def getHexagonSegments(self):
        self.nParts = 0
        ny_max = int(np.floor(float(self._radius)/self._cellSize*2/np.sqrt(3)))
        rsq = []
        self._parts = []
        self._grid = []
        self._angles = []
        x_shift = int(np.floor(self._cellSize+self._gap))
        x_shift -= np.mod(x_shift+1,2) 
        y_shift = int(np.floor((self._gap+self._cellSize)*np.sqrt(3)/2))-1
        y_shift += np.mod(y_shift+1,2) 
        for y in np.arange(-ny_max-1,ny_max+1):
            ind = np.mod(y,2)
            # one out of two lines is shifted
            add_shift = ind*int(0.5*(self._cellSize+self._gap))
            if (float(self._radius)/y_shift)**2-y**2 > 0:
                nx_max = int(np.floor(np.sqrt((float(self._radius)/y_shift)**2-y**2)))
            else:
                nx_max = 2
            for x in np.arange(-nx_max-1,nx_max+1):
                # The -np.sign(x)*x was added to compensate the previously odd hexagon size
                # pos = [int(self._center[0]+x*x_shift+add_shift ),int(self._center[1]+y*y_shift - np.sign(x)*x)]
                pos = [int(self._center[0]+x*x_shift+add_shift ),int(self._center[1]+y*y_shift )]
                _rsq = (pos[0]-self._center[0])**2 +  (pos[1]-self._center[1])**2
                if (_rsq <= self._radius**2):
                    self._parts.append((self._oneHex+np.array(pos)-[self._resCell[0]/2,self._resCell[0]/2]).astype(int))
                    self._grid.append(pos)
                    self.nParts += 1
                    rsq.append(_rsq)
                    self._angles.append(np.arctan2(1.*(pos[0]-self._center[0]),1.*(pos[1]-self._center[1])))

        self._parts = self._parts
        self._grid = np.array(self._grid)
        self._angles = np.array(self._angles)
        self._rparts = np.sqrt(rsq)
        
     
  
    

    


class Diamonds(Layout):
    def __init__(self,radius, cellSize, resolution, center = None, gap = 0,verbose = True):
        
        Layout.__init__(self)  
        self._cellSize = cellSize
        self._radius = radius
        self._res = resolution
        self._surfaces = []
        self._grid = []
        self._gap = gap
        
        if center == None:
            self._center = [float(self._res[0]-1.)/2,float(self._res[1]-1.)/2]
        else:
            self._center = center
        
        if verbose:
            logger.info('Creation of hexagonal layout.')

        self._parts = []
        if verbose:
            logger.info('Creation of the hexagons.' )
        self.getHexagonCell()
        if verbose:
            logger.info('Setting up the grid.')
        self.getDiamondSegments()
        if verbose:
            logger.info(('-> Number of segments = %g' % self.nParts))

        # If no gap, we need to check that the segments do not overlap
        # We recompile the segments so that there is no overlap
        # Small disparities in the segment surfaces arise
        if self._gap == 0:
            if verbose:
                logger.info('Removing overlaps.')
            self.removeOverlaps()
            
        self.calculateSurfaces()
        if self._gap == 0:
            if verbose:
                logger.info(('-> Maximum relative variation of segment surfaces = %0.3f' % (float(max(self._surfaces)-min(self._surfaces))/float(max(self._surfaces)))))

        
        if verbose:
            logger.info('Sorting segments.')
        self.sortSegments()
        
              

        
    def getDiamondCell(self):
        # need to draw two diamond manually 
        # they will be used to generate the whole pattern
        cellSize = int(2./np.sqrt(3)*self._cellSize)
        cellSize += np.mod(cellSize,2)
        self._resCell = [cellSize]*2
        r = 1./np.sqrt(3)*self._cellSize
        angles = np.arange(-np.pi/2,-2.5*np.pi,-np.pi/2)
        self._firstHex = []
        # First one
        vertices = [ [cellSize//2+r*np.cos(a),cellSize//2+r*np.sin(a)] for a in angles]
        self._oneHex = np.argwhere(createPolygon(self._resCell,vertices)).astype(int)

        


        

    
        

    def getDiamondSegments(self):
        self.nParts = 0
        ny_max = int(np.floor(float(self._radius)/self._cellSize*2/np.sqrt(3)))
        rsq = []
        self._parts = []
        self._grid = []
        self._angles = []
#        x_shift = int(np.floor(self._cellSize+self._gap))-1
        x_shift = int((2*self._gap + self._cellSize))-2
#        x_shift -= np.mod(x_shift+1,2) 
#        y_shift = int(np.floor((self._gap+self._cellSize)*np.sqrt(3)/2))-1
        y_shift = int(np.floor(0.5*self._cellSize+self._gap))-1
        y_shift += np.mod(y_shift+1,2) 
        for y in np.arange(-ny_max-1,ny_max+1):
            ind = np.mod(y,2)
            # one out of two line is shifted
#            add_shift = ind*int(0.5*(self._cellSize+self._gap))
            add_shift = ind * 0.5 * (self._cellSize + 2*self._gap-2)#-ind*int(0.5*(np.sqrt(2)*self._gap + self._cellSize))
            if (float(self._radius)/y_shift)**2-y**2 > 0:
                nx_max = int(np.floor(np.sqrt((float(self._radius)/y_shift)**2-y**2)))
            else:
                nx_max = 2
#            print '<<', nx_max, ny_max, '>>'
            for x in np.arange(-nx_max-1,nx_max+1):
                pos = [int(self._center[0]+x*x_shift+add_shift),int(self._center[1]+y*y_shift)] # Warning, x and y and inverted, as they are in many parts of this code.
                _rsq = (pos[0]-self._center[0])**2 +  (pos[1]-self._center[1])**2
                if (_rsq <= self._radius**2):
                    self._parts.append((self._oneHex+np.array(pos)-[self._resCell[0]/2,self._resCell[0]/2]).astype(int))
                    self._grid.append(pos)
                    self.nParts += 1
                    rsq.append(_rsq)
                    self._angles.append(np.arctan2(1.*(pos[0]-self._center[0]),1.*(pos[1]-self._center[1])))

        self._parts = self._parts
        self._grid = np.array(self._grid)
        self._angles = np.array(self._angles)
        self._rparts = np.sqrt(rsq)
        
     
  
          
        

        
class Pie(Layout):
    def __init__(self,radius,nParts, resolution, center = None):
        
       
        Layout.__init__(self) 
        
        self._res = resolution        
        self.nParts = nParts
        

        if center == None:
            center = (self._res[0]//2,self._res[1]//2)
        
        X,Y = np.meshgrid(np.arange(self._res[1]),np.arange(self._res[0]))
        X = X-center[1]; Y = Y-center[0]

        R = np.sqrt(X**2+Y**2)
        Th = np.arctan2(X,Y)
        
        self._parts = []
        thetaBnd = np.linspace(-np.pi,np.pi,nParts+1)
        for i in range(nParts):
            self._parts.append( (R < radius)*(Th < thetaBnd[i+1])*(Th > thetaBnd[i]))
