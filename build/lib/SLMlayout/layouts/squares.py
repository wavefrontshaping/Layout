from ..core import Layout, createPolygon, one_polygon_vertices, scale_coordinates
from ..core import get_logger
import numpy as np
import math
import matplotlib.pyplot as plt

logger = get_logger(__name__)

class Squares(Layout):
    def __init__(
        self,
        radius,
        cellSize,
        resolution,
        center = None,
        gap = 0,
        squareZone = False,
        verbose = True,
        ):
        
        Layout.__init__(self)  
        self._cellSize = (cellSize)
        self._radius = radius
        self._res = resolution
        self._surfaces = []
        self._grid = []
        self._gap = gap
        self._sqZone = squareZone
        
        if center == None:
            self._center = [float(self._res[0]-1.)/2,float(self._res[1]-1.)/2]
        else:
            self._center = center
        
        if verbose:
            logger.info('Creation of hexagonal layout.')
        self._parts = []
        if verbose:
            logger.info('Creation of the hexagons.' )
        self.getSquareCell()
        if verbose:
            logger.info('Setting up the grid.')
        self.getSquareSegments()
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
        ny_max = int(np.floor(float(self._radius)/(self._cellSize+ self._gap)))
        rsq = []
        self._parts = []
        self._grid = []
        self._angles = []
        x_shift = int(np.floor(self._cellSize+self._gap))
        y_shift = int(np.floor(self._gap+self._cellSize))
        for y in np.arange(-ny_max,ny_max+1):
            ind = np.mod(y,2)
            # one out of two line is shifted
            add_shift = ind*int(0.5*(self._cellSize+self._gap)) if not(self._sqZone) else 0
            if not(self._sqZone):
                if (float(self._radius)/y_shift)**2-y**2 > 0:
                    nx_max = int(np.floor(np.sqrt((float(self._radius)/y_shift)**2-y**2)))
                else:
                    nx_max = 2
            else:
                nx_max = int(np.floor(float(self._radius)/(self._cellSize+ self._gap)))
            for x in np.arange(-nx_max,nx_max+1):
                pos = [int(self._center[0]+x*x_shift+add_shift),int(self._center[1]+y*y_shift)]
                _rsq = (pos[0]-self._center[0])**2 +  (pos[1]-self._center[1])**2
                self._parts.append((self._oneSqua+np.array(pos)-[self._resCell[0]/2,self._resCell[0]/2]).astype(int))
                self._grid.append(pos)
                self.nParts += 1
                rsq.append(_rsq)
                self._angles.append(np.arctan2(1.*(pos[0]-self._center[0]),1.*(pos[1]-self._center[1])))

        self._parts = np.array(self._parts)
        self._grid = np.array(self._grid)
        self._angles = np.array(self._angles)
        self._rparts = np.sqrt(rsq)