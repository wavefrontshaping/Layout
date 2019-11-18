from .core import Layout
from .core import get_logger
import numpy as np
import math
from matplotlib.path import Path
import matplotlib.pyplot as plt

logger = get_logger(__name__)

def _scale_coordinates(generator, image_width, image_height, side_length, center, radius):
    scaled_width = int(image_width / side_length) + 2
    scaled_height = int(image_height / side_length) + 2
    scaled_center = [c / side_length for c in center]
    scaled_radius = radius / side_length

    for pos,coords in generator(scaled_width, scaled_height, scaled_center, scaled_radius):
        yield (pos[0]* side_length, pos[1]* side_length),[(x * side_length, y * side_length) for (x, y) in coords]

def createPolygon(shape, vertices):
    print('*'*10)
    print(shape)
    x, y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x,y)).T
    mask = Path(vertices).contains_points(points)
    return mask.reshape(shape)

def _one_polygon_vertices(x,y,radius,sides):
    angles = np.linspace(-np.pi/2, -2.5*np.pi,sides, endpoint=False)+np.pi/2
    # angles = np.arange(-np.pi/2,-2.5*np.pi,-np.pi/3)+np.pi/2
    # print(angles-angles2)
    return [ [x+radius*np.cos(a),y+radius*np.sin(a)] for a in angles]

def _one_hexagon_vertices(x,y,radius):
    return _one_polygon_vertices(x ,y, radius, sides = 6)

def _one_hexagon_path(x,y):
    h = math.sin(math.pi / 3)
    return [(x - .5,   (y-1) * h),
            (x + .5,   (y-1) * h),
            (x + 1,     y * h),
            (x + .5,   (y + 1) * h),
            (x - .5,   (y + 1) * h),
            (x - 1,     y * h)]

def _get_hexagon_cell(cellSize):
    '''
    Get one unit cell for 'equal' method that will be repeated identically across the modulated area.
    '''
    # need to draw two hexagons manually 
    # they will be used to generate the whole pattern
    size = int(2./np.sqrt(3)*cellSize)
    size += np.mod(size,2)+np.mod(cellSize,2)
    resCell = [size]*2
    r = 1./np.sqrt(3)*cellSize
    vertices = _one_hexagon_vertices(1.*size/2-.5, 1.*size/2-.5, radius = r)
    oneCell = np.argwhere(createPolygon(resCell,vertices)).astype(int)

    return resCell, oneCell


class Hexagons(Layout):
    '''
    Two methods: 
      method = 'equal' guarantees that, as long as there is no overlap, all macropixel will have the same exact shape.
      It creates a unique unit cell that is repeated across the modulated area. 
      However, even for gap = 0, there can be some unmodulated pixel.

      method = 'grid' guarantees that all the pixels in the modulated area are part of a macropixel. 
      However, macropixel may not all have the exact same shape.
      gap and checkOverlaps is disregarded when using the grid method. 
    '''
    def __init__(self,
                 radius, 
                 cellSize, 
                 resolution, 
                 method = 'equal',
                 center = None, 
                 gap = 0,
                 verbose = True, 
                 checkOverlaps = True):
        
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

        if method == 'equal':
            logger.info('Creation of the hexagonal layout')

            self._parts = []

            logger.debug('Creation of the unique cell')
            # Get one unit cell that will be repeated identically across the modulated area.
            self._resCell, self._oneCell = _get_hexagon_cell(self._cellSize)

            logger.debug('Setting up the grid')
            self._get_equal_hexagon_segments()

            # We need to check that the segments do not overlap
            if checkOverlaps and self.checkOverlaps():
            
                # We recompile the segments so that there is no overlap
                # Small disparities in the segment surfaces arise
                logger.debug('Removing overlaps')
                self.removeOverlaps()

        elif method == 'grid':

            self._generate_parts_grid()
        
        logger.info('-> Number of segments = %g' % self.nParts)

       
                
        self.calculateSurfaces()
        if self._gap == 0:
            logger.debug('-> Maximum relative variation of segment surfaces = %0.3f' % (float(max(self._surfaces)-min(self._surfaces))/float(max(self._surfaces))))
        
        logger.debug('Sorting segments accoring to distance from center (default)')
        self.sortSegments()

    def _generate_unit_hexagons_grid(self, width, height, center, radius):
        """
        Used only in the 'grid' method.
        Generate coordinates for a tiling of hexagons of unit size.
        """
        h = np.sin(np.pi / 3)

        # to ensure x = 0 exists
        x_max = width//2+3-(width//2 % 3)
        y_max = int(height / (2*h))
        
        for x in range(-x_max, x_max + 1, 3):
            for y in range(-y_max, y_max + 1):
                # Add the horizontal offset on every other row
                x_ = x if (y % 2 == 0) else x + 1.5   
                if (x_**2+(y*h)**2 < radius**2):
                    x_ += center[0]
                    y_ = y + center[1]/h
                    yield (x_,y_*h), _one_hexagon_path(x_,y_)
  


    def _generate_hexagons_mesh(self,*args, **kwargs):
        """
        Used only in the 'grid' method.
        Generates coordinates for a tiling of hexagons.
        """
        return _scale_coordinates(self._generate_unit_hexagons_grid, 
                                  *args, 
                                  image_width = self._res[0],
                                  image_height = self._res[1],
                                  side_length = self._cellSize/2*1./np.sin(np.pi/3),
                                  center = self._center,
                                  radius = self._radius,
                                  **kwargs)
            
    def _generate_parts_grid(self): 
        x, y = np.meshgrid(np.arange(self._res[1]), np.arange(self._res[0]))
        x, y = x.flatten(), y.flatten()
        
        points = np.vstack((x,y)).T
        
        self._grid = np.array([])
        self._rparts = np.array([])
        self._angles = np.array([])

        
        for ind,(pos,shape) in enumerate(self._generate_hexagons_mesh()):
            mask = Path(shape).contains_points(points)
            self._grid = np.append(self._grid,[pos])
            self._parts.append(np.array(np.where(mask.reshape(self._res))).transpose())
            #self._angles.append(np.arctan2(pos[0]-self._center[0],pos[1]-self._center[1]))
            self._angles = np.append(self._angles,[np.arctan2(pos[0]-self._center[0],pos[1]-self._center[1])])
            #self._rparts.append(np.sqrt((pos[0]-self._center[0])**2+(pos[1]-self._center[1])**2))
            self._rparts = np.append(self._rparts, np.sqrt((pos[0]-self._center[0])**2+(pos[1]-self._center[1])**2))

        self.nParts = ind+1 
              
    def getSurface(self):
        return self._surfaces
        
       

    def _get_equal_hexagon_segments(self):
        '''
        Used only in the 'equal' method. 
        Repeats the unit cell across the modulated area to create the different parts.
        '''
        self.nParts = 0
        # FIX ME
        ny_max = int(np.ceil(float(self._radius)/self._cellSize*2/np.sqrt(2)))
        rsq = []
        self._parts = []
        self._grid = []
        self._angles = []
        x_shift = int(np.floor(self._cellSize+self._gap))
        y_shift = int(np.floor((self._gap+self._cellSize)*np.sqrt(3)/2))+1
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
                pos = [int(self._center[1]+x*x_shift+add_shift ),int(self._center[0]+y*y_shift )]
                _rsq = (pos[0]-self._center[1])**2 +  (pos[1]-self._center[0])**2
                if (_rsq <= self._radius**2):
                    self._parts.append((self._oneCell+np.array(pos)-[self._resCell[0]/2,self._resCell[0]/2]).astype(int))
                    self._grid.append(pos)
                    self.nParts += 1
                    rsq.append(_rsq)
                    self._angles.append(np.arctan2(1.*(pos[0]-self._center[0]),1.*(pos[1]-self._center[1])))

        self._parts = self._parts
        self._grid = np.array(self._grid)
        self._angles = np.array(self._angles)
        self._rparts = np.sqrt(rsq)

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
        self._oneCell = np.argwhere(createPolygon(self._resCell,vertices)).astype(int)

        


        

    
        

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
                    self._parts.append((self._oneCell+np.array(pos)-[self._resCell[0]/2,self._resCell[0]/2]).astype(int))
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
