import numpy as np

def simulate_Fourier_filtering(img, 
             center_fft, 
             window_fft_width, 
             fft_zoom = 4):
    '''
    Simulates a pinhole in the Fourier plane of the DMD.

    Parameters
    ----------
    img : np.array
        The input field.

    center_fft : tuple of int
        The center of the filter in the Fourier plane.

    window_fft_width : int
        The width of the filter in the Fourier plane.

    fft_zoom : float, optional
        Coefficient for increasing the resolution if the Fourier plane,
        increase of resolution in the Fourier plane is usually required 
        for accurate results.
        If fft_zoom  > 1, add zeros in the Fourier tranform.
        If fft_zoom = 1, no zero padding.
        See help of ``numpy.fft.fft2``

    Returns
    -------
    array : np.array
        An array representing the optical field in the plane conjugated to 
        the input one after filtering in the Fourier plane.
        It has the same resolution the input array.
    '''
    def find_nearest(array, value):
        idx = (np.abs(array - value)).argmin()
        return idx

    Nx,Ny = img.shape
    Fh = np.fft.fftshift(np.fft.fft2(img,s=[Nx*fft_zoom,Ny*fft_zoom]))
      
    ## Create the vector of the spatial frequencies
    SfreqX = np.fft.fftshift(np.fft.fftfreq(Nx*fft_zoom, d=1.))
    SfreqY = np.fft.fftshift(np.fft.fftfreq(Ny*fft_zoom, d=1.))
 
    ## Frequency grid
    [Sy,Sx] = np.meshgrid(SfreqY,SfreqX)

    ## Now we want to do a filtering of the spatial frequencies to keep only the
    ## -1 order

    ## First create the mask, we want to conseve the spatial frequencies around minus 
    # the carrier frequency ('-freq') with a window of size 'width'
    Mask1 = (Sx-center_fft[0])**2+(Sy-center_fft[1])**2 < (window_fft_width/2)**2
        
    ## We get the field in the Fourier plane after filtering
    Fh2 = Fh*Mask1

    ## We shift the spatial frequencies around zero to remove the effect of the angular tilt due to -1 order
    I1 = np.nonzero(Mask1)                  # gets the indices corresponding to Mask1
    Fh3 = np.zeros(Fh.shape,dtype = np.complex)
    # We copy the -1 order centered to the zero frequency    
    center_ind = [find_nearest(SfreqX, -window_fft_width/2),
                  find_nearest(SfreqY, -window_fft_width/2)]

    Fh3[I1[0]-np.min(I1[0]).astype(int)+center_ind[0], 
        I1[1]-np.min(I1[1]).astype(int)+center_ind[1]] = Fh2[I1[0],I1[1]]


    ## Get the field after the second length
    finalField = np.fft.ifft2(np.fft.ifftshift(Fh3))[0:Nx,0:Ny]

    return finalField