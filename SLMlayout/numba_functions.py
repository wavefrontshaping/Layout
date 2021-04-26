from numba import jit, prange
import numpy as np

@jit(nopython = True)
def _getSingleProjection(parts, grid, array, method = 'sum'):
    out_vec = np.empty(len(parts), dtype = np.complex64)
    for idx,part in enumerate(parts):
        if method == 'complex' or  method == 'center':
            c_x, c_y = grid[idx][0], grid[idx][1]
            if method == 'complex':
                #angle_center = np.angle(array[coords_center[0],coords_center[1]])
                angle_center = np.angle(array[c_x, c_y])
                sum_abs_sub_array = 0
                for i_sub in range (len(part)):
                    sum_abs_sub_array += np.abs(array[part[i_sub,0],part[i_sub,1]])
                out_vec[idx] = sum_abs_sub_array * np.exp(1j * angle_center)
            elif method == 'center':
                out_vec[idx] = array[c_x, c_y]
        elif method == 'sum' or method == 'average':
            out_vec[idx] = 0
            for i_sub in range (len(part)):
                out_vec[idx] += array[part[i_sub,0],part[i_sub,1]]
        if method == 'average':
            out_vec[idx] = out_vec[idx]/len(part)
    return out_vec
 
@jit(nopython = True, parallel=True)
def _getStackProjections(parts, grid, array, method = 'sum'):
    n_array = array.shape[0]
    out = np.zeros((n_array, len(parts)), dtype = np.complex64)
    for i in prange(n_array):
        out[i] = _getSingleProjection(
                    parts, 
                    grid,
                    array[i], 
                    method
                )
    return out


@jit(nopython = True)
def _getBitPlaneFromVec(
        vec, 
        res,
        pos,
        lengths,
        parts,
        leePeriod, 
        angle,
        inversion
    ):
    
    dataSize = (res[0]*res[1]//8)
    tilt_y = np.cos(angle)
    tilt_x = np.sin(angle)
    vec_phi = np.angle(vec)
    vec_amp = np.abs(vec)



    
    
    partx = parts[:,0]
    party = parts[:,1]
    
    
    count_pos = 0

    img = np.zeros(dataSize, dtype = np.uint8)#(ct.c_ubyte*dataSize)()

    for i_vec in range(len(vec)):
        if vec_amp[i_vec] != 0:
            # encode the phase of the field in the spatial phase of the grating
            offset =  np.floor(vec_phi[i_vec]/(2.*np.pi)*leePeriod)

            for i_part in range(lengths[i_vec]):
                # first condition is the grating for encoding the phase
                # second condition encodes the amplitude by reducing removing
                # lines in the direction orthogonal to the first grating
                if ((tilt_y*party[count_pos+i_part] + offset + \
                     tilt_x*partx[count_pos+i_part]) % leePeriod) \
                         < (0.5*leePeriod) \
                and \
                   ((tilt_x*party[count_pos+i_part] - \
                     tilt_y*partx[count_pos+i_part]) % leePeriod) \
                        < (vec_amp[i_vec]*leePeriod):

                    img[pos[count_pos+i_part]] ^=(1<<(~party[count_pos+i_part]&7))

        count_pos += lengths[i_vec]

    if inversion:
        for i_part in range(dataSize):
            img[i_part] ^= 255
            
    return img


@jit(nopython = True)
def _getMaskFromBitPlane(bitPlane, res):
    
    img = np.empty((res[0], res[1]), dtype=np.uint8)

    for ix in range(res[0]):
        for iy in range(res[1]//8):
            for ibit in range(8):
                img[ix,iy*8+ibit] = (bitPlane[ix*res[1]//8+iy]>>7-ibit)&1

    return img

    