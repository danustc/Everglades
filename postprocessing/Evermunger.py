'''
This is the data munging class for Everglades.
'''
import numpy as np
import matplotlib.pyplot as plt
import glob
import tifffile as tf

global_datapath = '/home/sillycat/Programming/Python/data_test/Everglades/Nov27_2018/'


sq2 = np.sqrt(0.5)
deskew_mat = np.array([
    [0., -1, 0.],
    [sq2, 0., sq2],
    [-sq2, 0., sq2]
    ]) # the deskew matrix between the lab frame and the image frame

def deskew(stack, scan_range = 6.0, pxl = 0.103):
    '''
    img: the raw image acquired
    scan_range: the range of Z-scanning, unit micron.
    pxl: pixel size, unit micron.
    have a zero-padded stack.
    '''
    NZ, NY, NX = stack.shape
    yy = np.arange(NY)*pxl
    n_pad = int(np.ceil(scan_range/(2*np.sqrt(2.)*pxl)))# half number of padded points
    padded_stack = np.zeros((NZ, NY + 2*n_pad, NX))
    DXB = np.linspace(0., scan_range, NZ) # The scan range
    z_range = DXB/np.sqrt(2.)
    y_range = -z_range
    for nn in range(NZ):
        y_shift = int(z_range[NZ-nn-1]/pxl)
        print('y_shift:', y_shift)
        padded_stack[nn, y_shift:(y_shift + NY), :] = stack[nn]

    # next, shift the y back
    return padded_stack



def main():
    '''
    load and test the data
    '''
    SC_list = glob.glob(global_datapath + '*SC*.npy')
    SC_psf = np.load(SC_list[0])
    padded_stack = deskew(SC_psf).astype('uint16')

    tf.imsave('padded.tif', padded_stack)
    tf.imsave('original.tif', SC_psf.astype('uint16'))

if __name__ == '__main__':
    main()



