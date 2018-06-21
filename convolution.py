
import numpy as np
from scipy import signal
from scipy.ndimage import imread
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import morlet_test2 as morlet


sigma= [1, 3, 6]
theta = [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi*2/3, np.pi*3/4, np.pi*5/6]


b_img = imread("butterfly.jpg", True)
c_img= imread("noise.circle.jpg", True)


for i in range(0, 3, 1):
    for j in range(0, 6, 1):
        print(i)
        print (j)
        kernel = morlet.morletMatrix_real(-16, 16, sigma[i],theta[j])
        conv = signal.convolve2d(c_img, kernel, boundary='symm', mode='same')
        plt.figure()
        plt.imshow(np.absolute(conv), cmap = 'gray')
        plt.savefig("real_circle_sig{0}_theta{1}".format(sigma[i],j))




