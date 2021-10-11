# Credit to the GPCam package for the following functions.
import numpy as np

def kernel_l2_single_task(x1, x2, hyperparameters, obj):
    ################################################################
    ###standard anisotropic kernel in an input space with l2########
    ###########################for single task######################
    """
    x1: 2d numpy array of points
    x2: 2d numpy array of points
    obj: object containing kernel definition

    Return:
    -------
    Kernel Matrix
    """
    hps = hyperparameters
    distance_matrix = np.zeros((len(x1),len(x2)))
    #you can always get some help:
    #help(obj.squared_exponential_kernel)
    for i in range(len(x1[0])-1):
        distance_matrix += abs(np.subtract.outer(x1[:,i],x2[:,i])/hps[1+i])**2
    distance_matrix = np.sqrt(distance_matrix)
    #return   hps[0] * obj.squared_exponential_kernel(distance_matrix,1)
    #return   hps[0] * obj.exponential_kernel(distance_matrix,1) + obj.periodic_kernel(abs(x1[:,-1]-x2[:,-1]), hps[-1],hps[-2])
    return   hps[0] *  obj.matern_kernel_diff1(distance_matrix,1)
    #return   hps[0] * obj.matern_kernel_diff2(distance_matrix,1)

