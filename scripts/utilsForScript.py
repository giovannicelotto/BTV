import numpy as np
import math

def getPdgMask(GenPart_pdgId):
    ''''Given array of pdg values return a mask true if the pdgid is in the pdgList, false otherwise'''
    pdgList = [
                                                        # Bottom mesons
                                                            521,    # b+
                                                            511,    # b0
                                                            531,    # b_s 
                                                            541,    # b_c
                                                        # Charmed mesons
                                                            411,   # D+
                                                            421,   # D0
                                                            431,   # D_s+
                                                            #441 eta_c no displacement
                                                        
                                                             
                                                            3122,   #lambda
                                                            3222,   #sigma+
                                                            3212,   #sigma0
                                                            3112,   #sigma-
                                                            3322,   #csi0
                                                            3312,   #csi-
                                                            3334,   #omega-

                                                            4122,   #lambda_c  
                                                            4222,   #sigma_c++
                                                            4212,   #sigma_c+
                                                            4112,   #sigma_c0
                                                            4232,   # csi_c+
                                                            4132,   # csi'c0
                                                            4332,   # omega^0_c
                                                            4412,   #xi_cc^+
                                                            4422,   #xi_cc^++
                                                            4432,   #omega^+_cc
                                                            4444,   #omega^++_ccc

                                                            5122,   #lambda_b0
                                                            5112,   #sigma_b-
                                                            5212,   #sigma_b0
                                                            5222,   #sigma_b+
                                                            5132,   #xi_b-
                                                            5232,   #xi_b0
                                                            5332,   #omega_b-
                                                            5142,   #xi_bc0
                                                            5242,   #xi_bc+
                                                            5342,   #omega_bc0
                                                            5512,   #xi_bb-
                                                            5532,   #omega_bb-
                                                            5542,   #Omega_bbc0
                                                            5554,   #Omega_bbb-
                                                            ]
    pdgMask = np.isin(abs(GenPart_pdgId), pdgList)
    return pdgMask

def distance_3d(point1, point2):
    """
    Calculates the distance between two points in three-dimensional space.
    
    Args:
        point1 (tuple): The coordinates of the first point (x, y, z).
        point2 (tuple): The coordinates of the second point (x, y, z).
        
    Returns:
        float: The distance between the two points.
    """
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance
