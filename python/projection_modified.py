import numpy as np
from math import cos, sin, radians

def pointProject(X, yaw, pitch, dist, cameraTarget):
    p = radians(pitch)
    y = radians(yaw)
    target = np.array(cameraTarget)
    d  = target + np.array([cos(p)*sin(y), -cos(p)*cos(y), -sin(p)])*dist
    X_not_rot = X - d
    R = np.array([[cos(y),sin(y),0],[-sin(y)*sin(p), cos(y)*sin(p),-cos(p)]])
    projection = R @ X_not_rot
    # projection_2d = projection_3d[:2]
    return projection, R


def projectMatrix(yaw, pitch):
    p = radians(pitch)
    y = radians(yaw)
    R = np.array([[cos(y),sin(y),0],[-sin(y)*sin(p), cos(y)*sin(p),-cos(p)],[-sin(y)*cos(p),cos(y)*cos(p),sin(p)]])
    return R