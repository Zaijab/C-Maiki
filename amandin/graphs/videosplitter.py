#!/usr/bin/env python

import cv2

def FrameCapture(path):
    vidObj = cv2.VideoCapture(path)
    count = 0
    success = 1
    while success:
        success, image = vidObj.read()
        cv2.imwrite("frame%d.jpg" % count, image)
        count += 1

if __name__ == '__main__':
    FrameCapture("/home/zjabbar/code/C-Maiki/data/coral/Ant_Atoll_Deep_Wall_FSM_16th_July_2016.mp4")
