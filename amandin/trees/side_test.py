import numpy as np
import cv2
from matplotlib import pyplot as plt
import matplotlib as mpl
from numpy import linalg as LA
import json
import itertools
import os

#Settings
font = cv2.FONT_HERSHEY_SIMPLEX
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["axes.facecolor"] = "white"
mpl.rcParams["savefig.facecolor"] = "white"
mpl.rcParams["figure.figsize"]= (12,12)

# Import Text File
label_file = open('../data/tree/result2_2018_07_19.txt','r')
# Put text into a string
labels = label_file.read()
# Create dictionary of images to boxes {filename : list_of_boxes}
# Boxes of the form [LEFT,TOP,RIGHT,BOTTOM]
labelled_images = {}
for x in labels.split("Enter Image Path: ")[1:-1]:
    corners = []
    filename = x[:x.index(":")]
    for box in x.split("Bounding Box: ")[1:]:
        corners.append([int(corner[corner.index("=")+1:]) for corner in box.split(", ")])
    labelled_images[filename] = corners

for image, labels in labelled_images.items():
    print(len(labels))
