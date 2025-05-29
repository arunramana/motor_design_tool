#!/usr/bin/python3
# requires svg.path, install it like this: pip3 install svg.path

# converts a list of path elements of a SVG file to simple line drawing commands
from svg.path import parse_path
from svg.path.path import Line
from xml.dom import minidom

# read the SVG file
doc = minidom.parse('christmas-tree.svg')
#doc = minidom.parse('christmas-star.svg')

path_strings = [path.getAttribute('d') for path
                in doc.getElementsByTagName('path')]
doc.unlink()

import numpy as np

# recover the points coordinates and store in two lists:
_alpha = []
_beta  = []
for path_string in path_strings:
    path = parse_path(path_string)
    for e in path:
        _alpha.append(e.point(0).real)
        _beta.append(e.point(0).imag)

alpha = np.array(_alpha)
beta  = np.array(_beta)

import math

t    = np.linspace(0,0.1,len(alpha)) # we take 2 periods
#a = np.zeros(len(alpha))
#b = np.zeros(len(alpha))
#c = np.zeros(len(alpha))

a    = alpha
b = -0.5*alpha + math.sqrt(3)/2*beta
c = -0.5*alpha - math.sqrt(3)/2*beta

from dolomites import fnc
from PySide6.QtWidgets import QApplication
import sys, numpy as np

app = QApplication([])
win = fnc.fnc_widget()


win.add_series(t,a,b,c)

win.show()

sys.exit(app.exec())
