from dolomites import tiziano
import numpy
import math

draw = tiziano.drawing()

p0 = draw.add_point(0,0)
p1 = draw.add_point(0,-1)
p2 = draw.add_point(0,1)
p3 = draw.add_point(-1,0)
p4 = draw.add_point(1,0)

l1 = draw.add_line(p1,p2,2000)
l2 = draw.add_line(p3,p4,2000)

#for i in range(100):
    #for j in range(i):
        #draw.add_point(i,j)
    
draw.add_arc(p1,p3,math.pi/2,10*math.pi/180);
p5 = draw.add_point(2,2)
#p6 = draw.add_point(1,1)
#draw.add_line(p5,p6)
#draw.plot()

draw.save('test.tiz')


draw = tiziano.drawing()
draw.load('example.tiz')

from dolomites import tiziano_ui
from PySide2.QtWidgets import QApplication
import sys,os

# to fix PySide2 issue on BigSur
os.environ['QT_MAC_WANTS_LAYER'] = '1'

app = QApplication([])
win = tiziano_ui.tiziano_widget()
win.draw(draw)
win.show()

sys.exit(app.exec_())


