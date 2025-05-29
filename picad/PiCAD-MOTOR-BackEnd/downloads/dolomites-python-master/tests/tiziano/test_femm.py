from dolomites import tiziano
import numpy
import math

draw = tiziano.drawing()
draw.open_femm('SynRM1001_2poles.fem')

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


