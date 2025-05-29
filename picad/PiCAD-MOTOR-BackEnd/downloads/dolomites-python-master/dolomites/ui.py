from PySide2 import QtWidgets, QtCore, QtGui
from PySide2.QtWidgets import QWidget,QGraphicsObject, QGraphicsItem, QGraphicsView, QGraphicsScene, QVBoxLayout, QGraphicsSimpleTextItem
from PySide2.QtGui import  QPen, QColor, QPainterPath, QPainter
from PySide2.QtPrintSupport import QPrinter

from PySide2.QtCore import  Signal

import time,math

class dolomites_view (QGraphicsView):

    coord_changed = Signal(str)
    
    def __init__(self,scene):
        super(dolomites_view, self).__init__(scene)
        self.setMouseTracking(True);
        self.CTRL  = False;



    def wheelEvent(self,event):
        adj = (event.delta()/120) * 0.1
        if (self.CTRL == True):
            self.scale(1+adj,1+adj)
        else:
            QGraphicsView.wheelEvent(self,event);

    def keyPressEvent(self,e):

        scene = QGraphicsView.scene(self);

        if (e.key()==QtCore.Qt.Key_Control):
                self.CTRL  = True;
#            if (e.key()==QtCore.Qt.Key_Shift):
#                SHIFT == True;
        if (e.key()==QtCore.Qt.Key_Escape):
                scene.clearSelection();

    def keyReleaseEvent(self,e):

            if (e.key()==QtCore.Qt.Key_Control):
                self.CTRL  = False;
#            if (e->key()==QtCore.Qt.Key_Shift)   SHIFT = false;


    def mouseMoveEvent(self, event):
        
        p = self.mapToScene(event.x(),event.y())
        p = p*1000 #convert in mm the coordinates

        self.coord_changed.emit(str('({:6.3f}, {:6.3f}) (|{:6.3f}|,<{:5.2f})[mm,deg]'.format(p.x(),-p.y(),math.sqrt(p.x()**2+p.y()**2),math.atan2(-p.y(),p.x())*180/math.pi)))
        
        

        QGraphicsView.mouseMoveEvent (self, event );

    def mousePressEvent (self,event):

            if (event.button()==QtCore.Qt.LeftButton):
                self.setDragMode(QGraphicsView.ScrollHandDrag);
            p = self.mapToScene(event.x(),event.y())
            print(math.sqrt(p.x()**2+p.y()**2))
            QGraphicsView.mousePressEvent (self,event );

    def mouseReleaseEvent (self,event):

            if (event.button()==QtCore.Qt.LeftButton):
                self.setDragMode(QGraphicsView.NoDrag);
            QGraphicsView.mouseReleaseEvent (self, event );

    def save_PDF(self,file_name=''):
            homeLocation = "/Users/lalberti/luigi.pdf";
            if (file_name == ''):
                return False;
            printer = QPrinter(QPrinter.HighResolution);
            printer.setPaperSize(QPrinter.A4);
            printer.setOutputFormat(QPrinter.PdfFormat);
            printer.setOutputFileName(file_name);
#            printer.setFullPage(True)

            painter = QPainter(printer);
            self.scene().render(painter);
            painter.end()


class dolomites_scene (QGraphicsScene):

    def __init__(self):
        QGraphicsScene.__init__(self)


