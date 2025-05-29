from PySide6 import QtWidgets, QtCore, QtGui
from PySide6.QtWidgets import QWidget,QGraphicsObject, QGraphicsItem, QGraphicsView, QGraphicsScene, QVBoxLayout, QGraphicsSimpleTextItem, QGraphicsLineItem
from PySide6.QtGui import  QPen, QColor, QPainterPath, QPainter, QPainterPathStroker
from PySide6.QtPrintSupport import QPrinter
from PySide6.QtCharts import QChartView, QLineSeries


from PySide6.QtCore import  Signal, QPointF

import time,math
import numpy as np


class dolomites_view (QGraphicsView):

    coord_changed = Signal(str)

    def __init__(self,scene):
        super(dolomites_view, self).__init__(scene)
        self.setMouseTracking(True);
        self.CTRL  = False;
        self.setViewportUpdateMode(QGraphicsView.FullViewportUpdate)



    def wheelEvent(self,event):
        adj = (event.angleDelta().y())
        if (self.CTRL == True):
            if adj > 0:
                self.scale(1.05,1.05)
            if adj < 0:
                self.scale(0.95,0.95)

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

        QGraphicsView.keyPressEvent(self,e);

    def keyReleaseEvent(self,e):

            if (e.key()==QtCore.Qt.Key_Control):
                self.CTRL  = False;
#            if (e->key()==QtCore.Qt.Key_Shift)   SHIFT = false;

            QGraphicsView.keyReleaseEvent(self,e);

    def mouseMoveEvent(self, event):

        p = self.mapToScene(event.x(),event.y())
        #p = p*1000 #convert in mm the coordinates

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

# to draw a grid in the scene?
    #def drawBackground(painter, rect)

    #QPen pen;
    #painter.setPen(pen);
    #gridSize = 10;

    #left = int(rect.left()) - (int(rect.left()) % gridSize);
    #top = int(rect.top()) - (int(rect.top()) % gridSize);
    #QVector points;
    #for (qreal x = left; x < rect.right(); x += gridSize){
        #for (qreal y = top; y < rect.bottom(); y += gridSize){
            #points.append(QPointF(x,y));
        #}
    #}

    #painter.drawPoints(points.data(), points.size());
#}




class dolomites_item (QGraphicsObject):

    def __init__(self,color=QtCore.Qt.blue,parent=None):
        super().__init__(parent)

        info = ''
        self.myColor = color
        self.pen = QPen(self.myColor, 2, QtCore.Qt.SolidLine,
                        QtCore.Qt.RoundCap, QtCore.Qt.RoundJoin)

        self.normalColor   = self.pen.color().lighter(125);
        self.brushColor    = self.normalColor;
        self.hoverColor    = self.normalColor.lighter(150);
        self.selectedColor = self.normalColor.darker(200);


        self.path = QPainterPath();

        self.setFlag(QGraphicsItem.ItemIsSelectable, True );
        #self.setFlag(QGraphicsItem.ItemIsMovable, True);
        self.setAcceptHoverEvents(True)

    def hoverEnterEvent(self,event):

            if  QGraphicsObject.isSelected(self):
                self.brushColor = self.selectedColor;
            else:
                self.brushColor = self.hoverColor;
            QGraphicsItem.hoverEnterEvent(self,event);


    def hoverLeaveEvent(self,event):

            if QGraphicsItem.isSelected(self):
                self.brushColor = self.selectedColor;
            else:
                self.brushColor = self.normalColor;
            self.update();
            QGraphicsItem.hoverLeaveEvent(self,event);


    def itemChange(self, change, value):
       if change == QGraphicsObject.ItemSelectedHasChanged:
           if value == True:
               self.brushColor = self.selectedColor;
               self.pen.setColor(self.selectedColor);
               self.update();
           else:
               self.brushColor = self.normalColor;
               self.pen.setColor(self.normalColor);
               self.update();

       return QGraphicsObject.itemChange(self, change, value)

    def boundingRect(self):
        return QtCore.QRectF(self.path.boundingRect())

    def paint(self, painter, option, widget):
         painter.setPen(self.pen)
         painter.setBrush(self.brushColor)
         painter.drawPath(self.path)

    def shape(self):
         stroke = QPainterPathStroker();
         stroke.setWidth(3.);

         _path = QPainterPath(stroke.createStroke (  self.path )) ;
         _path.closeSubpath();
         return _path;
