from PySide2 import QtWidgets, QtCore, QtGui
from PySide2.QtWidgets import QWidget,QGraphicsObject, QGraphicsItem, QGraphicsView, QGraphicsScene, QVBoxLayout, QGraphicsSimpleTextItem, QLabel
from PySide2.QtGui import  QPen, QColor, QPainterPath, QPainter,QPainterPathStroker
from PySide2.QtPrintSupport import QPrinter
from PySide2.QtCore import Signal, Slot, QObject

from dolomites import ui,tiziano


class tiziano_widget (QWidget):


    def __init__(self, parent=None):
           QWidget.__init__(self)
           self.setWindowTitle("dolomites-tiziano")
           
           self.scene = ui.dolomites_scene()
           self.view  = ui.dolomites_view(self.scene)
           self.status_bar = QLabel()
           layout = QVBoxLayout();
           layout.addWidget(self.view)
           layout.addWidget(self.status_bar);
           self.setLayout(layout);
                    
           self.view.coord_changed.connect(self.func)
           
           return

       
    @Slot(str)   
    def func(self,msg):
        self.status_bar.setText(msg)
       
    def draw(self,_draw):
        """
        render a tiziano drawing
        """
        
        # check for input data
        if isinstance(_draw,tiziano.drawing) == False :
            print('Please provide a valid Tiziano drawing to plot!')
            return
        
        #x =[p.x for p in _draw.points]
        #y =[p.y for p in _draw.points]
        
        #for p in _draw.points:
            #self.scene.addItem(point_item(p.x,p.y))
        


        #fig, ax = plt.subplots()

        #ax.scatter(x,y, marker ="o")       
        for l in _draw.lines:
            _l = line_item(l.n1.x, l.n1.y, l.n2.x, l.n2.y)
            self.scene.addItem(_l)
           #if l.ph_tag > 0:
               #ax.text((l.n1.x+l.n2.x)/2, (l.n1.y+l.n2.y)/2, str(l.ph_tag),color='red')

        #for idx, l in enumerate(self.labels):
            #ax.text(l.x, l.y, l.ph_tag)
        
        self.view.fitInView(self.scene.sceneRect(),QtCore.Qt.KeepAspectRatio);
        return

            

class tiziano_item (QGraphicsObject):

    def __init__(self):
        QGraphicsObject.__init__(self)

        info = ''
        self.myColor = QtCore.Qt.blue
        self.pen = QPen(self.myColor, 1, QtCore.Qt.SolidLine,
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
   
   
class line_item(tiziano_item,tiziano.line):
    
    def __init__(self, x1,y1,x2,y2, parent=None):
        tiziano_item.__init__(self)

        self.pen.setWidthF(0.00001);
        self.setAcceptHoverEvents ( True );
        #self.setFlag(QGraphicsItem.ItemIgnoresTransformations, True );

        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        
        self.set_path();
        self.setZValue(1);
        #self.info = QString(tr("line length: %1")).arg(path.length());
        #physicalNumber = _physical;
        
        return
    
    def boundingRect(self):
        return QtCore.QRectF(self.path.boundingRect())

    def paint(self, painter, option, widget):
        painter.setPen(self.pen)
        painter.setBrush(self.brushColor)
        painter.drawPath(self.path)
        
        return
    
    def shape(self):
        stroke = QPainterPathStroker();
        stroke.setWidth(0.1);

        _path = QPainterPath(stroke.createStroke (  self.path )) ;
        _path.closeSubpath();
        return _path;

    def set_path(self):

   
        self.path.moveTo(self.x1, -self.y1);
        self.path.lineTo(self.x2, -self.y2);
        
        return 
    
    def paint(self, painter, option, widget):
        painter.setPen(self.pen)
        painter.setBrush(self.brushColor)
        painter.drawPath(self.path)
        
        return

class point_item(tiziano_item,tiziano.point):
    
    def __init__(self, x, y, ph=-1, parent=None):
        tiziano_item.__init__(self)
        tiziano.point.__init__(self,x,y,ph)

        self.pen.setWidthF(2);
        self.setAcceptHoverEvents ( True );

        self.dx = 5;
        self.dy = 5;

        

        # revert the y axis for the visualization in the graphicsview
        self.setX(x);
        self.setY(-y);

        #self.setFlags( ItemIsSelectable | ItemIgnoresTransformations);
        #self.setAcceptHoverEvents(true);
        self.setFlag(QGraphicsItem.ItemIsSelectable, True );
        self.setFlag(QGraphicsItem.ItemIgnoresTransformations, True );

        self.set_path();
        self.setZValue(1);
        #info = QString("x: %1; y: %2").arg(x).arg(y);

        
        
        
        return
    
    def boundingRect(self):
        return QtCore.QRectF(self.path.boundingRect())

    def paint(self, painter, option, widget):
        painter.setPen(self.pen)
        painter.setBrush(self.brushColor)
        painter.drawPath(self.path)
        
        return
    
    def shape(self):
        stroke = QPainterPathStroker();
        stroke.setWidth(3.);

        _path = QPainterPath(stroke.createStroke (  self.path )) ;
        _path.closeSubpath();
        return _path;

    def set_path(self):

   
        self.path =  QPainterPath();
        self.path.addRect(QtCore.QRectF(-self.dx/2,-self.dy/2,self.dx,self.dy));
        
        return 
    
    def paint(self, painter, option, widget):
        painter.setPen(self.pen)
        painter.setBrush(self.brushColor)
        painter.drawPath(self.path)
        
        return
