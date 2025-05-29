from PySide6 import QtWidgets, QtCore, QtGui
from PySide6.QtWidgets import QWidget,QGraphicsObject, QGraphicsItem, QGraphicsView, QGraphicsScene, QVBoxLayout, QGraphicsSimpleTextItem
from PySide6.QtGui import  QPen, QColor, QPainterPath, QPainter, QPainterPathStroker
from PySide6.QtPrintSupport import QPrinter
import time,math

from dolomites import ui_6 as ui, koil
import numpy as np

class koil_widget (QtWidgets.QWidget):
    """
    A widget to visualize a winding based on dolomites/koil
    No calculation is performed in this class, only visualization purposes
    """


    def __init__(self, parent=None):
           QWidget.__init__(self)
           self.setWindowTitle("dolomites - koil ui")

           self.slots = list()

           self.scene = ui.dolomites_scene()



           self.view = ui.dolomites_view(self.scene)
           layout = QVBoxLayout();
           layout.addWidget(self.view);
           self.setLayout(layout);
           #self.view.save_PDF()

           self.win = koil.m_phase_winding()




    def draw_slots(self, Q, angles = None, R_slot = 20, R = -1):
        """
        Draw the Q slots.
        """

        if angles is None:
             angles = np.linspace(0,2 * math.pi, Q+1)

        if (R<0):
            R  = 2*R_slot*Q / math.pi;


        for i in range (Q):
            slot  = slot_item(i+1,R_slot,angles[i]);
            slot.setPos(R*math.cos(angles[i]),-R*math.sin(angles[i]))
            self.scene.addItem(slot);
            slot.set_label()
            self.slots.append(slot)


    def draw_winding(self,m_phase_winding, colors = None):

           self.win = m_phase_winding

           if self.win.star.t < 1 :
               print('winding not feasible!')
               return
           if self.win.star.t == None:
               print('empty winding!')

           # let draw the slots
           self.draw_slots(self.win.star.Q)

           for w_ph in self.win.windings:
               if colors == None: # let select a valid color for the phase
                    c = list(np.random.choice(range(256), size=3))
                    color = QColor(c[0],c[1],c[2])
               # else:
                   # c=colors[]

               for coil in w_ph.coils:
                    begin = coil.begin
                    end   = coil.end
                    nc    = coil.nc
                    c = coil_item(color, self.slots[begin-1], self.slots[end-1],nc)
                    c.set_incidence(nc)
                    self.scene.addItem(c);

           self.view.fitInView(self.scene.sceneRect(),QtCore.Qt.KeepAspectRatio)


class slot_item (ui.dolomites_item):
    """
    A graphical representation of a slot.
    There are 2 layers, L1 the inner and L2 the outer
    id:    the index of the slot 1, 2, ..., Q
    R:     the radius at which the slot is placed
    angle: the angle where the slot is located.
    """


    def __init__(self,id,R,angle):
        super().__init__()
        # self.path.addEllipse(-R,-R,2*R,2*R);
        self.id = id
        self.angle = angle

        self.R = R

        self.L1 = slot_layer_item(self,0,R,id)
        self.L2 = slot_layer_item(self,0,R,id)

        self.L1.setPos(-R,0)
        self.L2.setPos(R,0)

        self.setRotation(-angle*180/math.pi)

    def set_label(self):
        # this must be invoked after the slot has been added to the scene
        label = QGraphicsSimpleTextItem(str(self.id),None);
        self.scene().addItem(label)
        center = label.boundingRect().center() # the center of the label

        # adjust the label position
        P = self.L1.scenePos()
        k = (P.manhattanLength()-self.R*2.)/P.manhattanLength()
        label.setPos(k*P.x()-center.x(),k*P.y()-center.y())



class slot_layer_item (ui.dolomites_item):
    """
    A graphical representation of one layer of the slot.
    The incidence describes if the coil is inword or outword in the layer
    """

    def __init__(self,parent,incidence,R,id):
        super().__init__(QtCore.Qt.black,parent)
        self.id = id
        self.R = R
        self.incidence = incidence
        self.set_path()

        self.setAcceptHoverEvents(False)
        self.normalColor   = self.pen.color().lighter(125);
        self.brushColor    = QtCore.Qt.white;
        # self.hoverColor    = self.normalColor.lighter(150);
        # self.selectedColor = self.normalColor.darker(200);


        self.setFlag(QGraphicsItem.ItemIsSelectable, False );
        self.setFlag(QGraphicsItem.ItemIsMovable, False);

    def set_path(self):

        self.path = QPainterPath();
        self.path.setFillRule(QtCore.Qt.WindingFill);
        self.path.addEllipse(-self.R,-self.R,2*self.R,2*self.R);

        if (self.incidence > 0):
            self.path.addEllipse(-self.R/4,-self.R/4,2*self.R/4,2*self.R/4);

        if (self.incidence < 0):
            self.path.moveTo(-self.R/4,-self.R/4);
            self.path.lineTo( self.R/4, self.R/4);
            self.path.moveTo(-self.R/4, self.R/4);
            self.path.lineTo( self.R/4,-self.R/4);

    def set_incidence(self,inc):

        self.incidence = inc
        self.set_path()


class coil_item (ui.dolomites_item):
        """
        A representation of a single coil of the winding
        """

        def __repr__(self):
            return "coil (%s, %s, %s)" % (self.begin.id, self.end.id, self.turns)

        def __init__(self, color, begin, end, turns = 1, x = -1, y = -1):
            super().__init__(color,None)

            # some variables to track mouse movement
            self.mouse_press = False;
            self.xold  = 0
            self.yold  = 0
            self.xmouse = 0
            self.ymouse = 0

            # the slots where the coil is inserted
            self.begin = begin
            self.end = end
            self.setZValue(-1); # The coil is behind the slots

            self.turns = turns;

            P1 = begin.L1.scenePos(); # recover the outer layer of slot 1
            xs1 = P1.x();
            ys1 = P1.y();
            P2 = end.L2.scenePos(); # recover the inner layer of slot 2
            xs2 = P2.x();
            ys2 = P2.y();

            self.width = math.sqrt((xs2-xs1)*(xs2-xs1)+(ys2-ys1)*(ys2-ys1));
            if (x == -1):
                 self.x = -self.width / 2; # set the default values
            if (y == -1):
                 self.y =  self.width / 2;

            self.cp = control_point_item(self,x,y);
            self.cp.setPos(self.x,self.y);

            self.cp.moved.connect(self.cp_moved);

            self.center = QtCore.QPointF( (xs1+xs2)/2 , (ys2+ys1)/2 )

            self.set_path();

            # move and rotate the coil on the two slots and rotate the label
            self.moveBy(self.center.x(),self.center.y());
            self.angle = math.atan2(ys2-ys1,xs2-xs1)*180/math.pi;
            self.setRotation(self.angle);

            self.itemChange(QGraphicsItem.ItemSelectedHasChanged, False);

            # set the default incidence for the coil
            # incidence = _pos;
            # info = tr("Coil with %1 turns").arg(turns);

        @QtCore.Slot(float,float)
        def cp_moved(self,x,y):
            # print(x,y)
            self.prepareGeometryChange();
            self.x = x;
            self.y = y;
            self.set_path();
            self.update();

            return

        def set_path(self):
            stroke = QPainterPathStroker();
            _path = QPainterPath();

            stroke.setWidth(8);
            _path.moveTo(-self.width/2, 0);
            _path.cubicTo(-self.width/2,0,self.x,self.y,0,self.y);
            _path.cubicTo(0,self.y,-self.x,self.y,self.width/2,0);

            self.path = QPainterPath(stroke.createStroke (  _path )) ;
            self.path.closeSubpath();


        def set_incidence(self,inc):
            self.begin.L1.set_incidence(inc)
            self.end.L2.set_incidence(-inc)


        def itemChange(self, change, value):

            if (change == QGraphicsItem.ItemSelectedHasChanged):
                if (value == True):
                    self.cp.setVisible(True);
                else:
                    self.cp.setVisible(False);
            self.update();

            return super().itemChange(change, value);


        def mousePressEvent (self, event ):
            if(event.button()==QtCore.Qt.LeftButton):
                self.mouse_press = True;
                self.xold  = self.cp.pos().x();
                self.yold  = self.cp.pos().y();
                self.xmouse = event.pos().x();
                self.ymouse = event.pos().y();
            return ui.dolomites_item.mousePressEvent(self,event);

        def mouseReleaseEvent(self,event):
            if(event.button()==QtCore.Qt.LeftButton):
                self.mouse_press = False;
            return ui.dolomites_item.mouseReleaseEvent(self,event);

        def mouseMoveEvent(self,event):
            x = event.pos().x();
            y = event.pos().y();
            if (self.mouse_press):
                if (self.xmouse <=0):
                     self.cp.setX(self.xold+x-self.xmouse);
                else:
                     self.cp.setX(self.xold-x+self.xmouse);
            self.cp.setY(self.yold+y-self.ymouse);

            return ui.dolomites_item.mouseMoveEvent(self,event);

class control_point_item (ui.dolomites_item):
        """
        An helper class to manage coil_item transformations
        """

        moved = QtCore.Signal(float, float)

        def __init__(self,parent, x, y):
            super().__init__(QtCore.Qt.red,parent)

            self.setFlags(QGraphicsItem.ItemIsMovable|QGraphicsItem.ItemSendsGeometryChanges);

            self.dx = 10;
            self.dy = 10;

            self.path.addEllipse(-self.dx/2,-self.dy/2,self.dx,self.dy);



        def itemChange(self, change, value):
             if change == QGraphicsItem.ItemPositionHasChanged:
                self.moved.emit(self.pos().x(),self.pos().y())
             return ui.dolomites_item.itemChange(self, change, value);
