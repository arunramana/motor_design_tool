from PySide6 import QtWidgets, QtCore, QtGui
from PySide6.QtWidgets import \
     QWidget,QGraphicsObject, QGraphicsItem,QGraphicsLineItem, QGraphicsView, \
     QGraphicsScene, QVBoxLayout, QGraphicsSimpleTextItem, QLabel, QGraphicsPathItem
from PySide6.QtGui import  QPen, QColor, QPainterPath, QPainter,QPainterPathStroker
from PySide6.QtPrintSupport import QPrinter
from PySide6.QtCore import Signal, Slot, QObject, QPointF
from PySide6.QtCharts import QChart, QChartView, QLineSeries, QScatterSeries

from dolomites import ui_6 as ui

import numpy as np
from math import e, pi
from cmath import *



class fnc_widget (QWidget):


    def __init__(self, parent=None):
           QWidget.__init__(self)
           self.setWindowTitle("dolomites-fnc")

           # space vector scene and view
           self.scene = ui.dolomites_scene()
           self.view  = ui.dolomites_view(self.scene)

           # status bar
           self.status_bar_sv = QLabel()
           self.status_bar_time = QLabel()


            # chart for time quantities

           chart = QChart()
           #chart_view = QChartView(chart, main_window)
           self.chart_view = fnc_chart_view(chart)
           self.chart_view.time_changed.connect(self.update_space_vector)

           chart.setTheme(QChart.ChartThemeDark)
           #chart.legend().hide()

           # setup layout
           self.view.setMinimumSize(600,350)
           self.chart_view.setMinimumSize(600,350)

           layout = QVBoxLayout();
           layout.addWidget(self.view)
           layout.addWidget(self.status_bar_sv);
           layout.addWidget(self.chart_view)
           layout.addWidget(self.status_bar_time);
           self.setLayout(layout);
           return

    @Slot(float,float,float,float)
    def add_series(self,t,a,b,c):
       self.chart_view.add_series(t,a,b,c)
       self.chart_view.chart().createDefaultAxes()


       # generate the sv trajectory
       sv = 2/3*(a + b*e**(2j*pi/3) + c*e**(4j*pi/3))

       trajectory =  QPainterPath();
       trajectory.moveTo(sv[0].real,-sv[0].imag)

       for p in sv[1:]:
          trajectory.lineTo(p.real,-p.imag)

       trajectory_path = QGraphicsPathItem(trajectory)
       trajectory_path.setPen(QPen(QColor("blue"), 2))

       self.scene.addItem(trajectory_path)

       r = point_item(sv[0].real,sv[0].imag)
       p = point_item(0,0,QtCore.Qt.darkRed) # origin
       self.a = arrow_item(r,p)


       self.scene.addItem(p)
       self.scene.addItem(r)
       self.scene.addItem(self.a)

       # initialize the time line at zero
       self.chart_view.update_time_line(0)



    def keyPressEvent(self,e):

       if e.key() == QtCore.Qt.Key_Space:
           self.chart_view.toggle_timer()

       QWidget.keyPressEvent(self,e);




    @Slot()
    def update_space_vector(self,t,a,b,c):


        sv = 2/3*(a + b*e**(2j*pi/3) + c*e**(4j*pi/3))

        self.a.set_tip(sv.real,sv.imag)
        msg = '<pre>({:-8.3f} + j {:-8.3f}) - ({:-8.3f}; {:-8.3f} rad [{:-8.3f} deg])</pre>'.format(sv.real, sv.imag, abs(sv), phase(sv),np.degrees(phase(sv)))
        self.status_bar_sv.setText(msg)

        msg = '(t:{:-8.3f}; a: {:-8.3f}; b: {:-8.3f}; c: {:-8.3f})'.format(t,a,b,c)
        self.status_bar_time.setText('<pre>' + msg + '</pre>')


        return



    @Slot(str)
    def print_msg(self,msg):
        self.status_bar_time.setText(msg)




class fnc_chart_view (QChartView):

    time_changed = Signal(float,float,float,float)


    def __init__(self,chart):
        super(fnc_chart_view, self).__init__(chart)
        self.setMouseTracking(True);
        #self.CTRL  = False;
        #self.setViewportUpdateMode(QGraphicsView.FullViewportUpdate)

        self.time_line = QGraphicsLineItem(0,0,0,0)
        self.time_line.setPen(QPen(QColor("white"), 2))
        #self.time_line.setFlag(QGraphicsItem.ItemIsSelectable, True );
        #self.time_line.setFlag(QGraphicsItem.ItemIsMovable, True );
        self.scene().addItem(self.time_line)

        self.max_y = 0
        self.max_t = 0
        self.min_y = 0
        self.min_t = 0

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.advance_time)
        self.timer.setInterval(50);

        self.time = 0
        self.dt = 0.5e-3; # the time step for the animation



    def toggle_timer(self):

        if self.timer.isActive():
            self.timer.stop()
        else:
             self.timer.start()

        return



    def add_series(self,t,y1,y2,y3):

           # delete previous drata
           self.chart().removeAllSeries()

           line1 = QLineSeries()
           line2 = QLineSeries()
           line3 = QLineSeries()

           line1.setName("a");
           line2.setName("b");
           line3.setName("c");
           #line.setPen(QPen(QColor("black"), 2))
           for _t, _y1, _y2, _y3 in np.nditer([t,y1,y2,y3]):
                line1.append(_t, _y1)
                line2.append(_t, _y2)
                line3.append(_t, _y3)
           #line.hovered.connect(on_line_hover)
           self.chart().addSeries(line1)
           self.chart().addSeries(line2)
           self.chart().addSeries(line3)

           self.max_y = np.max([np.max(y1),np.max(y2),np.max(y3)])
           self.min_y = np.min([np.min(y1),np.min(y2),np.min(y3)])
           self.max_t = np.max(t)
           self.min_t = np.min(t)

           # store a copy of the np arrays, QLineSeries is not suitable for interpolation
           self.t  = t
           self.y1 = y1
           self.y2 = y2
           self.y3 = y3


    def update_time_line(self,_t):


            P1 = self.chart().mapToPosition(QPointF(_t,self.min_y))
            P2 = self.chart().mapToPosition(QPointF(_t,self.max_y))

            self.time_line.setLine(P1.x(),P1.y(),P2.x(),P2.y())

            # recover the values of the three functions
            y1 = np.interp(_t, self.t, self.y1)
            y2 = np.interp(_t, self.t, self.y2)
            y3 = np.interp(_t, self.t, self.y3)

            self.time_changed.emit(_t, y1, y2, y3)


    def mousePressEvent (self,event):

            if (event.button()==QtCore.Qt.LeftButton):

                widgetPos = event.localPos();
                scenePos = self.mapToScene(QtCore.QPoint(int(widgetPos.x()), int(widgetPos.y()) ) );
                chartItemPos = self.chart().mapFromScene(scenePos);
                valueGivenSeries = self.chart().mapToValue(chartItemPos);
                #print("widgetPos:",    widgetPos);
                #print("scenePos:" ,    scenePos);
                #print("chartItemPos:", chartItemPos);
                #print("valSeries:",    valueGivenSeries);
                #print(event.x(),event.y())

                self.update_time_line(valueGivenSeries.x())
                self.time = valueGivenSeries.x()


            QGraphicsView.mousePressEvent (self,event );



    @QtCore.Slot()
    def advance_time(self):

        self.update_time_line(float(self.time))

        if self.time < self.max_t:
            self.time += self.dt
        else:
            self.time = 0
            self.timer.stop()

        return


    def mouseMoveEvent(self, event):

        if (event.buttons()==QtCore.Qt.LeftButton):
            widgetPos = event.localPos();
            scenePos = self.mapToScene(QtCore.QPoint(int(widgetPos.x()), int(widgetPos.y()) ) );
            chartItemPos = self.chart().mapFromScene(scenePos);
            valueGivenSeries = self.chart().mapToValue(chartItemPos);
            self.update_time_line(valueGivenSeries.x())
            self.time = valueGivenSeries.x()

        QGraphicsView.mouseMoveEvent (self, event );
        return



    def mouseReleaseEvent (self,event):

            if (event.button()==QtCore.Qt.LeftButton):
                self.setDragMode(QGraphicsView.NoDrag);
            QGraphicsView.mouseReleaseEvent (self, event );

class point_item(ui.dolomites_item):

    def __init__(self, x, y,color=QtCore.Qt.blue, parent=None):
        ui.dolomites_item.__init__(self,color,parent)

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
        self.setFlag(QGraphicsItem.ItemIsMovable, True );
        #self.setFlag(QGraphicsItem.ItemIgnoresTransformations, True );

        self.set_path();
        self.setZValue(1);
        #info = QString("x: %1; y: %2").arg(x).arg(y);

        #self.xChanged.connect(self.constrain_move)


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
        #self.path.addRect(QtCore.QRectF(-self.dx/2,-self.dy/2,self.dx,self.dy));
        self.path.addEllipse(QtCore.QRectF(-self.dx/2,-self.dy/2,self.dx,self.dy));

        return

    def paint(self, painter, option, widget):
        painter.setPen(self.pen)
        painter.setBrush(self.brushColor)
        painter.drawPath(self.path)

        return

class arrow_item(ui.dolomites_item):

    def __init__(self, P1, P2=0, parent=None):
        ui.dolomites_item.__init__(self)

        self.P1 = P1 # Final point
        self.P2 = P2 # Starting point, default the origin

        #self.P1.xChanged.connect(self.set_path)

        self.pen.setWidthF(2);
        self.setAcceptHoverEvents ( True );



        self.setFlag(QGraphicsItem.ItemIsSelectable, True );
        self.setZValue(1);
        self.set_path();

    @Slot()
    def set_tip(self,x,y):
        self.P1.setX(x)
        self.P1.setY(-y)
        self.set_path()


    @Slot()
    def set_path(self):


        if self.path == None:
             self.path =  QPainterPath();
        else:
            self.path.clear()

        self.path.moveTo(self.P2.x(),-self.P2.y())
        self.path.lineTo(self.P1.x(), self.P1.y())

        self.update(self.boundingRect())

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
        return _path.boundingRect();
