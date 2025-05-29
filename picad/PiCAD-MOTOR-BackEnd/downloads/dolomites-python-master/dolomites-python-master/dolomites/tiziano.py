import datetime
import logging
import numpy as np
# ATTENTION! Python 3.8 Anaconda dist. has not ezdxf package as default options on Windows 10
import ezdxf  # required by save_dxf function (conversion to .dxf file format)
import matplotlib.pyplot as plt
import scipy.optimize as opt
logging.basicConfig(filename='tiziano.log', level=logging.INFO)
import triangle

class point():
    """A point in the drawing, described by its coordinates."""

    def __init__(self, x, y, ph=-1):
        """Init point coordinates and physical tag."""
        self.x = x        # x coordinate
        self.y = y        # y coordinate
        self.ph_tag = ph  # optional physical tag

    def __repr__(self):
        """Point representation for printing purposes."""
        return "point at: (%s, %s, %s)" % (self.x, self.y, self.ph_tag)

    def distance(self, p=None):
        """
        Get the distance between this point (self) and a given point.

        If no second point is given, distance from origin (0, 0) is computed."""
        if p is None:
            return np.sqrt((self.x)**2 + (self.y)**2)
        if not isinstance(p, point):
            print('point.distance() ==> Invalid argument! Input must be a tiziano point.')
            return None

        return np.sqrt((self.x - p.x)**2 + (self.y - p.y)**2)


class label():
    """A label in the drawing, to assign regional (physical) tags."""

    def __init__(self, x, y, ph, area=-1):
        """Init label coordinates and physical tag."""
        self.x = x          # x coordinate
        self.y = y          # y coordinate
        self.ph_tag = ph    # physical tag
        self.area = area    # maximum size allowed for triangle mesh elements

    def __repr__(self):
        """Label representation for printing purposes."""
        return "label at: (%s, %s, %s)" % (self.x, self.y, self.ph_tag)

    def distance(self, p=None):
        """Get the distance between this label (self) and a given point."""
        if p is None:
            return np.sqrt((self.x)**2 + (self.y)**2)
        if not isinstance(p, point):
            print('point.distance() ==> Invalid argument! Input must be a tiziano point.')
            return None

        return np.sqrt((self.x - p.x)**2 + (self.y - p.y)**2)


class hole():
    """A hole label in the drawing, to have no mesh regions."""

    def __init__(self, x, y):
        """Init hole coordinates."""
        self.x = x  # x coordinate
        self.y = y  # y coordinate

    def __repr__(self):
        """Hole label representation for printing purposes."""
        return "hole at: (%s, %s)" % (self.x, self.y)

    def distance(self, p=None):
        """Get the distance between this hole label (self) and a given point."""
        if p is None:
            return np.sqrt((self.x)**2 + (self.y)**2)
        if not isinstance(p, point):
            print('point.distance() ==> Invalid argument! Input must be a tiziano point.')
            return None

        return np.sqrt((self.x - p.x)**2 + (self.y - p.y)**2)


class line():
    """A line in the drawing, described by its endpoints."""

    def __init__(self, n1, n2, ph=-1):
        """Init line endpoints and physical tag."""
        self.n1 = n1      # start point
        self.n2 = n2      # end point
        self.ph_tag = ph  # optional physical tag

    def __repr__(self):
        """Line representation for printing purposes."""
        if self.ph_tag == -1:
            return "line from (%s) to (%s)." % (self.n1, self.n2)

        return "line from (%s) to (%s). Physical tag: %s" % (self.n1, self.n2, self.ph_tag)

    def line_intersection(self, l):
        """Return the intersection point between this line (self) and a given line.

        If there is no intersection returns None.
        """
        tol = 1.e-10

        x = [self.n1.x, self.n2.x]
        y = [self.n1.y, self.n2.y]

        p = [l.n1.x, l.n2.x]
        q = [l.n1.y, l.n2.y]

        det = (p[1] - p[0])*(y[1] - y[0]) - (q[1] - q[0])*(x[1] - x[0])

        if abs(det) < tol:
            return None    # lines are parallel

        t = [((p[1]-p[0])*(q[0]-y[0]) - (q[1]-q[0])*(p[0]-x[0]))/det,
             ((x[1]-x[0])*(q[0]-y[0]) - (y[1]-y[0])*(p[0]-x[0]))/det]

        # check to see if the line segments intersect at a point sufficiently
        # far from the segment endpoints....
        z = True
        if t[0] < tol:
            z = False
        if t[1] < tol:
            z = False
        if t[0] > 1. - tol:
            z = False
        if t[1] > 1.-tol:
            z = False

        if z is False:
            return None

        _x = (x[0]+t[0]*(x[1]-x[0]))
        _y = (y[0]+t[0]*(y[1]-y[0]))
        return point(_x, _y)

    def shortest_distance(self, pt):
        """Get the shortest distance between this line (self) and a given point."""
        if not isinstance(pt, point):
            print('line.shortest_distance() ==> Invalid argument! Input must be a point.')
            return None

        p = pt.x
        q = pt.y
        x = [self.n1.x, self.n2.x]
        y = [self.n1.y, self.n2.y]

        t = ((p-x[0])*(x[1]-x[0]) + (q-y[0])*(y[1]-y[0])) / ((x[1]-x[0])**2 + (y[1]-y[0])**2)

        if t > 1.:
            t = 1.
        if t < 0.:
            t = 0.

        _x = x[0]+t*(x[1]-x[0])
        _y = y[0]+t*(y[1]-y[0])

        return np.sqrt((p - _x)**2 + (q - _y)**2)

    def manhattan_length(self):
        """Return the lenght of the line."""
        return np.sqrt((self.n1.x - self.n2.x)**2 + (self.n1.y - self.n2.y)**2)

    def is_end_point(self, p):
        """Return true if the given point is one of the two endpoints of this line."""
        return p in (self.n1,  self.n2)


class drawing():
    """
    A class to represent a drawing.

    Geometry entities:
    points: a list of points
    lines:  a list of lines (arcs are discretized as lines)
    holes:  a list of holes

    Mesh elements (from pytriangle):
    triangle_points
    triangle_edges
    triangle_triangles
    """

    def __init__(self):
        """Init drawing elements (list of geometry entities and mesh elements)."""
        # geometry entities
        self.points = []
        self.lines  = []
        self.labels = []
        self.holes  = []

        # geometry entities selection list
        self.selected = []

        # mesh elements
        self.triangle_points    = []
        self.triangle_edges     = []
        self.triangle_triangles = []

    def __repr__(self):
        """Drawing representation for printing purposes."""
        res =  """
               drawing with:
               %s points
               %s lines
               %s labels
               %s holes)
               """ % (len(self.points), len(self.lines), len(self.labels), len(self.holes))
        return res

    def plot(self):
        """Plot tiziano drawing using matplotlib.pyplot."""
        _, ax = plt.subplots()

        # points
        for p in self.points:
            if p in self.selected:
                ax.scatter(p.x, p.y, 10, color='orangered', marker="o")
                if p.ph_tag > 0:
                    ax.text(p.x, p.y, str(p.ph_tag), fontsize='small', color='orangered')
            else:
                ax.scatter(p.x, p.y, 10, color='blue', marker="o")
                if p.ph_tag > 0:
                    ax.text(p.x, p.y, str(p.ph_tag), fontsize='small', color='blue')

        # lines
        for l in self.lines:
            if l in self.selected:
                ax.plot([l.n1.x, l.n2.x], [l.n1.y, l.n2.y], color='orangered', lw=1)
                if l.ph_tag > 0:
                    ax.text((l.n1.x+l.n2.x)/2, (l.n1.y+l.n2.y)/2, str(l.ph_tag), color='orangered')
            else:
                ax.plot([l.n1.x, l.n2.x], [l.n1.y, l.n2.y], 'k', lw=1)
                if l.ph_tag > 0:
                    ax.text((l.n1.x+l.n2.x)/2, (l.n1.y+l.n2.y)/2, str(l.ph_tag), color='k')

        # labels
        for lab in self.labels:
            if lab in self.selected:
                ax.scatter(lab.x, lab.y, lab.area*82, marker='o', facecolors='none', edgecolors='orangered')
                ax.text(lab.x, lab.y, lab.ph_tag, color='orangered')
            else:
                ax.scatter(lab.x, lab.y, lab.area*82, marker='o', facecolors='none', edgecolor='darkgreen')
                ax.text(lab.x, lab.y, lab.ph_tag, color='darkgreen')

        # holes
        for h in self.holes:
            if h in self.selected:
                ax.scatter(h.x, h.y, 24, color='orangered', marker='o')
            else:
                ax.scatter(h.x, h.y, 24, color='dimgray', marker='o')

        ax.margins(0.05)
        ax.axis('equal')

        plt.show()

    def add_point(self, p, q=None, ph=-1, coherence=True, get_lns=False):
        """
        Add a point to the drawing.

        A tiziano point can be added to the drawing by using either an already
        defined tiziano point instance (p = point; q = None), or by using
        coordinates (p = x; q = y).
        An optional physical tag can be set.
        If coherence checks are done (default), intersections with other lines
        already in the drawing are searched and if they occur, each intersecting
        line is splitted into two new lines.

        Parameters
        ----------
        p : tiziano point or int or float
            - tiziano point --> if one wants to add a defined tiziano point to
                                the drawing (default behaviour)
            - int or float  --> if one wants to add a tiziano point to the
                                drawing by coordinates definition (p = x)
        q : None or int or float, default None
            - None         --> if one wants to add a defined tiziano point to
                               the drawing (default)
            - int or float --> if one wants to add a tiziano point to current
                               drawing by coordinates definition (q = y)
        ph : int, deafult -1
            Optional physical tag.
        coherence : bool, default True
            Coherence check (intersections with other points and lines).
            True  ==> check coherence (stronger but slower method) (default)
            False ==> do not check coherence (weaker but faster method)
        get_lns : bool, default False
            True  ==> returns new splitted line instances when the point is on
                      a line already in the drawing and coherence check is done
            False ==> returns only the point instance newly generated in the
                      drawing or already on it (default), even if coherence is
                      checked and intersections with lines in the drawing occur

        Returns
        -------
        _p : tiziano point
            Tiziano point instance newly added to the drawing or already in it.
        new_lns : list of tiziano lines (optional)
            Tiziano lines newly generated if coherence is checked and
            intersections occurs with lines already in the drawing.
            Returned only when 'get_lns' is set to True.
        """
        tol = 1e-10  # minimum distance allowed between two points

        # if we have just one argument we are expecting a point instance
        if q is None:
            if not isinstance(p, point):
                print("drawing.add_point() ==> Invalid argument! A point must be passed as single input when q=None.")
                if get_lns is True:
                    return None, []
                return None
            cp = self.closest_point(p.x, p.y)
            if cp[1] < tol:
                print("drawing.add_point() --> ", p, " already in the drawing.")
                if get_lns is True:
                    return cp[0], []
                return cp[0]

            self.points.append(p)
            _p = p

        # if we have two arguments we are expecting coordinates
        else:
            cp = self.closest_point(p, q)  # p = x and q = y
            if cp[1] < tol:
                print("drawing.add_point() --> point at (", p, ',', q, ") already in the drawing.")
                if get_lns is True:
                    return cp[0], []
                return cp[0]

            _p = point(p, q, ph)
            self.points.append(_p)

        # COHERENCE CHECKS =====================================================
        # skip coherence checks
        if coherence is False:
            if get_lns is True:
                return _p, []
            return _p
        # INTERSECTIONS WITH LINES ---------------------------------------------
        # test to see if the point is on an existing line;
        # if so, break intersecting line into two lines;
        to_remove = []  # list for lines to be removed from the drawing
        new_lns = []    # list for new lines to be added to the drawing
        for l in self.lines:
            # for each intersecting line create a couple of splitting lines
            if l.shortest_distance(_p) < tol:
                print("drawing.add_point() --> ", _p, " is on a line. Splitting ", l, "...")
                l1 = line(l.n1, _p, l.ph_tag)
                l2 = line(_p, l.n2, l.ph_tag)
                new_lns.append(l1)
                new_lns.append(l2)
                to_remove.append(l)
        # intersecting lines splitting is performed in two steps
        for l in to_remove:      # first: intersecting lines are removed
            self.lines.remove(l)
        for l in new_lns:        # second: new couples of lines are added
            self.lines.append(l)
        # returns also new splitting lines, if desired
        if get_lns is True:
            return _p, new_lns

        return _p

    def add_line(self, n1, n2, ph=-1, coherence=True, get_pts=False):
        """
        Add a line to the drawing.

        Starting and ending points should be already added to the drawing;
        otherwise, no line is generated.
        An optional physical tag can be set.
        If coherence checks are done (default), intersections with other points
        and lines already in the drawing are searched and if they occur, the
        line is splitted into multiple lines, as well as the other intersecting
        lines.

        Parameters
        ----------
        n1, n2 : tiziano point
            Starting and ending point of the line, respectively.
            They must be already added to the drawing.
        ph : int, optional
            Physical tag (default -1).
        coherence : bool, default True
            Coherence check (intersections with other points and lines).
            True  ==> check coherence (stronger but slower method) (default)
            False ==> do not check coherence (weaker but faster method)
        get_pts : bool, default False
            Flag to return newly generated intersection points
            (from intersections with other lines).
            True  ==> return intersection points with other lines
            False ==> do not return intersection points with other lines (default)

        Returns
        -------
        new_lns : list of tiziano lines
            List of newly generated tiziano lines.
        new_pts : list of tiziano points (optional)
            List of newly generated tiziano points (intersection points with
            other lines). Returned only if coherence is checked and 'get_pts' is
            set to True.
        """
        # check the input type
        if isinstance(n1, point) is False or isinstance(n2, point) is False:
            print("drawing.add_line() ==> Invalid arguments! To add a line please provide two tiziano points as input.")
            if get_pts is True:
                return None, None
            return None

        if n1 not in self.points or n2 not in self.points:
            print("drawing.add_line() ==> Invalid arguments! Endpoints must be already in the drawing.")
            if get_pts is True:
                return None, None
            return None

        # skip dummy lines
        if n1 == n2:
            print("drawing_add_line() --> The two input endpoints are identical. A line cannot be added to drawing.")
            if get_pts is True:
                return [], []
            return []

        # skip lines already added
        for l in self.lines:
            if (l.n1 == n1 and l.n2 == n2) or (l.n1 == n2 and l.n2 == n1):
                print("drawing.add_line() --> ", l, " already in the drawing.")
                # this check is needed for coherence in geom. transformation
                # method: traslation, rotation, mirroring
                if get_pts is True:
                    return [], []
                return []

        _l = line(n1, n2, ph)  # create a line instance from given endpoints
        new_lns = [_l]         # put it into the returned new lines list
        # add line instance to the drawing
        self.lines.append(_l)

        # COHERENCE CHECKS =====================================================
        # skip coherence checks
        if coherence is False:
            if get_pts is True:
                return new_lns, []
            return new_lns

        new_pts = []  # list for returned points to add at line intersections
        tol = 1e-10       # tolerance for intersections computations
        intsc_found = 0   # intersections founded flag
        # INTERSECTIONS WITH POINTS --------------------------------------------
        # check to see if proposed line passes through other points;
        # if so, delete the line and create lines that link intermediate points;
        # does this by recursive use of add_line;
        # a copy of drawing points must be done to get the iterable of this loop
        # (again, this is due to recursive calls to add_line method)
        iter_pts = self.points.copy()
        for p in iter_pts:
            d = _l.shortest_distance(p)
            # skip line endpoints (or points that are too close to them)
            if p.distance(_l.n1) < tol:
                d = 2*tol
            if p.distance(_l.n2) < tol:
                d = 2*tol
            if d < tol:
                print("drawing.add_line() --> ", p, " already in the drawing is on ", _l, ". Splitting the line...")
                intsc_found = 1  # flag to 1 if there is at least one intersect
                self.lines.remove(_l)
                new_lns.remove(_l)
                lns1, pts1 = self.add_line(n1, p, ph, get_pts=True)
                lns2, pts2 = self.add_line(p, n2, ph, get_pts=True)
                new_lns.extend(lns1+lns2)
                new_pts.extend(pts1+pts2)
                break
        # if at least one intersect with points is founded,
        # then searching intersections with lines is useless,
        # because they are already checked inside the previous for loop
        # where the add_line method is called in a recursive way.
        if intsc_found == 1:
            if get_pts is True:
                return new_lns, new_pts
            return new_lns

        # INTERSECTIONS WITH LINES ---------------------------------------------
        # check to see if there are intersections with other lines in the list
        for l in self.lines:
            int_p = l.line_intersection(_l)
            if int_p is not None:
                intsc_found = 1  # flag to 1 if there is at least one intersect
                new_pts.append(int_p)
        # if so, add points at the intersections and split intersecating lines;
        # do this with the help of add_point method
        to_check = []
        for int_p in new_pts:
            _, lns = self.add_point(int_p, get_lns=True)
            new_lns.extend(lns)
            to_check.extend(lns[2:])
        # check to see if some lines in output list are intersecated by new points
        # if so, remove them from output list
        # the follow check can be done using lines shortest_distance and points
        # distance methods as above, but I think that simple comparisons are faster
        for ll in to_check:
            max_x = max(ll.n1.x, ll.n2.x)
            min_x = min(ll.n1.x, ll.n2.x)
            max_y = max(ll.n1.y, ll.n2.y)
            min_y = min(ll.n1.y, ll.n2.y)
            for pp in new_pts:
                if min_x < pp.x < max_x or min_y < pp.y < max_y:
                    new_lns.remove(ll)
                    break
        if intsc_found == 1:    # if an intersection with a line is found
            new_lns.remove(_l)  # original line is removed from returned list

        if get_pts is True:
            return new_lns, new_pts
        return new_lns

    def add_arc(self, n1, n2, angle, step_angle, ph=-1, coherence=True):
        """
        Add an arc to the drawing.

        The arc is immediately discretized.

        n1:         begin
        n2:         end
        angle:      the span angle of the arc in radians
        step_angle: the discretization angle in radians
                    (should be smaller than angle)
        ph:         optional physical tag
        coherence: True  ==> check coherence (default)
                             [intersections with other points and lines]
                   False ==> do not check coherence
        """
        # check the input type
        if step_angle <= 0:
            print("drawing.add_arc() ==> Invalid arguments! To add an arc please provide a positive integer step angle.")
            return None, None

        if isinstance(n1, point) is False or isinstance(n2, point) is False:
            print("drawing.add_arc() ==> Invalid arguments! To add an arc please provide two points as input.")
            return None, None

        if n1 not in self.points or n2 not in self.points:
            print("drawing.add_arc() ==> Invalid arguments! Endpoints must be already in the drawing.")
            return None, None

        # skip dummy arcs
        if n1 == n2:
            print("drawing.add_arc() ==> The two points are identical. An arc cannot be added to the drawing.")
            return [], []

        # only positive angles are considered. If a negative angle is requested
        # we swap the ending points
        if angle < 0:
            angle = -angle
            n1, n2 = n2, n1

        # check if the arc has to be discretized or not
        # if not a line is added
        if step_angle >= abs(angle):
            print("drawing.add_arc() ==> Dummy arc inserted. Please check angles.")
            return [n1, n2], self.add_line(n1, n2, ph, coherence=coherence)

        points = []  # arc points list
        lines  = []  # arc lines list

        # arc chord with respect to origin (0, 0)
        d = np.sqrt((n1.x - n2.x)**2 + (n1.y - n2.y)**2)
        # arc radius with respect to origin (0, 0)
        radius  = 0.5*d/np.sin(angle/2)

        t = [(n2.x - n1.x)/d, (n2.y - n1.y)/d]
        p = [-t[1], t[0]]

        # arc centre coordinates
        center = [n1.x + t[0] * d/2 + p[0] * np.sqrt(radius**2 - 0.25*d**2),
                  n1.y + t[1] * d/2 + p[1] * np.sqrt(radius**2 - 0.25*d**2)]
        # number of segments (lines!) to discretize the arc
        num_seg = int(abs(angle)/step_angle)
        # angle spanned by each segment
        d_angle = angle / num_seg

        # append first point (already in the drawing!) in arc points list
        points.append(n1)
        # add intermediate points to the drawing
        start_angle = np.arctan2(n1.y - center[1], n1.x - center[0]) + d_angle
        for _ in range(num_seg - 1):
            x = radius*np.cos(start_angle) + center[0]
            y = radius*np.sin(start_angle) + center[1]
            pt, lns = self.add_point(x, y, coherence=coherence, get_lns=True)
            points.append(pt)
            lines.extend(lns)
            start_angle += d_angle
        # append last point (already in the drawing!) in arc points list
        points.append(n2)

        # add lines to make the discretized arc
        for i in range(len(points)-1):
            lns, pts = self.add_line(points[i], points[i+1], ph,
                                     coherence=coherence, get_pts=True)
            points.extend(pts)
            lines.extend(lns)

        return points, lines

    def add_label(self, x, y, ph, area):
        """
        Add a label in the current tiziano drawing.

        - x:    x coordinate
        - y:    y coordinate
        - ph:   REQUIRED physical tag
        - area: REQUIRED maximum mesh elements size for triangle mesher
                (for mesh control and refinement)
        """
        l = label(x, y, ph, area)
        self.labels.append(l)

        return l

    def add_hole(self, x, y):
        """
        Add a hole in the current tiziano drawing.

        - x coordinate
        - y coordinate
        """
        h = hole(x, y)
        self.holes.append(h)
        return h

    def closest_point(self, x, y):
        """Return the closest point in the drawing for the given (x, y) coordinates."""
        distance_min = 1e10
        _pp = None

        for _p in self.points:
            distance = np.sqrt(((_p.x-x)**2)+((_p.y-y)**2))
            # print('checking point', _p,distance)
            if distance < distance_min:
                distance_min = distance
                _pp = _p
        return _pp, distance_min

    def select(self, item_list):
        """Add geometry entities to selection list."""
        if not isinstance(item_list, list):
            print('drawing.select() ==> Invaid argument! A list of tiziano geometry entities must be passed.')
            return None

        geom_ent = self.points + self.lines + self.labels + self.holes
        to_select = []
        for item in item_list:
            if isinstance(item, (point, line, label, hole)):
                if item in geom_ent:
                    to_select.append(item)
                    if item not in self.selected:
                        self.selected.append(item)
        return to_select

    def select_all(self, entity=-1):
        """
        Select all geometry entities in drawing of a specified type.

        - entity: geometry entity types to select in the drawing
                  * -1: all (default)
                  * 0:  all points
                  * 1:  all lines (without their endpoints!)
                  * 2:  all labels and holes
        """
        to_select = []
        if entity in [-1, 0]:
            to_select.extend(self.select(self.points))
        if entity in [-1, 1]:
            to_select.extend(self.select(self.lines))
        if entity in [-1, 2]:
            to_select.extend(self.select(self.labels + self.holes))

        return to_select

    def select_rect(self, xb, yb, a, b, entity=-1, mode='<'):
        """Select geometry entities inside/outside a defined rectangle."""
        xb_1 = xb + a
        yb_1 = yb + b
        x_max = max(xb, xb_1)
        y_max = max(yb, yb_1)
        x_min = min(xb, xb_1)
        y_min = min(yb, yb_1)
        to_select = []
        if mode == '<':
            if entity in [-1, 0]:
                for pp in self.points:
                    if x_min <= pp.x <= x_max and y_min <= pp.y <= y_max:
                        to_select.append(pp)
            if entity in [-1, 1]:
                for ll in self.lines:
                    if (
                        x_min <= ll.n1.x <= x_max and y_min <= ll.n1.y <= y_max and
                        x_min <= ll.n2.x <= x_max and y_min <= ll.n2.y <= y_max
                       ):
                        to_select.append(ll)
            if entity in [-1, 2]:
                for lab in self.labels:
                    if x_min <= lab.x <= x_max and y_min <= lab.y <= y_max:
                        to_select.append(lab)
                for hh in self.holes:
                    if x_min <= hh.x <= x_max and y_min <= hh.y <= y_max:
                        to_select.append(hh)
        elif mode == '>':
            if entity in [-1, 0]:
                for pp in self.points:
                    if (pp.x < x_min or pp.x > x_max) and (pp.y < y_min or pp.y > y_max):
                        to_select.append(pp)
            if entity in [-1, 1]:
                for ll in self.lines:
                    if (
                        ((ll.n1.x < x_min or ll.n1.x > x_max) and (ll.n1.y < y_min or ll.n1.y > y_max))
                        and
                        ((ll.n2.x < x_min or ll.n2.x > x_max) and (ll.n2.y < y_min or ll.n2.y > y_max))
                       ):
                        to_select.append(ll)
            if entity in [-1, 2]:
                for lab in self.labels:
                    if (lab.x < x_min or lab.x > x_max) and (lab.y < y_min or lab.y > y_max):
                        to_select.append(lab)
                for hh in self.holes:
                    if (hh.x < x_min or hh.x > x_max) and (hh.y < y_min or hh.y > y_max):
                        to_select.append(hh)
        else:
            print('drawing.select_rect() ==> Invalid argument. Mode must be "<" or ">".')

        self.selected.extend([item for item in to_select if item not in self.selected])
        return to_select

    def select_circ(self, r, xc=0, yc=0, R=-1, entity=-1, mode='<'):
        """Select geometry entities inside/outside a defined circle."""
        c = point(xc, yc)
        to_select = []
        if mode == '<':
            if entity in [-1, 0]:
                for pp in self.points:
                    if pp.distance(c) <= r:
                        to_select.append(pp)
            if entity in [-1, 1]:
                for ll in self.lines:
                    if ll.n1.distance(c) <= r and ll.n2.distance(c) <= r:
                        to_select.append(ll)
            if entity in [-1, 2]:
                for lab in self.labels:
                    if lab.distance(c) <= r:
                        to_select.append(lab)
                for hh in self.holes:
                    if hh.distance(c) <= r:
                        to_select.append(hh)
        elif mode == '>':
            if entity in [-1, 0]:
                for pp in self.points:
                    if pp.distance(c) > r:
                        to_select.append(pp)
            if entity in [-1, 1]:
                for ll in self.lines:
                    if ll.n1.distance(c) > r and ll.n2.distance(c) > r:
                        to_select.append(ll)
            if entity in [-1, 2]:
                for lab in self.labels:
                    if lab.distance(c) > r:
                        to_select.append(lab)
                for hh in self.holes:
                    if hh.distance(c) > r:
                        to_select.append(hh)
        elif mode == '<>':
            if R <= r:
                print("drawing.select_circ() --> Please specify a valid range to select points with mode '<>'.")
            if entity in [-1, 0]:
                for pp in self.points:
                    if r < pp.distance(c) < R:
                        to_select.append(pp)
            if entity in [-1, 1]:
                for ll in self.lines:
                    if r < ll.n1.distance(c) < R and r < ll.n2.distance(c) < R:
                        to_select.append(ll)
            if entity in [-1, 2]:
                for lab in self.labels:
                    if r < lab.distance(c) < R:
                        to_select.append(lab)
                for hh in self.holes:
                    if r < hh.distance(c) < R:
                        to_select.append(hh)
        else:
            print("drawing.select_circ() ==> Invalid value for kwarg mode.")

        self.selected.extend([item for item in to_select if item not in self.selected])
        return to_select

    def clear_selected(self):
        """Clear geometry entities selection list of the drawing."""
        self.selected.clear()

    def clear(self, entity=-1):
        """Remove geometry elements from the drawing.

        - entity: geometry entity types to remove from drawing
                  * -1: all (default)
                  * 0:  all points and lines
                  * 1:  only all lines (without their endpoints!)
                  * 2:  all labels and holes
        """
        if entity in [-1, 0]:
            self.points.clear()
            self.lines.clear()
        if entity in [-1, 1]:
            self.lines.clear()
        if entity in [-1, 2]:
            self.labels.clear()
            self.holes.clear()

    def remove_point(self, p, mode='<', r=-1):
        """
        Remove a point and all the associated lines from the drawing.
        
        All the lines that have as endpoint one of the given points are
        removed.

        Parameters
        ----------
        p: list of tiziano points or int or float
            If a list of point() instances is given the points are removed
            with their associated lines.
            If a real number is given (int or float), it is assumed to be
            a radius of a circle used to decide which points to remove.
        mode : '<'  (default) the point inside the radius are removed
               '>'  the point outside the radius are removed
               '<>' the point within p and r are removed
        r : int or float
            External radius in case of mode='<>'
        """
        _to_remove = []

        # if we get a list of points
        if isinstance(p, list):
            for pp in p:
                if isinstance(pp, point):
                    to_remove = []
                    # check if there are lines connected to the point and delete them
                    for l in self.lines:
                        if l.is_end_point(pp):
                            to_remove.append(l)
                    for l in to_remove:
                        self.lines.remove(l)
                        if l in self.selected:
                            self.selected.remove(l)
                    # remove the point from current drawing
                    self.points.remove(pp)
                    if pp in self.selected:
                        self.selected.remove(pp)

        # if we get a number we consider it as a radius
        elif isinstance(p, (float, int)):
            if mode == '<':
                # select all the points inside the radius
                for _p in self.points:
                    if _p.distance() <= p:
                        _to_remove.append(_p)
            elif mode == '>':
                # select all the points outside the radius
                for _p in self.points:
                    if _p.distance() >= p:
                        _to_remove.append(_p)
            elif mode == '<>':
                # select all the points within p and r
                if r <= p:
                    print("drawing.remove_point() --> Please specify a valid range to remove points with mode '<>'.")
                for _p in self.points:
                    if p <= _p.distance() <= r:
                        _to_remove.append(_p)

            # remove selected points with recursive use of remove_point function
            self.remove_point(_to_remove)

    def copy_traslate(self, items, X, Y, N=1, entity=-1, coherence=True):
        """
        Copy traslate tiziano geometry entities in the drawing.

        A number of copies different from one (deafult) can be specified with
        the 'N' kwarg, as well as which type of geom. entities in 'items' to
        copy traslate with the 'entity' kwarg.
        Coherence can be checked (default) or not during copy traslations.

        Parameters
        ----------
        items : list of tiziano geometry entities [points, lines, labels, holes]
            List of geom. entities to copy traslate in the drawing.
        X, Y : int or float
            Offsets along x and y axis respectively for traslation
        N : int, default 1
            Number of copies.
        entity : int, default -1
            Geom. entity types in items to copy traslate
            -1; all (default)
             0; only points
             1; only lines with thei endpoints
             2; only labels and holes
        coherence : bool, default True
            Coherence check (intersections with other points and lines).
            True  ==> check coherence (stronger but slower method) (default)
            False ==> do not check coherence (weaker but faster method)
        """
        # check items parameter input type
        if not isinstance(items, list):
            print("drawing.copy_traslate() ==> Invalid argument. A list of tiziano geom. entities must be used.")
            return None, None, None, None
        # check input type for entity parameter
        if entity not in [-1, 0, 1, 2]:
            print('drawing.copy_traslate() ==> Invalid argument! Entity must be -1, 0, 1, 2.')
            return None, None, None, None

        # initialize lists of geom. entities to return--------------------------
        points = []
        lines  = []
        labels = []
        holes  = []
        # get number of points and lines already in the drawing-----------------
        # Counter variables are used inside the main loop to keep track of those
        # objects that are newly created and so they have to be returned.
        # This check is applied only to points and lines because in the add_
        # methods of labels and holes no checks for existance in current drawing
        # are implemented yet.
        count_pts = len(self.points)
        count_lns = len(self.lines)
        for n in range(1, N+1):  # loop on copies ------------------------------
            # compute offsets for each copy
            xx = n*X
            yy = n*Y
            for item in items:   # loop on tiziano geometry entities -----------
                if entity in [-1, 0]:  # copy points
                    if isinstance(item, point):
                        pt, lns = self.add_point(item.x + xx, item.y + yy,
                                                 item.ph_tag,
                                                 coherence=coherence,
                                                 get_lns=True)
                        # add pt to output points list only if pt is newly created
                        new_count_pts = len(self.points)
                        if new_count_pts > count_pts:
                            points.append(pt)
                            count_pts = new_count_pts
                        # add lns to output lines list
                        lines.extend(lns)
                        count_lns = len(self.lines)
                if entity in [-1, 1]:  # copy lines
                    if isinstance(item, line):
                        new_n1, lns1 = self.add_point(item.n1.x+xx, item.n1.y+yy,
                                                      item.n1.ph_tag,
                                                      coherence=coherence,
                                                      get_lns=True)
                        # add new line endpoints to output points list
                        # only if they are newly created
                        new_count_pts = len(self.points)
                        if new_count_pts > count_pts:
                            points.append(new_n1)
                            count_pts = new_count_pts
                        new_n2, lns2 = self.add_point(item.n2.x+xx, item.n2.y+yy,
                                                      item.n2.ph_tag,
                                                      coherence=coherence,
                                                      get_lns=True)
                        new_count_pts = len(self.points)
                        if new_count_pts > count_pts:
                            points.append(new_n2)
                            count_pts = new_count_pts
                        # add splitted lines to output list
                        lines.extend(lns1 + lns2)
                        count_lns = len(self.lines)
                        lns, pts = self.add_line(new_n1, new_n2, item.ph_tag,
                                                 coherence=coherence,
                                                 get_pts=True)
                        # add lines to ouput list only if they are newly created
                        new_count_lns = len(self.lines)
                        if new_count_lns > count_lns:
                            lines.extend(lns)
                            count_lns = new_count_lns
                        # add intersecting points to output list
                        points.extend(pts)
                        count_pts = len(self.points)
                if entity in [-1, 2]:  # copy labels and holes
                    if isinstance(item, label):
                        labels.append(self.add_label(item.x + xx, item.y + yy,
                                                     item.ph_tag, item.area))
                    elif isinstance(item, hole):
                        holes.append(self.add_hole(item.x + xx, item.y + yy))
                if not isinstance(item, (point, line, label, hole)):
                    print('drawing.copy_traslate() --> Skip element in list. Element is not a tiziano geometry entity.')

        # final check not needed if coherence was not checked when new
        # points and lines are added to the drawing
        if coherence is False:
            return points, lines, labels, holes

        # FINAL COHERENCE CHECK
        # to return lines actually generated after copy-transformation
        tol = 1e-10
        to_check = lines.copy()
        for ll in to_check:
            for p in points:
                d = ll.shortest_distance(p)
                # skip line endpoints (or points that are too close to them)
                if p.distance(ll.n1) < tol:
                    d = 2*tol
                if p.distance(ll.n2) < tol:
                    d = 2*tol
                if d < tol:
                    lines.remove(ll)
                    break

        return points, lines, labels, holes

    def copy_rotate(self, items, angle, deg=False, xc=0, yc=0, N=1, entity=-1, coherence=True):
        """Copy rotate tiziano geometry entities in the drawing.

        Positive angles refer to counterclockwise rotations.
        Rotation angles should be specified in radians as default behaviour, but
        inputs in degrees are also possible by setting the kwarg 'deg' to True.
        Rotation is performed around origin point (0, 0) as default behaviour.
        Rotation can be specified around a point different from origin by
        acting on 'xc' and 'yc' kwargs.
        A number of copies different from one (deafult) can be specified with
        the 'N' kwarg, as well as which type of geom. entities in 'items' to
        copy rotate with the 'entity' kwarg.
        Coherence can be checked (default) or not during copy rotations.

        Parameters
        ----------
        items : list of tiziano geom. entities [points, lines, labels, holes]
            List of tiziano geom. entities to copy rotate
        angle : int or float
            Rotation angle (in radians as default behaviour)
        deg : bool, default False
            True; 'angle' is assumed to be in degrees
            False; (deafult) 'angle' is assumed to be in radians
        xc, yc : int or float, default 0
            Coordinates of rotation centre (default origin).
        N : int, default 1
            Number of copies.
        entity : int, default -1
            Geom. entity types in 'items' to copy rotate
                  -1; all (default)
                   0;  only points
                   1;  only lines with their endpoints
                   2;  only labels and holes
        coherence : bool, default True
            Coherence check (intersections with other points and lines).
            True  ==> check coherence (stronger but slower method) (default)
            False ==> do not check coherence (weaker but faster method)

        Returns
        -------
        points : list of tiziano points
            New points generated by copy rotation action.
        lines : list of tiziano lines
            New lines generated by copy rotation action.
        labels : list of tiziano labels
            New labels generated by copy rotation action.
        holes : list of tiziano holes
            New holes generated by copy rotation action.
        """
        # check items parameter input type
        if not isinstance(items, list):
            print("drawing.copy_rotate() ==> Invalid argument. A list of tiziano geom. entities must be used.")
            return None, None, None, None
        # check input type for entity parameter
        if entity not in [-1, 0, 1, 2]:
            print('drawing.copy_rotate() ==> Invalid argument! Entity must be -1, 0, 1, 2.')
            return None, None, None, None

        if deg is True:                 # if rotation angle is given in degs
            angle = angle*np.pi/180     # convert it to rads
        elif deg is False:              # if rotation angle is given in rads
            pass                        # maintain it in rads
        else:
            print('drawing.copy_rotate() ==> Invalid argument! deg must be True or False.')
            return None, None, None, None

        xy_c = np.array([xc, yc])  # array of rotation centre coordinates
        # initialize lists of geom. entities to return--------------------------
        points = []
        lines  = []
        labels = []
        holes  = []
        # get number of points and lines already in the drawing-----------------
        # Counter variables are used inside the main loop to keep track of those
        # objects that are newly created and so they have to be returned.
        # This check is applied only to points and lines because in the add_
        # methods of labels and holes no checks for existance in current drawing
        # are implemented yet.
        count_pts = len(self.points)
        count_lns = len(self.lines)
        for n in range(1, N+1):  # loop on copies ------------------------------
            # compute rotation matrix for each copy
            rot_mat = np.array([[np.cos(n*angle), -np.sin(n*angle)],
                                [np.sin(n*angle),  np.cos(n*angle)]])
            # loop over geom. entities to copy rotate --------------------------
            # New coordinates of each geom. element to rotate are computed using
            # rotation transformation. Then, add_ methods are called and newly
            # generated objects are appended to output points and lines lists.
            for item in items:
                if entity in [-1, 0]:  # copy points
                    if isinstance(item, point):
                        xy = np.array([item.x, item.y])
                        new_xy = (rot_mat @ xy) - (rot_mat @ xy_c) + xy_c
                        pt, lns = self.add_point(new_xy[0], new_xy[1],
                                                 item.ph_tag,
                                                 coherence=coherence,
                                                 get_lns=True)
                        # add pt to output points list only if pt is newly created
                        new_count_pts = len(self.points)
                        if new_count_pts > count_pts:
                            points.append(pt)
                            count_pts = new_count_pts
                        # add splitted lines to output list
                        lines.extend(lns)
                        count_lns = len(self.lines)
                if entity in [-1, 1]:  # copy lines
                    if isinstance(item, line):
                        xy1 = np.array([item.n1.x, item.n1.y])
                        new_xy1 = (rot_mat @ xy1) - (rot_mat @ xy_c) + xy_c
                        xy2 = np.array([item.n2.x, item.n2.y])
                        new_xy2 = (rot_mat @ xy2) - (rot_mat @ xy_c) + xy_c
                        new_n1, lns1 = self.add_point(new_xy1[0], new_xy1[1],
                                                      item.n1.ph_tag,
                                                      coherence=coherence,
                                                      get_lns=True)
                        # add new line endpoints to output points list
                        # only if they are newly created
                        new_count_pts = len(self.points)
                        if new_count_pts > count_pts:
                            points.append(new_n1)
                            count_pts = new_count_pts
                        new_n2, lns2 = self.add_point(new_xy2[0], new_xy2[1],
                                                      item.n2.ph_tag,
                                                      coherence=coherence,
                                                      get_lns=True)
                        new_count_pts = len(self.points)
                        if new_count_pts > count_pts:
                            points.append(new_n2)
                            count_pts = new_count_pts
                        # add splitted lines to output list
                        lines.extend(lns1 + lns2)
                        count_lns = len(self.lines)
                        lns, pts = self.add_line(new_n1, new_n2, item.ph_tag,
                                                 coherence=coherence,
                                                 get_pts=True)
                        # add lns to output list only if ln is newly created
                        new_count_lns = len(self.lines)
                        if new_count_lns > count_lns:
                            lines.extend(lns)
                            count_lns = new_count_lns
                        # add intersecting points to output list
                        points.extend(pts)
                        count_pts = len(self.points)
                if entity in [-1, 2]:  # copy labels and holes
                    if isinstance(item, label):
                        xy = np.array([item.x, item.y])
                        new_xy = (rot_mat @ xy) - (rot_mat @ xy_c) + xy_c
                        labels.append(self.add_label(new_xy[0], new_xy[1],
                                                     item.ph_tag, item.area))
                    elif isinstance(item, hole):
                        xy = np.array([item.x, item.y])
                        new_xy = (rot_mat @ xy) - (rot_mat @ xy_c) + xy_c
                        holes.append(self.add_hole(new_xy[0], new_xy[1]))
                if not isinstance(item, (point, line, label, hole)):
                    print('drawing.copy_rotate() --> Skip element in list. Element is not a tiziano geometry entity.')

        # final check not needed if coherence was not checked when new
        # points and lines are added to the drawing
        if coherence is False:
            return points, lines, labels, holes

        # FINAL COHERENCE CHECK
        # to return lines actually generated after copy-transformation
        tol = 1e-10
        to_check = lines.copy()
        for ll in to_check:
            for p in points:
                d = ll.shortest_distance(p)
                # skip line endpoints (or points that are too close to them)
                if p.distance(ll.n1) < tol:
                    d = 2*tol
                if p.distance(ll.n2) < tol:
                    d = 2*tol
                if d < tol:
                    lines.remove(ll)
                    break

        return points, lines, labels, holes

    def move_traslate(self, items, X, Y, entity=-1, coherence=True):
        """
        Move traslate tiziano geometry entities in the drawing.

        Geometry entities type in 'items' to move traslate can be specified with
        the 'entity' kwarg.
        Coherence can be checked (default) or not during move traslations.

        Parameters
        ----------
        items : list of tiziano geometry entities [points, lines, labels, holes]
            List of geom. entities to move traslate in the drawing.
        X, Y : int or float
            Offsets along x and y axis respectively for traslation
        entity : int, default -1
            Geom. entity types in items to move traslate
            -1; all (default)
             0; only points
             1; only lines with thei endpoints
             2; only labels and holes
        coherence : bool, default True
            Coherence check (intersections with other points and lines).
            True  ==> check coherence (stronger but slower method) (default)
            False ==> do not check coherence (weaker but faster method)
        """
        # check items parameter input type
        if not isinstance(items, list):
            print("drawing.move_traslate() ==> Invalid argument. A list of tiziano geom. entities must be used.")
            return None, None, None, None
        # check input type for entity parameter
        if entity not in [-1, 0, 1, 2]:
            print('drawing.move_traslate() ==> Invalid argument! Entity must be -1, 0, 1, 2.')
            return None, None, None, None

        # initialize lists of geom. entities to return--------------------------
        points = []
        lines  = []
        labels = []
        holes  = []
        
        # init lists for points and lines to be moved
        pts_tomove = []      # points to be moved
        lns_tomove_n1 = []   # lines wit only first endp to be moved
        lns_tomove_n2 = []   # lines with only second endp to be moved
        lns_tomove_n12 = []  # lines with both endp to be moved
        # loop on input tiziano geometry entities to collect points and lines
        # to be moved
        for item in items:
            if entity in [-1, 0]:  # consider points in items
                if isinstance(item, point):
                    if item not in pts_tomove:
                        pts_tomove.append(item)  # store point to be moved
                    for ll in self.lines:        # store associated lines
                        if ll.n1 == item:        # only first endp
                            if ll in lns_tomove_n2:
                                lns_tomove_n2.remove(ll)
                                lns_tomove_n12.append(ll)
                            elif ll in lns_tomove_n12:
                                pass
                            elif ll not in lns_tomove_n1:
                                lns_tomove_n1.append(ll)
                            else:
                                pass
                        if ll.n2 == item:       # only second endp
                            if ll in lns_tomove_n1:
                                lns_tomove_n1.remove(ll)
                                lns_tomove_n12.append(ll)
                            elif ll in lns_tomove_n12:
                                pass
                            elif ll not in lns_tomove_n2:
                                lns_tomove_n2.append(ll)
                            else:
                                pass

            if entity in [-1, 1]:  # consider lines in items
                if isinstance(item, line) and item not in lns_tomove_n12:
                    lns_tomove_n12.append(item)
                    if item.n1 not in pts_tomove:
                        pts_tomove.append(item.n1)
                        for ll in self.lines:        # store associated lines
                            if ll.n1 == item.n1:        # only first endp
                                if ll in lns_tomove_n2:
                                    lns_tomove_n2.remove(ll)
                                    lns_tomove_n12.append(ll)
                                elif ll in lns_tomove_n12:
                                    pass
                                elif ll not in lns_tomove_n1:
                                    lns_tomove_n1.append(ll)
                                else:
                                    pass
                            if ll.n2 == item.n1:       # only second endp
                                if ll in lns_tomove_n1:
                                    lns_tomove_n1.remove(ll)
                                    lns_tomove_n12.append(ll)
                                elif ll in lns_tomove_n12:
                                    pass
                                elif ll not in lns_tomove_n2:
                                    lns_tomove_n2.append(ll)
                                else:
                                    pass
                    if item.n2 not in pts_tomove:
                        pts_tomove.append(item.n2)
                        for ll in self.lines:        # store associated lines
                            if ll.n1 == item.n2:        # only first endp
                                if ll in lns_tomove_n2:
                                    lns_tomove_n2.remove(ll)
                                    lns_tomove_n12.append(ll)
                                elif ll in lns_tomove_n12:
                                    pass
                                elif ll not in lns_tomove_n1:
                                    lns_tomove_n1.append(ll)
                                else:
                                    pass
                            if ll.n2 == item.n2:       # only second endp
                                if ll in lns_tomove_n1:
                                    lns_tomove_n1.remove(ll)
                                    lns_tomove_n12.append(ll)
                                elif ll in lns_tomove_n12:
                                    pass
                                elif ll not in lns_tomove_n2:
                                    lns_tomove_n2.append(ll)
                                else:
                                    pass

            if entity in [-1, 2]:  # consider labels and holes in items
                # labels and holes can be moved right here
                # (no problem of coherence checks)
                if isinstance(item, label):
                    self.labels.remove(item)
                    labels.append(self.add_label(item.x + X, item.y + Y,
                                                 item.ph_tag, item.area))
                elif isinstance(item, hole):
                    self.holes.remove(item)
                    holes.append(self.add_hole(item.x + X, item.y + Y))
            if not isinstance(item, (point, line, label, hole)):
                print('drawing.move_traslate() --> Skip element in list. Element is not a tiziano geometry entity.')

        # remove points that have to be moved with their associated lines
        self.remove_point(pts_tomove)
        # remove lines 
        for ll in lns_tomove_n1 + lns_tomove_n2 + lns_tomove_n12:
            if ll in self.lines:
                self.lines.remove(ll)
        # get number of points and lines now in the drawing --------------------
        # Counter variables to keep track of those objects that are newly
        # created and so they have to be returned in the end.
        # This check is applied only to points and lines because in the add_
        # methods of labels and holes no checks for existance in current drawing
        # are implemented yet.
        count_pts = len(self.points)
        count_lns = len(self.lines)
        # move points
        for pp in pts_tomove:
            pt, lns = self.add_point(pp.x + X, pp.y + Y, pp.ph_tag,
                                     coherence=coherence, get_lns=True)
            # add pt to output points list only if pt is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(pt)
                count_pts = new_count_pts
            # add lns to output lines list
            lines.extend(lns)
            count_lns = len(self.lines)

        # move lines
        for ll in lns_tomove_n1:  # lines with only first endp to move
            new_n1, lns = self.add_point(ll.n1.x+X, ll.n1.y+Y, ll.n1.ph_tag,
                                         coherence=coherence, get_lns=True)
            # add new line endpoint to output points list
            # only if it is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(new_n1)
                count_pts = new_count_pts
            # add lns to output lines list
            lines.extend(lns)
            count_lns = len(self.lines)
            lns, pts = self.add_line(new_n1, ll.n2, ll.ph_tag,
                                     coherence=coherence, get_pts=True)
            # add lines to ouput list only if they are newly created
            new_count_lns = len(self.lines)
            if new_count_lns > count_lns:
                lines.extend(lns)
                count_lns = new_count_lns
            # add intersecting points to output list
            points.extend(pts)
            count_pts = len(self.points)
        for ll in lns_tomove_n2:  # lines with only second endp to move
            new_n2, lns = self.add_point(ll.n2.x+X, ll.n2.y+Y, ll.n2.ph_tag,
                                         coherence=coherence, get_lns=True)
            # add new line endpoint to output points list
            # only if it is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(new_n2)
                count_pts = new_count_pts
            # add lns to output lines list
            lines.extend(lns)
            count_lns = len(self.lines)
            lns, pts = self.add_line(ll.n1, new_n2, ll.ph_tag,
                                     coherence=coherence, get_pts=True)
            # add lines to ouput list only if they are newly created
            new_count_lns = len(self.lines)
            if new_count_lns > count_lns:
                lines.extend(lns)
                count_lns = new_count_lns
            # add intersecting points to output list
            points.extend(pts)
            count_pts = len(self.points)
        for ll in lns_tomove_n12:  # lines with both endp to move
            new_n1, lns1 = self.add_point(ll.n1.x+X, ll.n1.y+Y, ll.n1.ph_tag,
                                       coherence=coherence, get_lns=True)
            # add new line endpoint to output points list
            # only if it is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(new_n1)
                count_pts = new_count_pts
            new_n2, lns2 = self.add_point(ll.n2.x+X, ll.n2.y+Y, ll.n2.ph_tag,
                                       coherence=coherence, get_lns=True)
            # add new line endpoint to output points list
            # only if it is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(new_n2)
                count_pts = new_count_pts
            # add splitted lines to output list
            lines.extend(lns1 + lns2)
            count_lns = len(self.lines)
            lns, pts = self.add_line(new_n1, new_n2, ll.ph_tag,
                                     coherence=coherence, get_pts=True)
            # add lines to ouput list only if they are newly created
            new_count_lns = len(self.lines)
            if new_count_lns > count_lns:
                lines.extend(lns)
                count_lns = new_count_lns
            # add intersecting points to output list
            points.extend(pts)
            count_pts = len(self.points)

        # final check not needed if coherence was not checked when new
        # points and lines are added to the drawing
        if coherence is False:
            return points, lines, labels, holes

        # FINAL COHERENCE CHECK
        # to return lines actually generated after move-transformation
        tol = 1e-10
        to_check = lines.copy()
        for ll in to_check:
            for p in points:
                d = ll.shortest_distance(p)
                # skip line endpoints (or points that are too close to them)
                if p.distance(ll.n1) < tol:
                    d = 2*tol
                if p.distance(ll.n2) < tol:
                    d = 2*tol
                if d < tol:
                    lines.remove(ll)
                    break

        return points, lines, labels, holes

    def move_rotate(self, items, angle, deg=False, xc=0, yc=0, entity=-1, coherence=True):
        """Move rotate tiziano geometry entities in the drawing.

        Positive angles refer to counterclockwise rotations.
        Rotation angles should be specified in radians as default behaviour, but
        inputs in degrees are also possible by setting the kwarg 'deg' to True.
        Rotation is performed around origin point (0, 0) as default behaviour.
        Rotation can be specified around a point different from origin by
        acting on 'xc' and 'yc' kwargs.
        Coherence can be checked (default) or not during move rotation.

        Parameters
        ----------
        items : list of tiziano geom. entities [points, lines, labels, holes]
            List of tiziano geom. entities to move rotate
        angle : int or float
            Rotation angle (in radians as default behaviour)
        deg : bool, default False
            True; 'angle' is assumed to be in degrees
            False; (deafult) 'angle' is assumed to be in radians
        xc, yc : int or float, default 0
            Coordinates of rotation centre (default origin).
        entity : int, default -1
            Geom. entity types in 'items' to move rotate
                  -1; all (default)
                   0;  only points
                   1;  only lines with their endpoints
                   2;  only labels and holes
        coherence : bool, default True
            Coherence check (intersections with other points and lines).
            True  ==> check coherence (stronger but slower method) (default)
            False ==> do not check coherence (weaker but faster method)

        Returns
        -------
        points : list of tiziano points
            New points generated by move rotation action.
        lines : list of tiziano lines
            New lines generated by move rotation action.
        labels : list of tiziano labels
            New labels generated by move rotation action.
        holes : list of tiziano holes
            New holes generated by move rotation action.
        """
        # check items parameter input type
        if not isinstance(items, list):
            print("drawing.move_rotate() ==> Invalid argument. A list of tiziano geom. entities must be used.")
            return None, None, None, None
        # check input type for entity parameter
        if entity not in [-1, 0, 1, 2]:
            print('drawing.move_rotate() ==> Invalid argument! Entity must be -1, 0, 1, 2.')
            return None, None, None, None

        if deg is True:                 # if rotation angle is given in degs
            angle = angle*np.pi/180     # convert it to rads
        elif deg is False:              # if rotation angle is given in rads
            pass                        # maintain it in rads
        else:
            print('drawing.move_rotate() ==> Invalid argument! deg must be True or False.')
            return None, None, None, None

        # initialize lists of geom. entities to return--------------------------
        points = []
        lines  = []
        labels = []
        holes  = []
        
        # init lists for points and lines to be moved
        pts_tomove = []      # points to be moved
        lns_tomove_n1 = []   # lines wit only first endp to be moved
        lns_tomove_n2 = []   # lines with only second endp to be moved
        lns_tomove_n12 = []  # lines with both endp to be moved
        # compute rotation matrix
        rot_mat = np.array([[np.cos(angle), -np.sin(angle)],
                            [np.sin(angle),  np.cos(angle)]])
        xy_c = np.array([xc, yc])  # array of rotation centre coordinates
        # loop on input tiziano geometry entities to collect points and lines
        # to be moved
        for item in items:
            if entity in [-1, 0]:  # consider points in items
                if isinstance(item, point):
                    if item not in pts_tomove:
                        pts_tomove.append(item)  # store point to be moved
                    for ll in self.lines:        # store associated lines
                        if ll.n1 == item:        # only first endp
                            if ll in lns_tomove_n2:
                                lns_tomove_n2.remove(ll)
                                lns_tomove_n12.append(ll)
                            elif ll in lns_tomove_n12:
                                pass
                            elif ll not in lns_tomove_n1:
                                lns_tomove_n1.append(ll)
                            else:
                                pass
                        if ll.n2 == item:       # only second endp
                            if ll in lns_tomove_n1:
                                lns_tomove_n1.remove(ll)
                                lns_tomove_n12.append(ll)
                            elif ll in lns_tomove_n12:
                                pass
                            elif ll not in lns_tomove_n2:
                                lns_tomove_n2.append(ll)
                            else:
                                pass

            if entity in [-1, 1]:  # consider lines in items
                if isinstance(item, line) and item not in lns_tomove_n12:
                    lns_tomove_n12.append(item)
                    if item.n1 not in pts_tomove:
                        pts_tomove.append(item.n1)
                        for ll in self.lines:        # store associated lines
                            if ll.n1 == item.n1:        # only first endp
                                if ll in lns_tomove_n2:
                                    lns_tomove_n2.remove(ll)
                                    lns_tomove_n12.append(ll)
                                elif ll in lns_tomove_n12:
                                    pass
                                elif ll not in lns_tomove_n1:
                                    lns_tomove_n1.append(ll)
                                else:
                                    pass
                            if ll.n2 == item.n1:       # only second endp
                                if ll in lns_tomove_n1:
                                    lns_tomove_n1.remove(ll)
                                    lns_tomove_n12.append(ll)
                                elif ll in lns_tomove_n12:
                                    pass
                                elif ll not in lns_tomove_n2:
                                    lns_tomove_n2.append(ll)
                                else:
                                    pass
                    if item.n2 not in pts_tomove:
                        pts_tomove.append(item.n2)
                        for ll in self.lines:        # store associated lines
                            if ll.n1 == item.n2:        # only first endp
                                if ll in lns_tomove_n2:
                                    lns_tomove_n2.remove(ll)
                                    lns_tomove_n12.append(ll)
                                elif ll in lns_tomove_n12:
                                    pass
                                elif ll not in lns_tomove_n1:
                                    lns_tomove_n1.append(ll)
                                else:
                                    pass
                            if ll.n2 == item.n2:       # only second endp
                                if ll in lns_tomove_n1:
                                    lns_tomove_n1.remove(ll)
                                    lns_tomove_n12.append(ll)
                                elif ll in lns_tomove_n12:
                                    pass
                                elif ll not in lns_tomove_n2:
                                    lns_tomove_n2.append(ll)
                                else:
                                    pass

            if entity in [-1, 2]:  # consider labels and holes in items
                # labels and holes can be moved right here
                # (no problem of coherence checks)
                if isinstance(item, label):
                    self.labels.remove(item)
                    xy = np.array([item.x, item.y])
                    new_xy = (rot_mat @ xy) - (rot_mat @ xy_c) + xy_c
                    labels.append(self.add_label(new_xy[0], new_xy[1],
                                                 item.ph_tag, item.area))
                elif isinstance(item, hole):
                    self.holes.remove(item)
                    xy = np.array([item.x, item.y])
                    new_xy = (rot_mat @ xy) - (rot_mat @ xy_c) + xy_c
                    holes.append(self.add_hole(new_xy[0], new_xy[1]))
            if not isinstance(item, (point, line, label, hole)):
                print('drawing.move_rotate() --> Skip element in list. Element is not a tiziano geometry entity.')

        # remove points that have to be moved with their associated lines
        self.remove_point(pts_tomove)
        # remove lines 
        for ll in lns_tomove_n1 + lns_tomove_n2 + lns_tomove_n12:
            if ll in self.lines:
                self.lines.remove(ll)
        # get number of points and lines now in the drawing --------------------
        # Counter variables to keep track of those objects that are newly
        # created and so they have to be returned in the end.
        # This check is applied only to points and lines because in the add_
        # methods of labels and holes no checks for existance in current drawing
        # are implemented yet.
        count_pts = len(self.points)
        count_lns = len(self.lines)
        # move points
        for pp in pts_tomove:
            xy = np.array([pp.x, pp.y])
            new_xy = (rot_mat @ xy) - (rot_mat @ xy_c) + xy_c
            pt, lns = self.add_point(new_xy[0], new_xy[1], pp.ph_tag,
                                     coherence=coherence, get_lns=True)
            # add pt to output points list only if pt is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(pt)
                count_pts = new_count_pts
            # add lns to output lines list
            lines.extend(lns)
            count_lns = len(self.lines)

        # move lines
        for ll in lns_tomove_n1:  # lines with only first endp to move
            xy1 = np.array([ll.n1.x, ll.n1.y])
            new_xy1 = (rot_mat @ xy1) - (rot_mat @ xy_c) + xy_c
            new_n1, lns = self.add_point(new_xy1[0], new_xy1[1], ll.n1.ph_tag,
                                         coherence=coherence, get_lns=True)
            # add new line endpoint to output points list
            # only if it is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(new_n1)
                count_pts = new_count_pts
            # add lns to output lines list
            lines.extend(lns)
            count_lns = len(self.lines)
            lns, pts = self.add_line(new_n1, ll.n2, ll.ph_tag,
                                     coherence=coherence, get_pts=True)
            # add lines to ouput list only if they are newly created
            new_count_lns = len(self.lines)
            if new_count_lns > count_lns:
                lines.extend(lns)
                count_lns = new_count_lns
            # add intersecting points to output list
            points.extend(pts)
            count_pts = len(self.points)
        for ll in lns_tomove_n2:  # lines with only second endp to move
            xy2 = np.array([ll.n2.x, ll.n2.y])
            new_xy2 = (rot_mat @ xy2) - (rot_mat @ xy_c) + xy_c
            new_n2, lns = self.add_point(new_xy2[0], new_xy2[1], ll.n2.ph_tag,
                                         coherence=coherence, get_lns=True)
            # add new line endpoint to output points list
            # only if it is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(new_n2)
                count_pts = new_count_pts
            # add lns to output lines list
            lines.extend(lns)
            count_lns = len(self.lines)
            lns, pts = self.add_line(ll.n1, new_n2, ll.ph_tag,
                                     coherence=coherence, get_pts=True)
            # add lines to ouput list only if they are newly created
            new_count_lns = len(self.lines)
            if new_count_lns > count_lns:
                lines.extend(lns)
                count_lns = new_count_lns
            # add intersecting points to output list
            points.extend(pts)
            count_pts = len(self.points)
        for ll in lns_tomove_n12:  # lines with both endp to move
            xy1 = np.array([ll.n1.x, ll.n1.y])
            new_xy1 = (rot_mat @ xy1) - (rot_mat @ xy_c) + xy_c
            xy2 = np.array([ll.n2.x, ll.n2.y])
            new_xy2 = (rot_mat @ xy2) - (rot_mat @ xy_c) + xy_c
            new_n1, lns1 = self.add_point(new_xy1[0], new_xy1[1], ll.n1.ph_tag,
                                          coherence=coherence, get_lns=True)
            # add new line endpoint to output points list
            # only if it is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(new_n1)
                count_pts = new_count_pts
            new_n2, lns2 = self.add_point(new_xy2[0], new_xy2[1], ll.n2.ph_tag,
                                          coherence=coherence, get_lns=True)
            # add new line endpoint to output points list
            # only if it is newly created
            new_count_pts = len(self.points)
            if new_count_pts > count_pts:
                points.append(new_n2)
                count_pts = new_count_pts
            # add splitted lines to output list
            lines.extend(lns1 + lns2)
            count_lns = len(self.lines)
            lns, pts = self.add_line(new_n1, new_n2, ll.ph_tag,
                                     coherence=coherence, get_pts=True)
            # add lines to ouput list only if they are newly created
            new_count_lns = len(self.lines)
            if new_count_lns > count_lns:
                lines.extend(lns)
                count_lns = new_count_lns
            # add intersecting points to output list
            points.extend(pts)
            count_pts = len(self.points)

        # final check not needed if coherence was not checked when new
        # points and lines are added to the drawing
        if coherence is False:
            return points, lines, labels, holes

        # FINAL COHERENCE CHECK
        # to return lines actually generated after move-transformation
        tol = 1e-10
        to_check = lines.copy()
        for ll in to_check:
            for p in points:
                d = ll.shortest_distance(p)
                # skip line endpoints (or points that are too close to them)
                if p.distance(ll.n1) < tol:
                    d = 2*tol
                if p.distance(ll.n2) < tol:
                    d = 2*tol
                if d < tol:
                    lines.remove(ll)
                    break

        return points, lines, labels, holes
        
    def mirror(self, items, x1, y1, x2, y2, entity=-1, coherence=True):
        """Generate symmetric objects with respect to a generic symmetry line.

        This function generates copies of tiziano geometry entities with respect
        to a generic symmetry line.
        Coherence can be checked or not (intersections with other points and
        lines already in the drawing can be investigated or skipped).
        If an element in 'items' is not a tiziano geom. entity, it is skipped
        without raising any errors.

        Parameters
        ----------
        items : list of tiziano geom entities [points, lines, labels, holes]
            List of tiziano geom. entities to mirror (i.e. a symmetric copy of
            each geom. element in 'items' will be created in the drawing).
        x1, y1, x2, y2 : int or float
            Coordinates of first and second points used to identify the symmetry
            line (the symmetry line passes through first and second point).
        entity : int, default -1
            Geom. entity type in 'items' to be mirrored.
            -1; all (default)
             0; only points
             1; only lines with their endpoints
             2; labels and holes
        coherence : bool, default True
            Coherence check (intersections with other points and lines).
            True  ==> check coherence (stronger but slower method) (default)
            False ==> do not check coherence (weaker but faster method)

        Returns
        -------
        points : list of tiziano points
            New points generated by mirror action.
        lines : list of tiziano lines
            New lines generated by mirror action.
        labels : list of tiziano labels
            New labels generated by mirror action.
        holes : list of tiziano holes
            New holes generated by mirror action.
        """
        # check items parameter input type
        if not isinstance(items, list):
            print("drawing.mirror() ==> Invalid argument. A list of tiziano geom. entities must be used.")
            return None, None, None, None
        # check input type for entity parameter
        if entity not in [-1, 0, 1, 2]:
            print('drawing.mirror() ==> Invalid argument! Entity must be -1, 0, 1, 2.')
            return None, None, None, None

        # compute transformation (symmetry) matrix
        if x1 != x2:  # symmetry line finite slope (non-vertical symmetry line)
            m = (y1 - y2)/(x1 - x2)
            q = y1 - m*x1
            den = 1 + m**2
            a = (1 - m**2)/den
            b = 2*m/den
            sym_mat = np.array([[a,  b],
                                [b, -a]])
            sym_vec = np.array([-2*m*q/den, 2*q/den])
        else:       # symmetry line with infinte slope (vertical symmetry line)
            sym_mat = np.array([[-1, 0], [0, 1]])
            sym_vec = np.array([2*x1, 0])

        # initialize lists of geom. entities to return--------------------------
        points = []
        lines  = []
        labels = []
        holes  = []
        # get number of points and lines already in the drawing-----------------
        # Counter variables are used inside the main loop to keep track of those
        # objects that are newly created and so they have to be returned.
        # This check is applied only to points and lines because in the add_
        # methods of labels and holes no checks for existance in current drawing
        # are implemented yet.
        count_pts = len(self.points)
        count_lns = len(self.lines)
        # loop over geom. entities to mirror -----------------------------------
        # In this main loop new coordinates of each geom. element to mirror are
        # computed using symmetry transformation. Then, add_ methods are called
        # and newly generated objects are appended to output points and lines
        # lists.
        for item in items:
            if entity in [-1, 0]:  # mirror points
                if isinstance(item, point):
                    xy = np.array([item.x, item.y])
                    new_xy = (sym_mat @ xy) + sym_vec
                    pt, lns = self.add_point(new_xy[0], new_xy[1], item.ph_tag,
                                             coherence=coherence,
                                             get_lns=True)
                    # add pt to output points list only if pt is newly created
                    new_count_pts = len(self.points)
                    if new_count_pts > count_pts:
                        points.append(pt)
                        count_pts = new_count_pts
                    # add splitted lines to output list
                    lines.extend(lns)
                    count_lns = len(self.lines)
            if entity in [-1, 1]:  # mirror lines
                if isinstance(item, line):
                    xy1 = np.array([item.n1.x, item.n1.y])
                    new_xy1 = (sym_mat @ xy1) + sym_vec
                    xy2 = np.array([item.n2.x, item.n2.y])
                    new_xy2 = (sym_mat @ xy2) + sym_vec
                    new_n1, lns1 = self.add_point(new_xy1[0], new_xy1[1],
                                                  item.n1.ph_tag,
                                                  coherence=coherence,
                                                  get_lns=True)
                    # add new line endpoints to output points list
                    # only if they are newly created
                    new_count_pts = len(self.points)
                    if new_count_pts > count_pts:
                        points.append(new_n1)
                        count_pts = new_count_pts
                    new_n2, lns2 = self.add_point(new_xy2[0], new_xy2[1],
                                                  item.n2.ph_tag,
                                                  coherence=coherence,
                                                  get_lns=True)
                    new_count_pts = len(self.points)
                    if new_count_pts > count_pts:
                        points.append(new_n2)
                        count_pts = new_count_pts
                    # add splitted lines to output list
                    lines.extend(lns1 + lns2)
                    count_lns = len(self.lines)
                    lns, pts = self.add_line(new_n1, new_n2, item.ph_tag,
                                             coherence=coherence,
                                             get_pts=True)
                    # add lns to output lines list only if ln is newly created
                    new_count_lns = len(self.lines)
                    if new_count_lns > count_lns:
                        lines.extend(lns)
                        count_lns = new_count_lns
                    # add intersecting points to output list
                    points.extend(pts)
                    count_pts = len(self.points)
            if entity in [-1, 2]:  # mirror labels and holes
                if isinstance(item, label):
                    xy = np.array([item.x, item.y])
                    new_xy = (sym_mat @ xy) + sym_vec
                    labels.append(self.add_label(new_xy[0], new_xy[1],
                                                 item.ph_tag, item.area))
                elif isinstance(item, hole):
                    xy = np.array([item.x, item.y])
                    new_xy = (sym_mat @ xy) + sym_vec
                    holes.append(self.add_hole(new_xy[0], new_xy[1]))
            if not isinstance(item, (point, line, label, hole)):
                print('drawing.mirror() --> Skip element in list. Element is not a tiziano geometry entity.')

        # final check not needed if coherence was not checked when new
        # points and lines are added to the drawing
        if coherence is False:
            return points, lines, labels, holes

        # FINAL COHERENCE CHECK
        # to return lines actually generated after copy-transformation
        tol = 1e-10
        to_check = lines.copy()
        for ll in to_check:
            for p in points:
                d = ll.shortest_distance(p)
                # skip line endpoints (or points that are too close to them)
                if p.distance(ll.n1) < tol:
                    d = 2*tol
                if p.distance(ll.n2) < tol:
                    d = 2*tol
                if d < tol:
                    lines.remove(ll)
                    break

        return points, lines, labels, holes

    def corner_arc_solver1(self, l1, l2, r):
        """Compute arc parameters to round a corner defined by two tiziano lines
        with a given radius.

        Parameters
        ----------
        l1, l2 : tiziano lines
            Tiziano lines that define a corner.
            These lines can or cannot be in the drawing.
        r : int or float
            Radius to round the given corner.

        Returns
        -------
        start, end : tiziano points
            Starting and ending points of rounding arc, respectivly.
        ip : tiziano point
            Corner point (common endpoint between the two given lines).
        arc_angle : float
            Rounding arc angle in rad.
        """
        if r == 0:
            print("drawing.corner_arc_solver1() ==> Null radius specified. No corner arc can be found.")
            return None, None, None, None

        # find common corner endpoint (intersection point ip)
        if l1.n1 == l2.n1:
            ip = l1.n1
            n1 = l1.n2
            n2 = l2.n2
        elif l1.n1 == l2.n2:
            ip = l1.n1
            n1 = l1.n2
            n2 = l2.n1
        elif l1.n2 == l2.n1:
            ip = l1.n2
            n1 = l1.n1
            n2 = l2.n2
        elif l1.n2 == l2.n2:
            ip = l1.n2
            n1 = l1.n1
            n2 = l2.n1
        else:
            print("drawing.corner_arc_solver1() ==> Input lines do not make a corner.")
            return None, None, None, None

        tol = 1e-10               # [u]   tolerance for linear quantities
        tol_ang = 0.01*np.pi/180  # [rad] tolerance for angles (0.01 deg)

        # define lines directing vectors
        u1 = np.array([n1.x - ip.x,
                       n1.y - ip.y])
        u2 = np.array([n2.x - ip.x,
                       n2.y - ip.y])

        # face numerical cancellation in directing vectors definition
        if abs(n1.x - ip.x) < tol:
            u1[0] = 0
        if abs(n1.y - ip.y) < tol:
            u1[1] = 0
        if abs(n2.x - ip.x) < tol:
            u2[0] = 0
        if abs(n2.y - ip.y) < tol:
            u2[1] = 0

        # compute norm of lines' directing vectors
        u1_n = np.linalg.norm(u1)
        u2_n = np.linalg.norm(u2)

        # lines' known term vector (common for both intersecting lines)
        p = np.array([ip.x,
                      ip.y])

        # compute angle between the two line using directing vectors
        theta = np.arccos((u1 @ u2)/(u1_n*u2_n))
        # if angle is (near) 180 deg, there is no corner to be rounded
        if np.pi-tol_ang < theta < np.pi+tol_ang:
            print("drawing.corner_arc_solver1() ==> Given corner cannot be rounded.")
            return None, None, None, None

        # compute distance between corner point and new points
        d = r/np.tan(theta/2)
        # if such distance is greater than one corner segment,
        if d > l1.manhattan_length() or d > l2.manhattan_length():
            # then corner cannot be rounded with the given radius
            print("drawing.corner_arc_solver1() ==> Radius too big to round the given corner.")
            return None, None, None, None

        # compute coordintes of new points using lines in parametric form
        nw_n1 = p + (d/u1_n)*u1
        nw_n2 = p + (d/u2_n)*u2
        # define tiziano point objects (but not added to the drawing yet!)
        new_n1 = point(nw_n1[0], nw_n1[1])
        new_n2 = point(nw_n2[0], nw_n2[1])

        # compute directing vectors angles
        # lines angle with positive x-axis in [0 - 2*pi] rad
        n1_angle = -abs(np.sign(u1[1]))*(-1 + np.sign(u1[1]))*np.pi + np.arctan2(u1[1], u1[0])
        n2_angle = -abs(np.sign(u2[1]))*(-1 + np.sign(u2[1]))*np.pi + np.arctan2(u2[1], u2[0])
        # find arc starting and ending points by looking at lines' angles
        # (in tiziano, arcs are defined positive counterclockwise)
        if n1_angle > n2_angle:
            if (n1_angle - n2_angle) < np.pi:
                start = new_n1
                end = new_n2
            else:
                start = new_n2
                end = new_n1
        else:
            if (n2_angle - n1_angle) < np.pi:
                start = new_n2
                end = new_n1
            else:
                start = new_n1
                end = new_n2
        # compute arc rounding angle
        arc_angle = np.pi - theta

        return start, end, ip, arc_angle

    def corner_arc_solver2(self, l, p, r):
        """Compute arc parameters to round a corner defined by a tiziano line and
        a tiziano point with a given radius.

        Parameters
        ----------
        l : tiziano line
            This line can or cannot be in the drawing.
        p : tiziano point
            This point can or cannot be in the drawing.
        r : int or float
            Radius to round the given corner.

        Returns
        -------
        start, end : tiziano points
            Starting and ending points of rounding arc, respectively.
        cp : tiziano point
            Arc center point.
        arc_angle : float
            Rounding arc angle in rad.
        """
        if r == 0:
            print("drawing.corner_arc_solver2() ==> Null radius specified. No corner arc can be found.")
            return None, None, None, None

        tol = 1e-10               # [u]   tolerance for linear quantities

        tol_rel = 1e-6            # relative tolerance for Newthon-Raphson
        max_iter = 20             # maximum iteration for Newthon-Raphson
        # starting point for Newthon-Raphson non-linear solver
        if line(p, l.n1).manhattan_length() > line(p, l.n2).manhattan_length():
            x0 = np.array([(p.x + l.n1.x)/2, (p.y + l.n1.y)/2])
        else:
            x0 = np.array([(p.x + l.n2.x)/2, (p.y + l.n2.y)/2])

        ux = l.n1.x - l.n2.x
        uy = l.n1.y - l.n2.y

        # if line is not vertical (ux =! 0 with a tolerance)
        if abs(ux) > tol:
            m = uy/ux
            if abs(uy) < tol:  # face numerical cancellation
                m = 0
            q = l.n1.y - m*l.n1.x
            F = lambda xc: np.array([
            (2*m*q - 2*xc[0] - 2*m*xc[1])**2 - 4*(1+m**2)*(q**2 - 2*q*xc[1] + xc[0]**2 + xc[1]**2 - r**2),
            p.x**2 + p.y**2 - 2*p.x*xc[0] - 2*p.y*xc[1] + xc[0]**2 + xc[1]**2 - r**2])
            J = lambda xc: np.array([
            [-4*(2*m*q - 2*xc[0] - 2*m*xc[1]) - 8*(1+m**2)*xc[0],         -2*p.x + 2*xc[0]],
            [-4*m*(2*m*q - 2*xc[0] - 2*m*xc[1]) - 8*(1+m**2)*(xc[1] - q), -2*p.y + 2*xc[1]]
            ])
            xc, _, ier, mesg = opt.fsolve(F, x0, fprime=J, full_output=True,
                                             col_deriv=True, xtol=tol_rel,
                                             maxfev=max_iter)
        else:
            # k param is computed as follow to inimize error when line is not
            # perfectly vertical
            k = (l.n1.x + l.n2.x)/2
            F = lambda xc: np.array([
            xc[0]**2 - 2*k*xc[0] - r**2 + k**2,
            p.x**2 + p.y**2 - 2*p.x*xc[0] - 2*p.y*xc[1] + xc[0]**2 + xc[1]**2 - r**2])
            J = lambda xc: np.array([
            [2*xc[0] - 2*k, -2*p.x + 2*xc[0]],
            [0            , -2*p.y + 2*xc[1]]
            ])
            xc, _, ier, mesg = opt.fsolve(F, x0, fprime=J, full_output=True,
                                             col_deriv=True, xtol=tol_rel,
                                             maxfev=max_iter)

        print("drawing.corner_arc_solver2() ==> ", mesg)
        if ier != 1:
            print("drawing.corner_arc_solver2() ==> Solution not found.")
            return None, None, xc, None

        # find second arc point (point in input line) by line intersection
        # - intersection between input line and its normal passing through xc -
        A = np.array([[ux,  uy],
                      [uy, -ux]])
        b = np.array([xc[0] - l.n1.x,
                      xc[1] - l.n1.y])
        tk = np.linalg.solve(A, b)

        # define tiziano points
        lp = point(l.n1.x + tk[0]*ux, l.n1.y + tk[0]*uy)  # line point
        if l.shortest_distance(lp) > tol:  # if line point is not on input line
            print("drawing.corner_arc_solver2() ==> Radius too big or too small for the given input line and point.")
            return None, None, None, None

        cp = point(xc[0], xc[1])                          # arc center point
        print("drawing.corner_arc_solver2() ==> Arc center ", cp)

        # define radius lines directing vectors
        u1 = np.array([p.x - cp.x,
                       p.y - cp.y])
        u2 = np.array([lp.x - cp.x,
                       lp.y - cp.y])

        # face numerical cancellation in directing vectors definition
        if abs(u1[0]) < tol:
            u1[0] = 0
        if abs(u1[1]) < tol:
            u1[1] = 0
        if abs(u2[0]) < tol:
            u2[0] = 0
        if abs(u2[1]) < tol:
            u2[1] = 0

        # compute norm of radius lines' directing vectors
        u1_n = np.linalg.norm(u1)
        u2_n = np.linalg.norm(u2)

        # compute arc angle between the two radius line using directing vectors
        arc_angle = np.arccos((u1 @ u2)/(u1_n*u2_n))

        # compute directing vectors angles
        # lines angle with positive x-axis in [0 - 2*pi] rad
        p_angle = -abs(np.sign(u1[1]))*(-1 + np.sign(u1[1]))*np.pi + np.arctan2(u1[1], u1[0])
        lp_angle = -abs(np.sign(u2[1]))*(-1 + np.sign(u2[1]))*np.pi + np.arctan2(u2[1], u2[0])
        # find arc starting and ending points by looking at lines' angles
        # (in tiziano, arcs are defined positive counterclockwise)
        if p_angle > lp_angle:
            if (p_angle - lp_angle) < np.pi:
                start = lp
                end = p
            else:
                start = p
                end = lp
        else:
            if (lp_angle - p_angle) < np.pi:
                start = p
                end = lp
            else:
                start = lp
                end = p

        return start, end, cp, arc_angle

    def round_corner1(self, l1, l2, r, step_angle, ph=-1, coherence=True):
        """Round corner between two lines with a given radius.

        Common corner point is removed from the drawing AFTER rounding arc
        creation.

        Parameters
        ----------
        l1, l2 : tiziano lines
            Tiziano lines that define a corner.
            These lines can or cannot be in the drawing.
        r : int or float
            Radius to round the given corner.
        step_angle : int or float
            Step_angle in radians to discretized the rounding arc.
        ph : int, defualt -1
            Optional physical tag for lines that discretized the rounding arc.
        coherence : bool, default True
            Check coherence or not when drawing rounding arc.

        Returns
        -------
        points : list of tiziano points
            Newly created points to add rounding arc.
        lines : list of tiziano lines
            Newly created lines while discretizing rounding arc.
        """
        start, end, ip, arc_angle = self.corner_arc_solver1(l1, l2, r)

        if start is None:  # returns empty lists when the corner cannot be rounded
            return [], []

        # add new endpoints to the drawing
        # (coherence must be always checked here to split corner lines)
        start, lns1 = self.add_point(start, get_lns=True)
        end, lns2 = self.add_point(end, get_lns=True)

        # add round arc to the drawing
        points, lines = self.add_arc(start, end, arc_angle, step_angle, ph=ph,
                                     coherence=coherence)
        # add to output list splitted corner lines
        lines.extend(lns1 + lns2)

        # remove corner point (and connected lines)
        self.remove_point([ip])

        # remove deleted lines from output list
        lines = [ll for ll in lines if not ll.is_end_point(ip)]

        return points, lines

    def round_corner2(self, l, p, r, step_angle, rp=None, ph=-1, coherence=True):
        """Round corner between a line and a point with a given radius.

        A specified line endpoint can be removed from the drawing AFTER rounding
        arc creation.

        Parameters
        ----------
        l : tiziano line
            This line can or cannot be in the drawing.
        p : tiziano point
            This point can or cannot be in the drawing.
        r : int or float
            Radius to round the given corner.
        step_angle : int or float
            Step_angle in radians to discretized the rounding arc.
        rp : tiziano point, default None
            Line endpoint to be removed from the drawing after arc creation.
            This point must be in the drawing.
        ph : int, defualt -1
            Optional physical tag for lines that discretized the rounding arc.
        coherence : bool, default True
            Check coherence or not when drawing rounding arc.

        Returns
        -------
        points : list of tiziano points
            Newly created points to add rounding arc.
        lines : list of tiziano lines
            Newly created lines while discretizing rounding arc.
        """
        start, end, _, arc_angle = self.corner_arc_solver2(l, p, r)

        if start is None:  # returns empty lists when solver failed
            return [], []

        # add new endpoints to the drawing
        # (coherence must be always checked here to split input line)
        start, lns1 = self.add_point(start, get_lns=True)
        end, lns2 = self.add_point(end, get_lns=True)

        # add round arc to the drawing
        points, lines = self.add_arc(start, end, arc_angle, step_angle, ph=ph,
                                     coherence=coherence)

        # add to output list splitted corner lines
        lines.extend(lns1 + lns2)

        # remove corner point (and connected lines)
        self.remove_point([rp])

        # remove deleted lines from output list
        lines = [ll for ll in lines if not ll.is_end_point(rp)]

        return points, lines

    def open_femm(self, filename, scale=None, Coherence=False):
        """
        Open a Femm file and load its geometry in the drawing.

        An optional scale parameter can be given.
        Arcs are discretized immediately.
        """
        if scale is None:
            scale = 1.0

        with open(filename, 'r') as fp:
            while True:
                ll = fp.readline()
                if not ll:
                    break
                linelist = ll.split()
                if linelist[0] == "[NumPoints]":
                    print('Importing '+ str(linelist[2]) +' points')
                    for _ in range(int(linelist[2])):
                        ll = fp.readline()
                        x, y = ll.split()[:2]
                        self.add_point(float(x)*scale, float(y)*scale, coherence=Coherence)
                if linelist[0] == "[NumSegments]":
                    print('Importing '+ str(linelist[2]) +' lines')
                    for _ in range(int(linelist[2])):
                        ll = fp.readline()
                        # n1, n2, MaxSideLength, ph, hidden, group
                        n1, n2, _, ph = ll.split()[:4]
                        if ph != 0:
                            self.add_line(self.points[int(n1)], self.points[int(n2)], int(ph), coherence=Coherence)
                        else:
                            self.add_line(self.points[int(n1)], self.points[int(n2)], coherence=Coherence)
                if linelist[0] == "[NumArcSegments]":
                    print('Importing '+ str(linelist[2]) +' arcs')
                    for _ in range(int(linelist[2])):
                        #print(str(_))
                        ll = fp.readline()
                        # n1, n2, ArcLength, MaxSideLength, ph, hidden, group
                        # femm use deg for angles for default
                        n1, n2, angle, step, ph = ll.split()[:5]
                        if ph != 0:
                            self.add_arc(self.points[int(n1)], self.points[int(n2)], float(angle)*np.pi/180, float(step)*np.pi/180, int(ph), coherence=Coherence)
                        else:
                            self.add_arc(self.points[int(n1)], self.points[int(n2)], float(angle)*np.pi/180, float(step)*np.pi/180, coherence=Coherence)
                if linelist[0] == "[NumBlockLabels]":
                    print('Importing '+ str(linelist[2]) +' labels')
                    for _ in range(int(linelist[2])):
                        ll = fp.readline()
                        x, y  = ll.split()[0:2]
                        area = ll.split()[3]
                        ph   = ll.split()[6]
                        self.add_label(float(x)*scale, float(y)*scale, int(ph), np.pi*((float(area)*scale)/2)**2)
                if linelist[0] == "[NumHoles]":
                    print('Importing '+ str(linelist[2]) +' holes')
                    for _ in range(int(linelist[2])):
                        ll = fp.readline()
                        x, y = ll.split()[:2]
                        self.add_hole(float(x)*scale, float(y)*scale)

    def save(self, filename):
        """Save the current tiziano drawing on a file."""
        now = datetime.datetime.now()

        with open(filename, 'w') as fp:
            fp.write("# dolomites-tiziano drawing\n")
            fp.write("# file generated: " + now.strftime("%Y-%m-%d_%H:%M") + '\n')
            # write points
            fp.write("points: %s\n" % (len(self.points)))
            for p in self.points:
                fp.write("%s %s %s\n" % (p.x, p.y, p.ph_tag))
            # write lines
            fp.write("lines: %s\n" % (len(self.lines)))
            for l in self.lines:
                n1 = self.points.index(l.n1)
                n2 = self.points.index(l.n2)
                fp.write("%s %s %s\n" % (n1, n2, l.ph_tag))
            # write labels
            fp.write("labels: %s\n" % (len(self.labels)))
            for ll in self.labels:
                fp.write("%s %s %s %s\n" % (ll.x, ll.y, ll.ph_tag, ll.area))
            # write holes
            fp.write("holes: %s\n" % (len(self.holes)))
            for h in self.holes:
                fp.write("%s %s\n" % (h.x, h.y))

    def load(self, filename, cl=True):
        """
        Open a tiziano draw previously saved.

        Note that coherence is not checked during this operation, the topology
        is assumed to be robust.
        """
        if cl is True:
            # clear all elements in the drawing, if any
            self.clear()

        # offset point index in the case of loading over a non empty drawing
        offset = len(self.points)
        # load a tiziano draw previously saved
        with open(filename, 'r') as fp:
            while True:
                _line = fp.readline()
                if not _line:
                    break
                linelist = _line.split()
                if linelist[0] == "#":
                    continue

                if linelist[0] == "points:":
                    for _ in range(int(linelist[1])):
                        _line    = fp.readline()
                        x, y, ph = _line.split()
                        _p = point(float(x), float(y), int(ph))
                        self.points.append(_p)
                if linelist[0] == "lines:":
                    for _ in range(int(linelist[1])):
                        _line     = fp.readline()
                        n1, n2, ph  = _line.split()
                        _l = line(self.points[int(n1)+offset], self.points[int(n2)+offset], int(ph))
                        self.lines.append(_l)
                if linelist[0] == "labels:":
                    for _ in range(int(linelist[1])):
                        _line     = fp.readline()
                        x, y, ph, area  = _line.split()
                        self.add_label(float(x), float(y), int(ph), float(area))
                if linelist[0] == "holes:":
                    for _ in range(int(linelist[1])):
                        _line     = fp.readline()
                        x, y  = _line.split()
                        self.add_hole(float(x), float(y))


    def get_points(self):
        """Export a list with points suitable to be used with pytriangle."""
        _list = []
        for p in self.points:
            _list.append((p.x, p.y))
        return _list

    def get_lines(self):
        """Export a list with lines suitable to be used with pytriangle (segments)."""
        _list = []
        _marks = []
        for l in self.lines:
            _list.append((self.points.index(l.n1), self.points.index(l.n2)))
            _marks.append(l.ph_tag)
        return _list, _marks

    def get_regions(self):
        """Export a list with regions suitable to be used with pytriangle."""
        _list = []
        for _, l in enumerate(self.labels):
            _list.append((float(l.x), float(l.y), float(l.ph_tag), float(l.area)))
        return _list

    def get_holes(self):
        """Export a list with holes suitable to be used with pytriangle."""
        _list = []
        for _, l in enumerate(self.holes):
            _list.append((l.x, l.y))
        return _list

    def mesh_triangle(self, flags='qpzaeA'):
        """Call triangle to mesh the drawing with the given flags."""
        points       = self.get_points()
        segs, marks  = self.get_lines()
        regs         = self.get_regions()
        hols         = self.get_holes()

        t = triangle.Triangle()
        t.set_points(points)
        t.set_segments(segs, marks)
        t.set_regions(regs)
        t.set_holes(hols)

        t.triangulate(mode=flags)

        self.triangle_points    = t.get_points()
        self.triangle_edges     = t.get_edges()
        self.triangle_triangles = t.get_triangles()
        return t

    def area(self, ph_list):
        """Compute area of regions with given physical tags.

        Areas of multiple separeted regions are summed together!.

        Parameters
        ----------
        ph_list : list of int
            List of physical tags of mesh triangles.
            Ph tags of mesh elements in the regions whose area want to be
            computed.

        Returns
        -------
        area : float
            Computed overall area of regions with given ph tags.
        """
        # check if triangle has been called previously
        if len(self.triangle_triangles) == 0:
            logging.error('Mesh has no triangle, please call mesh_triangle() first')
            return

        area = 0  # init area value as 0
        for tt in self.triangle_triangles:  # loop over mesh elements (triangles)
            if tt[2][0] in ph_list:         # if an element has a a given ph tag
                # define element node as tiziano points
                p1 = point(self.triangle_points[tt[0][0]][0][0],  # x1
                           self.triangle_points[tt[0][0]][0][1])  # y1
                p2 = point(self.triangle_points[tt[0][1]][0][0],  # x2
                           self.triangle_points[tt[0][1]][0][1])  # y2
                p3 = point(self.triangle_points[tt[0][2]][0][0],  # x3
                           self.triangle_points[tt[0][2]][0][1])  # y3
                # define arrays for two triangle edges
                v = np.array([p2.x - p1.x, p2.y - p1.y])
                u = np.array([p3.x - p1.x, p3.y - p1.y])
                V = line(p1, p2).manhattan_length()  # ||v|| length of v
                U = line(p1, p3).manhattan_length()  # ||u|| length of u
                # sum up triangle element area:
                # Area = (1/2)*||v||*||u||*sin(alpha),
                # with alpha angle between v and u
                # alpha = arccos((v @ u)/(||v||*||u||))
                # Since sin(arccos(x)) = sqrt(1 - x^2),
                # then: Area = (1/2)*||v||*||u||*sqrt(1 - ((v @ u)/(V*U))^2) 
                area += 0.5*V*U*np.sqrt(1 - ((v @ u)/(V*U))**2)

        return area

    def save_mesh(self, filename, _format='gmsh'):
        """
        Save the mesh in the given format filename.

        The output file _format: for the moment only gmsh 2.2 .
        """
        # Example -------------------------------------------------------------
        # $MeshFormat
        # 2.2 0 8
        # $EndMeshFormat
        # $Nodes
        # 5
        # 1 0 0 0
        # 2 0.001 0 0
        # 3 0.001 0.001 0
        # 4 -0 0.001 0
        # 5 0.0005 0.0005 0
        # $EndNodes
        # $Elements
        # 2
        # 1   1   2  10 181 320 319
        # 2   2   2  1000 1 1 2 3
        # $EndElements
        # ---------------------------------------------------------------------

        # check the correct mesh format
        if _format != 'gmsh':
            logging.error(('mesh format %s not supported') % _format)
            return

        # check if triangle has been called
        if len(self.triangle_triangles) == 0:
            logging.error('Mesh has no triangle, please call mesh_triangle() first')
            return

        # recover the mesh edges to be exported (i.e. with a physical tag)
        _list = []
        for e in self.triangle_edges:
            if e[1] > 0:
                # print (e)
                _list.append(e)

        with open(filename, 'w') as fp:

            # write out file heading
            fp.write('$MeshFormat\n')
            fp.write('2.2 0 8\n')
            fp.write('$EndMeshFormat\n')

            # write out nodes
            fp.write('$Nodes\n')
            fp.write(str(len(self.triangle_points)) + '\n')
            for i, p in enumerate(self.triangle_points):
                fp.write(('%s %s %s %s\n') % (i+1, p[0][0], p[0][1], 0))
            fp.write('$EndNodes\n')

            # write out elements
            fp.write('$Elements\n')
            fp.write(str(len(_list) + len(self.triangle_triangles)) + '\n')

            # write out liness
            for i, l in enumerate(_list):             # n 1 2     10    n        320         319
                fp.write(('%s 1 2 %s %s %s %s \n') % (i+1, int(l[1]), i+1, l[0][0]+1, l[0][1]+1))

            # write out triangles
            for i, t in enumerate(self.triangle_triangles):
                # ([n1, n2, n3], (), [regiona_tag])
                p = t[0]
                fp.write(('%s 2 2 %s %s %s %s %s\n') % (i+1+len(_list), int(t[2][0]), i+1, p[0]+1, p[1]+1, p[2]+1))
            fp.write('$EndElements\n')

    def save_dxf(self, filename):
        """Export a dxf copy of the drawing."""
        # if filename == '':
        #     logging.info('Please specify a valid filename to export a DXF of the drawing!')
        #     return None
        doc = ezdxf.new()
        msp = doc.modelspace()

        for l in self.lines:
            msp.add_line((l.n1.x, l.n1.y), (l.n2.x, l.n2.y), dxfattribs={'layer': 'Layer'})

        doc.saveas(filename)


# bool FemmppScene::writeGmshGeo(QString filename, bool exportSurface, qreal scale){
#
#
#     QDir Dir(qApp->applicationDirPath());
#
#     if (filename.isEmpty())
#         filename = QFileDialog::getSaveFileName(NULL,tr("Save geo file"), Dir.path(), tr("Gmsh file (*.geo)"));
#     if (filename.isEmpty()) {
#         printMessage(tr("No new file has been loaded"));
#         return false;}
#
#
#     QFile file(filename);
#     if (!file.open(QIODevice::WriteOnly | QIODevice::Text))   {
#         printMessage(tr("Unable to open file")+filename);
#         return false;
#     }
#     QTextStream out(&file);
#     out<<tr("// File generated by f2g plugin of Dolomites")<<endl;
#     out<<tr("// https://sourceforge.net/projects/dolomites/")<<endl<<endl;
#     out<<QString("scale = %1; // a global scale factor for the model").arg(scale)<<endl<<endl;
#
#
#     if (exportSurface) {
#         // This works only if each region is simply connected.
#         // Each region must have a different physical tag
#         QList<int> tags;  // a list to recover all the physical tags defined in the model
#
#         foreach (LabelItem * label, labellist) {
#
#             if (!tags.contains(label->physicalNumber)) tags.append(label->physicalNumber);
#         }
#
#
#         QMultiMap<int, int> physicalGroups;
#         QMultiMap<int, int> physicalLines;
#
#         int id_point = 1;
#         int id_line  = 1;
#         int id_loop  = 1;
#         int start_id_point;
#
#         for(int i =0; i< labellist.size(); ++i){
#
#
# //            QPolygonF region = mesh->getRegion(tag); // recover the union of triangles with the given tag
#             emit printMessage(tr("Writing surface %1 of %2").arg(i).arg(labellist.size()));
#             QPolygonF region = mesh->getRegion(i+1); // recover the union of triangles with the given tag
#
#             start_id_point = id_point - 1;
#
#     //        QList<int> points;
#             QList<int> lines;
#
#
#             for (int i=0; i<region.size()-1;++i) {// write out the list of nodes for the region
#                 out<<"Point("<< id_point <<") = {" << region[i].x() <<"*scale, "<<region[i].y()<<"*scale, 0};"<<endl;
#                 id_point++;
#             }
#
#             for (int i=0; i<region.size()-2;++i) {// write the boundary of the region
#                 out<<"Line("<< id_line <<") = {" << start_id_point+i+1 <<", " << start_id_point+i+2<<"};"<<endl;
#                 lines.append(id_line);
#                 id_line++;
#             }
#             out<<"Line("<< id_line <<") = {" << start_id_point+region.size()-1 <<", " << start_id_point+1<<"};"<<endl;
#             lines.append(id_line);
#             id_line++;
#
#             out<<"Line Loop ("<< id_loop<<") = {";
#             id_loop++;
#             for (int i=0; i<lines.size()-1;++i)
#                 out<<lines[i]<<",";
#             out<<lines[lines.size()-1]<<"};"<<endl;
#
#             out<<"Plane Surface("<<id_loop<<") = {"<<id_loop-1<<"};"<<endl<<endl;
#
#
#
#             // check to which physical group belong the current surface
#             int phy_id = labellist[i]->physicalNumber;
#             if (phy_id>0)
#                           physicalGroups.insert(phy_id,id_loop);
#
#
#         }
#
#         for(int i =0; i< linelist.size(); ++i){
#             int phy_id = linelist[i]->physicalNumber;
#             qDebug()<<"Physical Lines "<<phy_id;
#             if (phy_id>0){
#                 LineItem * line = linelist[i];
#                 out<<"Point("<< id_point <<") = {" << line->P1().x() <<"*scale, "<<line->P1().y()<<"*scale, 0};"<<endl;
#                 id_point++;
#                 out<<"Point("<< id_point <<") = {" << line->P2().x() <<"*scale, "<<line->P2().y()<<"*scale, 0};"<<endl;
#                 id_point++;
#                 out<<"Line("<< id_line <<") = {" << id_point-2 <<", " << id_point-1<<"};"<<endl;
#                 id_line++;
#                 physicalLines.insert(phy_id,id_line-1);
#             }
#         }
#
#         // check the physical tag on arcs
#         for(int i =0; i< arclist.size(); ++i){
#             ArcItem * arc = arclist[i];
#             if (arc->physicalNumber>0){
#                 QList<int> arcPoints;
#
#                 // add the first node of the arc
#                 out<<"Point("<< id_point <<") = {" << arc->N1->getActualPos().x() <<"*scale, "
#                                                    << arc->N1->getActualPos().y() <<"*scale, 0};"<<endl;
#                 id_point++;
#                 arcPoints.append(id_point-1);
#
#                 // recover the nodes in the middle from the arc and add them
#                 QList<QPointF> arcnodes = arclist[i]->discretize();
#                 if(arcnodes.size()>0){ // check if the arc has been discretize or not
#                     // add the arc nodes to the geo file
#                     foreach (QPointF P, arcnodes) {
#                         out<<"Point("<< id_point <<") = {" << P.x() <<"*scale, "<<P.y()<<"*scale, 0};"<<endl;
#                         id_point++;
#                         arcPoints.append(id_point-1);
#                     }
#                 }
#                // add the last node of the arc
#                 out<<"Point("<< id_point <<") = {" << arc->N2->getActualPos().x() <<"*scale, "
#                                                    << arc->N2->getActualPos().y()<<"*scale, 0};"<<endl;
#                 id_point++;
#                 arcPoints.append(id_point-1);
#
#                // add the arc segments
#                 for(int j=0; j<arcPoints.size()-2; ++j){
#                     out<<"Line("<< id_line <<") = {" << arcPoints[j]<<", " << arcPoints[j+1] <<"};"<<endl;
#                     id_line++;
#                     physicalLines.insert(arc->physicalNumber,id_line-1);
#                   }
#                 // add the last arc segments
#                 out<<"Line("<< id_line <<") = {" << arcPoints[arcPoints.size()-1]-1<<", " << id_point-1 <<"};"<<endl;
#                 id_line++;
#                 physicalLines.insert(arc->physicalNumber,id_line-1);
#
#           }
#
#         }
#
#
#
#
#
#         // print physical lines
#         QList<int> phy_list = physicalLines.uniqueKeys();
#
#         foreach (int ii, phy_list) {
#             QList<int> list = physicalLines.values(ii);
#             QString ids = "";
#             for (int i=0; i<list.size();++i){
#                 ids += QString("%1").arg(list[i]);
#                 if(i<list.size()-1) ids += ",";
#
#             }
#             out <<"Physical Line("<<ii<< ")= {" <<ids<<"};"<<endl;
#         }
#
#         // print physical sourface
#          phy_list = physicalGroups.uniqueKeys();
#
#         foreach (int ii, phy_list) {
#             QList<int> list = physicalGroups.values(ii);
#             QString ids = "";
#             for (int i=0; i<list.size();++i){
#                 ids += QString("%1").arg(list[i]);
#                 if(i<list.size()-1) ids += ",";
#
#             }
#             out <<"Physical Surface("<<ii<< ")= {" <<ids<<"};"<<endl;
#         }
#
#         out<<"Coherence;"<<endl;
#     }
#
#     // we just export points and line of the model.
#     // lineloops and surfaces must be then defined in gmsh
#     else{ // if not exportSurface
#
#         int id_point = 1;
#         int id_line  = 1;
#
#         foreach (LineItem * line, linelist){
#
#             out<<"Point("<< id_point <<") = {" << line->N1->actualPos.x() <<"*scale, "
#                                                << line->N1->actualPos.y() <<"*scale, 0};"<<endl;
#             id_point++;
#             out<<"Point("<< id_point <<") = {" << line->N2->actualPos.x() <<"*scale, "
#                                                << line->N2->actualPos.y() <<"*scale, 0};"<<endl;
#             id_point++;
#             out<<"Line("<< id_line <<") = {" << id_point-2 <<", " << id_point-1<<"};"<<endl;
#             id_line++;
#         }
#         foreach (ArcItem * arc, arclist){
#             out<<"Point("<< id_point <<") = {" << arc->N1->actualPos.x() <<"*scale, "
#                                                << arc->N1->actualPos.y() <<"*scale, 0};"<<endl;
#             id_point++;
#             out<<"Point("<< id_point <<") = {" << arc->center.x() <<"*scale, "
#                                                << arc->center.y() <<"*scale, 0};"<<endl;
#             id_point++;
#             out<<"Point("<< id_point <<") = {" << arc->N2->actualPos.x() <<"*scale, "
#                                                << arc->N2->actualPos.y() <<"*scale, 0};"<<endl;
#             id_point++;
#             out<<"Circle("<< id_line <<") = {" << id_point-3 <<", "
#                                                << id_point-2 <<", "
#                                                << id_point-1 <<"};"<<endl;
#             id_line++;
#         }
#
#         out<<"Coherence;"<<endl;
#
#     }
#
#
#     file.close();
#     return true;
# }
