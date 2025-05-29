# %% Geo
# Geo module improved capabilities:
# -) separated surfaces with the same ph tags are identified correctly
# -) lines with ph tags are identified correctly
# -) physical lines inside surfaces are treated properly

# %% Import useful modules
import numpy as np
from dolomites import tiziano

# %%
print('Module to generate geo files from (py)triangle mesh imported successfully.')

# %% Define functions to create a geo file from a (py)triangle mesh (dolomites/tiziano)
def loop_area(x, y):
    """
    Compute the area within a generic closed loop given points coordinates.
    
    The result is positive if the points are sorted counterclockwise, otherwise
    it is negative. The absolute value is the area within the loop.
    
    See https://www.wikihow.com/Calculate-the-Area-of-a-Polygon:
        --> Part 3 - Finding the area of irregular polygons.
    """
    np.append(x, x[0])
    np.append(y, y[0])
    return (np.sum(x[:-1]*y[1:]) - np.sum(y[:-1]*x[1:]))/2

def identify_loops(surf_bounds, mesh_nodes):
    """Identify curve loops of plane surfaces in Gmsh style."""
    # sort surfaces in ascending order according to number of bound edges
    # (outer surfaces with holes should be processed at the end)
    fsort = lambda sb: len(sb[1])  # sorting function
    surf_bounds.sort(key=fsort)
    curve_loops = []  # [ [[geo lines], [incidences]], ...]
    plane_surfs = []  # [ [ph_tag, [curve loops indeces]]]
    for sb in surf_bounds:  # loop over all the plane surfaces
        plane_surfs.append([sb[0], []])  # append ph tag to plane_surfaces list
        pls = plane_surfs[-1]
        edges = sb[1]         # <-- deepcopy here, maybe ??
        to_remove = []
        # if lines forming an already defined curve loop are founded inside boundary edges,  
        for loop_idx, loop in enumerate(curve_loops):
            if all([edg in edges for edg in loop[0]]):
                # then add the curve loop idx to plane
                pls[1].append(loop_idx)
                to_remove.extend(loop[0])
        for edg in to_remove:  # and remove those lines from boundary edges
            edges.remove(edg)
    
        while len(edges)>0:
            to_remove = []
            loop = []       # current loop edges (ordered)
            loop_inc = []   # current loop lines incidences (i.e. +1 or -1)
            x = []          # x-coordinates of loop nodes
            y = []          # y-coordinates of loop nodes
             # add the first element to current loop and remove it from surface edges
            next_edg = edges[0]
            loop.append(next_edg)
            to_remove.append(next_edg)
            begin = next_edg[0][0]  # the starting point idx of the current loop
            end = next_edg[0][1]    # end point idx of the current edge in the loop
            loop_inc.append(1)
            x.append(mesh_nodes[next_edg[0][0]][0][0])  # x begin
            x.append(mesh_nodes[next_edg[0][1]][0][0])  # x end (first)
            y.append(mesh_nodes[next_edg[0][0]][0][1])  # y begin
            y.append(mesh_nodes[next_edg[0][1]][0][1])  # y end (first)
    
            # build the loop (reorder edges)
            while end != begin:  # until the current end point idx is not equal to begin idx
                # find all the edges that have a point idx equal to current end point idx 
                jj = [ii for ii in edges if end in ii[0]]
                if len(jj) == 0:  # if no one is found, break the loop (critical situation!)
                    print('==> Something is wrong. Critical geometry encountered: unable to solve topology.')
                    break
                if len(jj) > 1:  # if there is more than one edges
                    try:
                        # take the first one that is different from the previous one
                        next_edg = next(filter(lambda ee: ee not in to_remove, jj))
                    except:
                        # take the first one that has an incidence opposite than before
                        next_edg = [kk for kk in jj if (kk[0][0]==end and loop_inc[loop.index(kk)]==-1) or (kk[0][1]==end and loop_inc[loop.index(kk)]==1)][0]
                else:  # else take the very first one as new edge (i.e. the previous one)
                    next_edg = jj[0]
                loop.append(next_edg)
                to_remove.append(next_edg)
                if next_edg[0][0] == end:  # and update current end point idx
                    loop_inc.append(1)
                    end = next_edg[0][1]
                else:
                    loop_inc.append(-1)
                    end = next_edg[0][0]
                x.append(mesh_nodes[end][0][0])
                y.append(mesh_nodes[end][0][1])
            
            area = loop_area(np.array(x), np.array(y))  # compute loop area
            # When the area within the loop results negative, the edges are sorted
            # clockwise. This leads to a normal with opposite direction, that may
            # causes troubles in a FE formulation.
            # Thus the loop direction is inverted by changing sign to edge incidences.
            if area < 0:
                loop_inc = ((-1)*np.array(loop_inc)).tolist()
                area = -area
            # append the built curve loop to curve_loops list
            curve_loops.append([loop, loop_inc, area])
            # and append its index to the last element of plane_surfaces list
            pls[1].append(len(curve_loops)-1)
            # remove all the edges in the just identified loop from surface boundaries
            to_remove = list(set(to_remove))  # get rid of double elements
            for ee in to_remove:
                edges.remove(ee)

        if len(pls[1]) > 1:  # is a plane surf is defined by more than one curve loop
            # put the curve loop that has the bigger area in first position
            loop_areas = [curve_loops[cl_idx][2] for cl_idx in pls[1]]
            max_area_idx = loop_areas.index(max(loop_areas))
            pls[1][0], pls[1][max_area_idx] = pls[1][max_area_idx], pls[1][0]

    return curve_loops, plane_surfs

def identify_surfbounds(triangles, edges):
    """Ïndentify boundary lines of surfaces."""
    # collect lines to be exported in geo file
    # ==> for the moment all mesh edges except the ones with ph tag = 0
    #     mesh edges with ph tag = 0 are for sure internal edges, i.e. not part
    #     of a surface boundary and without any physical meaning
    geo_lines = [edg for edg in edges if edg[1] != 0]
    
    mesh_elems = []  # preallocate 2D mesh elements list (it will be a list of lists)
    for tt in triangles:
        mesh_elem_edgs = []  # preallocation of edges of mesh element
        for edg in edges:
            if all([n_idx in tt[0] for n_idx in edg[0]]):
                mesh_elem_edgs.append(edg)
        # each mesh element has a ph tag and a list of three edges
        # [ph_tag, [((0, 1), -1), ((1, 2), 0), ((0, 2), 10)]]
        mesh_elems.append([tt[2][0], mesh_elem_edgs])

    ph_tags = []      # preallocation of ph tags list
    surf_bounds = []  # preallocation of surface boundaries list of lists
    # surf_bounds = [ [ph_tag, [boundary edges]], ...]
    for ee in mesh_elems:
        to_remove = []  # preallocation of suface boundaries to be removed
        ph_tag = ee[0]
        surf_bounds.append(ee)  # actual element is the last surface
        # if the ph tag of the current mesh element is not in the list yet,
        # then the element is separated surface for the moment
        if ph_tag not in ph_tags:
            ph_tags.append(ph_tag)
        else:  # else ...
            for sb in surf_bounds[:-1]:  # we have to check continuity
                if sb[0] == ph_tag:      # with all the other surfaces with the same ph tag
                    new_edgs = [ii for ii in sb[1] if ii not in surf_bounds[-1][1]]
                    comm_edgs = [ii for ii in sb[1] if ii in surf_bounds[-1][1] and ii[1]<=0]
                    # if there are some common edges (with ph tags 0 or -1 !)
                    if len(comm_edgs)>0:
                        # the last added surface incorporates the other 'óld'
                        # surface
                        surf_bounds[-1][1].extend(new_edgs)
                        # common edges are internal edges:
                        # they are removed from boundary list
                        for edg in comm_edgs:
                            surf_bounds[-1][1].remove(edg)
                            # and also from geo lines if they have no physical meaning
                            if edg[1] <= 0 and edg in geo_lines:
                                geo_lines.remove(edg)
                        to_remove.append(sb)
        for sb in to_remove:  # the old incorporated surfaces must be removed
            surf_bounds.remove(sb)

    return ph_tags, surf_bounds, geo_lines


def save_geo(draw, filename, lc, call_triangle=True, gmsh_mesh=None):
    """Write and save geo file from triangle mesh."""
    if call_triangle is True:  # call (py)triangle to generate a mesh
        # low quality flags --> minimum number of triangles
        draw.mesh_triangle(flags='pzeA')

    # check if triangle has been called previously
    if len(draw.triangle_triangles) == 0:
        print('Mesh has no triangle, please call mesh_triangle() first')
        return

    ph_tags, surf_bounds, geo_lns = identify_surfbounds(draw.triangle_triangles,
                                                          draw.triangle_edges)
    curve_loops, plane_surfs = identify_loops(surf_bounds, draw.triangle_points)
    
    with open(filename + '.geo', 'w') as fp:
        # write out file heading
        fp.write('// Geo file generated by dolomites/tiziano\n')
        fp.write('// --------------------------------------------------------\n')
        # local mesh size close to geo points
        fp.write(('lc = %s;\n')%(lc))

        fp.write('\n// Points and Lines ---------------------------------------\n')
        geo_pts = []
        ph_tags_lns = []
        for ll_idx, ll in enumerate(geo_lns):
            if ll[1] not in ph_tags_lns and ll[1] > 0:
                ph_tags_lns.append(ll[1])
            pt1 = draw.triangle_points[ll[0][0]]
            if pt1 not in geo_pts:
                geo_pts.append(pt1)
                fp.write(('Point(%s)={%s, %s, 0, lc};\n') % (int(ll[0][0]+1), pt1[0][0], pt1[0][1]))
            pt2 = draw.triangle_points[ll[0][1]]
            if pt2 not in geo_pts:
                geo_pts.append(pt2)
                fp.write(('Point(%s)={%s, %s, 0, lc};\n') % (int(ll[0][1]+1), pt2[0][0], pt2[0][1]))
            fp.write(('Line(%s)={%s, %s};\n') % (ll_idx+1, int(ll[0][0]+1), int(ll[0][1]+1)))

        fp.write('\n// Curve Loops --------------------------------------------\n ')
        for loop_idx, loop in enumerate(curve_loops):
            fp.write(('\nCurve Loop(%s)={') % (loop_idx+1))
            for ii, inc in enumerate(loop[1]):
                fp.write(('%s') % (int(inc*(geo_lns.index(loop[0][ii])+1))) )
                if ii < len(loop[1]) -1:
                    fp.write(', ')
            fp.write('};\n')

        fp.write('\n// Plane Surfaces -----------------------------------------\n')
        for s_idx, surf in enumerate(plane_surfs):
            fp.write(('Plane Surface(%s)={') % (s_idx+1))
            for jj, cl_idx in enumerate(surf[1]):
                fp.write(('%s') % (cl_idx+1))
                if jj < len(surf[1]) -1:
                    fp.write(', ')
            fp.write('};\n')

        fp.write('\n// Physical Surfaces --------------------------------------\n')
        for tag in ph_tags:
            ph_surf = [s_idx for s_idx, surf in enumerate(plane_surfs) if surf[0]==tag]
            fp.write(('Physical Surface(%s)={') % (int(tag)))
            for jj, s_idx in enumerate(ph_surf):
                fp.write(('%s') % (s_idx+1))
                if jj < len(ph_surf) - 1:
                    fp.write(', ')
            fp.write('};\n')

        fp.write('\n// Physical Lines -----------------------------------------\n')
        for tag in ph_tags_lns:
            ph_line = [l_idx for l_idx, ln in enumerate(geo_lns) if ln[1]==tag]
            fp.write(('Physical Line(%s)={') % (int(tag)))
            for ii, l_idx in enumerate(ph_line):
                fp.write(('%s') % (l_idx+1))
                if ii < len(ph_line) - 1:
                    fp.write(', ')
            fp.write('};\n')

        if gmsh_mesh is not None:
            fp.write('\n// Mesh -----------------------------------------------\n')
            fp.write('Mesh 2;\n')
            fp.write(('Mesh.MshFileVersion = %s;\n')%(gmsh_mesh))  # 2.0, 2.2, 4.1, ...
            fp.write('Save "' + filename + '.msh";\n')
