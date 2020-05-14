import numpy as np
import math
import trimesh
from stl import mesh
import timeit
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from scipy.signal import savgol_filter
from scipy.spatial import distance as distancelist
from scipy.signal import argrelextrema
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go
import sys, os
import plotly.offline
from PyQt5.QtCore import QUrl
#from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QApplication

equidistant_step_size = 0.5
percentile_pc = 5
max_points_in_pc = 10000

#Comment_DB: Based off of Tape_EA_Edit5 and stl_prepEdit2 (Starting point defined at starting edge) --> Goal: Automate bends


#Berechnet den Abstand von 2 Punkten
def distance(p1,p2):
    distance = np.linalg.norm(p2-p1)
    return distance

#Dreiecksnormalen(sortiert)
def tri_normals(ID_list,triangles,stl_normals):
    normals=[]
    #We generate our own normals with the id_list. Notice that the triangle order of stl_normals and the vertices
    #(triangles) is not the same

    for i in range(len(ID_list)):

        v1=triangles[int(ID_list[i])][0]-triangles[int(ID_list[i])][1]
        v2= triangles[int(ID_list[i])][0]-triangles[int(ID_list[i])][2]
        n=np.cross(v1,v2)
        n=(1 / np.linalg.norm(n)) * n
        normals.append(n)

    normals=np.asarray(normals)

    # the following average stl_normal always point at the outside of the object:
    avg_stl_normal = sum(stl_normals) / len(stl_normals)
    # average of the created normals:
    avg_sorted_normal=sum(normals) / len(normals)

    #Normalenvektoren werden immer in positive z-Richtung ausgerichtet
    if avg_sorted_normal[0]*avg_stl_normal[0]<0 or avg_sorted_normal[2]<0:
        normals=np.negative(normals)
    return normals

#Winkel zwischen den Dreiecksflächen
def normal_angles(tri_normals):
    normal_angles=[]
    for i in range(len(tri_normals)-1):
        n1=tri_normals[i]

        n2=tri_normals[i+1]

        dot=np.dot(n1,n2)

        if (1-dot)<0.0002:
            dot=1
        angle = np.arccos(dot)*180/math.pi
        normal_angles.append(angle)
    normal_angles=np.asarray(normal_angles)
    return normal_angles

#Dreiecksmittelpunkt
def tri_centerpoints (sorted_list):
    tri_centerpoints=[]
    for i in range(len(sorted_list)):
        center=np.array([(sorted_list[i][0][0]+sorted_list[i][1][0]+sorted_list[i][2][0])/3,(sorted_list[i][0][1]+sorted_list[i][1][1]+sorted_list[i][2][1])/3,(sorted_list[i][0][2]+sorted_list[i][1][2]+sorted_list[i][2][2])/3])
        tri_centerpoints.append(center)

    tri_centerpoints=np.asarray(tri_centerpoints)
    return tri_centerpoints

def tri_areas(patchvectors):
    #Berechnung der Dreiecksflächen und Speichern in einer Liste
    tri_surface = []
    for i in range(len(tri_centerpoints(patchvectors))):
        tri_surface.append(0.5 * (
            np.linalg.norm((patchvectors[i][0] - patchvectors[i][1]) - (patchvectors[i][0] - patchvectors[i][2]))))
    tri_surface = np.asarray(tri_surface)
    return tri_surface

#Punktwolke aus den Dreieckspunkten und Dreiecksmittelpunkten
def patch_pointcloud_weighted_by_area(patchvectors):
    #
    centerpoints_weights_area_tri = weights_center_points_by_percentil_area(patchvectors)

    # für jedes Dreieck werden die Mittelpunkte so oft gewertet (in die Punktwolke getan) wie es in centerpoints_weights_area_tri steht
    num_of_triangles = len(patchvectors)
    pointcloud_weighted=[]

    for i in range(num_of_triangles):
        for j in range(centerpoints_weights_area_tri[i]):
            pointcloud_weighted.append(tri_centerpoints(patchvectors)[i])

    pointcloud_weighted=np.asarray(pointcloud_weighted)

    return pointcloud_weighted


def weights_center_points_by_percentil_area(patchvectors):
    # In dieser Funktion wird eine Punktwolke aus den Dreiecksmittelpunkten der stl Datei generiert. Da mit dieser
    # Punktwolke später die Hauptwertzerlegung erfolgt, müssen große Dreiecke gegenüber kleinen Dreiecken stärker
    # gewichtet werden. Dies geschieht indem das Verhältnis der Dreiecksflächen zueinander berechnet wird und die
    # Mittelpunkte großer Dreiecke öfter zur Punktwolke hinzugefügt werden als die Mittelpunkte kleiner Dreiecke.
    ####Gewichtung großer Dreiecke
    # 1) Berechnung der Dreiecksflächen und Speichern in einer Liste
    tri_area_list = tri_areas(patchvectors)

    area_whole_patch = sum(tri_area_list)
    # 2) Die Dreiecksflächen werden miteinander verglichen, um später große Dreiecke gegenüber kleinen stärker zu
    # gewichten. Als "kleinste" Referenzfläche wird das 10-er Quantil der Flächen genommen (90% aller anderen Dreiecke
    # sind größer). In centerpoints_weights_area_tri wird für jedes Dreieck berechnet, um welchen Faktor es größer als das
    # Referenzdreicek ist (der Faktor wird aufgerundet). Zur Pointcloud wird von jedem Dreieck der Mittelpunkt so oft
    # hinzugefügt wie hoch der Flächenfaktor aus der centerpoints_weights_area_tri ist (mindestens einmal).
    # Um zu große Rechenzeiten zu vermeiden wird im Vorhinein abgeschätzt wie viele Punkte die Punktwolke
    # haben wird und bei Überschreiten eines Grenzwerts( max_points_in_pc) wird das Programm abgebrochen.
    #
    lower_percentil_area = np.percentile(tri_area_list, percentile_pc)
    estimated_number_points_in_pc = math.ceil(area_whole_patch / lower_percentil_area)

    # Abbruchbedingung: in der Pointcloud sollen maximal "max_points_in_pc" berechnet werden.
    if max_points_in_pc < estimated_number_points_in_pc:
        print("ERROR: Please use a .stl-object with reduced resolution ")
        print("Number of triangles: ", len(patchvectors))
        print("Estimated number of points in pointcloud:", estimated_number_points_in_pc)
        print("Allowed number of points in pointcloud:", max_points_in_pc)
        exit(1)

    # Im Folgenden wird jedes Dreieck mit dem kleinsten Dreieck verglichen und in centerpoints_weights_area_tri festgehalten, wie
    # oft das kleinste Dreieck in das jeweilige Dreieck hineinpasst.
    centerpoints_weights_area_tri = []
    for i in range(len(tri_area_list)):
        centerpoints_weights_area_tri.append(math.ceil(tri_area_list[i] / lower_percentil_area))

    return centerpoints_weights_area_tri


#Gerade durch Punktwolke mit kleinstem Abstand zu allen Punkten (COMMENT_DB: pc --> point cloud)
def pc_trendline(pointcloud):
    # ->Hauptkomponentenanalyse: https://de.wikipedia.org/wiki/Hauptkomponentenanalyse/Principal Component Analysis(PCA)
    # aka Singulärwertzerlegung / SVD
    # Generate some data that lies along a line
    #SOURCE: https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    #Explanation: https://www.tutorialspoint.com/scipy/scipy_linalg.htm
    maxvals = np.amax(pointcloud, axis=0)
    data = pointcloud
    # Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean = data.mean(axis=0)

    # Do an SVD(Singular Value Decomposition) on the mean-centered data.
    uu, dd, vv = np.linalg.svd(data - datamean)

    # Now vv[0] contains the first principal component, i.e. the direction
    # vector of the 'best fit' line in the least squares sense.

    # Now generate some points along this best fit line, for plotting.

    # It's a straight line, so we only need 2 points.
    linepts = vv[0] * np.mgrid[-max(maxvals):max(maxvals):2j][:, np.newaxis]

    # shift by the mean to get the line in the right place
    linepts += datamean #COMMENT_DB: linepts = linepts + datamean (shifting it in the same direction as the vv[0] direction vector)

    return linepts

#Hauptachsen der Punktwolke
def pc_axes(pointcloud):
    # Generate some data that lies along a line
    #SOURCE: https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    #Explanation: https://www.tutorialspoint.com/scipy/scipy_linalg.htm
    maxvals = np.amax(pointcloud, axis=0)

    data = pointcloud

    # Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean = data.mean(axis=0)

    # Do an SVD(Singular Value Decomposition) on the mean-centered data.
    uu, dd, vv = np.linalg.svd(data - datamean)

    axpts=[]
    for i in range(3):
        linepts = vv[i] * np.mgrid[-2*max(maxvals):2*max(maxvals):2j][:, np.newaxis]
        linepts += datamean
        axpts.append(linepts)

    axpts=np.asarray(axpts)

    return axpts

def pc_trendline_projection(pointcloud_weighted, triangle_centerpoints):
    # ->Hauptkomponentenanalyse: https://de.wikipedia.org/wiki/Hauptkomponentenanalyse/Principal Component Analysis(PCA)
    # aka Singulärwertzerlegung / SVD
    # Generate some data that lies along a line
    # SOURCE: https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    # Explanation: https://www.tutorialspoint.com/scipy/scipy_linalg.htm
    max_x_y_z_values = np.amax(pointcloud_weighted, axis=0)

    # Calculate the mean of the points, i.e. the 'center' of the cloud
    center_point_of_cloud = pointcloud_weighted.mean(axis=0)  #Mean of x,y,z-Values
    # Do an SVD on the mean-centered data.
    uu, dd, vv = np.linalg.svd(pointcloud_weighted - center_point_of_cloud)
    # Now vv[0] contains the first principal component, i.e. the direction

    # 2 points on line
    linepts = vv[0] * np.mgrid[-max(5*max_x_y_z_values):max(5*max_x_y_z_values):2j][:, np.newaxis]
    # shift by the mean to get the line in the right place
    linepts += center_point_of_cloud

    #SOURCE: https://gamedev.stackexchange.com/questions/72528/how-can-i-project-a-3d-point-onto-a-3d-line
    trendline_projection = [linepts[0]] #Comment_DKu_Wenzel: Warum linepts0 in der Liste?

    for i in range(len(triangle_centerpoints)):
        trendline_projection.append(project_pointtoline(triangle_centerpoints[i], linepts[0], linepts[1]))

    trendline_projection=np.asarray(trendline_projection)

    return trendline_projection

#Dreiecks-ID Reihenfolge entlang der Trendlinie
def sort_tri_id_by_trendline(trendline):
    # Sortiert die projizierten Punkte der Trendlinie entlang und gibt so eine Reihenfolge der
    # Dreiecke entlang der Trendlinienrichtung an. Es wird eine geordnete Liste der
    # Dreiecks-IDs zurückgegeben
    distance_list = []
    for i in range(1, len(trendline)):
        distance = np.array([np.linalg.norm((trendline[i] - trendline[0])), int(i - 1)])
        distance_list.append(distance)
    distance_list = np.asarray(distance_list)
    sorted_trendline_list = distance_list[distance_list[:, 0].argsort()]
    id_list=sorted_trendline_list[:,1]
    sorted_triangle_id_by_trendline=[]

    for i in range(len(id_list)):
        sorted_triangle_id_by_trendline.append(int(id_list[i]))

    return sorted_triangle_id_by_trendline

#Liste mit den sortierten auf die Trendline projizierten Punkte
def sorted_trendline_projection(unsorted_trendline, sorted_tri_ids):
    #Die auf die Trendline projizierten Mittelpunkte der Dreiecke werden entlang der Trendline sortiert
    sorted_trendline = []
    for i in range(len(sorted_tri_ids)):
        sorted_trendline.append(unsorted_trendline[sorted_tri_ids[i]+1])
    sorted_trendline=np.asarray(sorted_trendline)
    return sorted_trendline

#SAVITZKY-GOLAY-GLÄTTUNG
def smooth_savgol(equidistant_data_set,polynom_order,savgol_window_quotient):
    savgol_window = int(len(equidistant_data_set) / savgol_window_quotient)
    #polynom_order = 3 # Comment_DKu_Wenzel: Warum?
    # window must be impair:
    if savgol_window % 2 == 0:
        savgol_window += 1
    return savgol_filter(equidistant_data_set[:,1], savgol_window, polynom_order) # (data,window size, polynomial order)

#find the index of the closest value in an array to a given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def project_pointtoplane(Point_to_project,plane_normal,plane_point):
    P = Point_to_project
    S = plane_point
    n = plane_normal
    proj_point = (P - np.dot((P - S), n) / np.dot(n, n) * n)
    return proj_point

def project_pointtoline(Point_to_project,linept1,linept2):
    P = Point_to_project
    A = linept1
    B = linept2
    AB = B - A
    AP = P - A
    proj_point = A + np.dot(AP, AB) / np.dot(AB, AB) * AB
    # SOURCE: https://gamedev.stackexchange.com/questions/72528/how-can-i-project-a-3d-point-onto-a-3d-line

    return proj_point

####MAIN
def startparam(input_file,poly_order,savgol_window_quotient,max_distance):
    print(
        "Preprocessing läuft... Bei hochaufgelösten stl-Dateien und großen Flächengrößenunterschieden kann es zu "
        "längeren Berechnungszeiten kommen")

    #Read in the file:
    testpatch_vector = mesh.Mesh.from_file(input_file) #Comment_DB: stl mesh
    triangles = testpatch_vector.vectors #Comment_DB: triangle edges (wireframe)
    stl_normals = testpatch_vector.normals

    #Creating pointcloud:
    patch_pc_weighted = patch_pointcloud_weighted_by_area(triangles) #!!!!!Comment_DB: The blue points!

    #Creating trendline:
    trendline = pc_trendline_projection(patch_pc_weighted, tri_centerpoints(triangles))

    # Comment_DKu_Wenzel: 2 mal sorted
    #Sorted list of triangle points projected on trendline:
    sorted_projection_points = sorted_trendline_projection(trendline, sort_tri_id_by_trendline(trendline))
    #Sorted triangle normals:
    Normlist = tri_normals(sort_tri_id_by_trendline(trendline), triangles, stl_normals)

    #Average normal to determine the plane which we will use to create the 2D stripe:
    tri_surfaces = tri_areas(triangles)
    sorted_tri_surfaces = []
    for i in range(len(tri_surfaces)):
        sorted_tri_surfaces.append(tri_surfaces[sort_tri_id_by_trendline(trendline)[i]])
    total_surface=sum(sorted_tri_surfaces)
    weighted_norms = []
    for i in range(len(sorted_tri_surfaces)):
        weighted_norms.append(sorted_tri_surfaces[i]*Normlist[i]/total_surface)

    avg_tri_norm=sum(weighted_norms)

    # Singular value decomposition of the pointcloud to determine the main axis / directions of the patch
    datamean = patch_pc_weighted.mean(axis=0)
    uu, dd, vv = np.linalg.svd(patch_pc_weighted - datamean)

    # Definition der Hauptachsen, in denen die Winkel berechnet werden sollen
    trendline_x_axis = vv[0]
    trendline_x_axis = (1 / np.linalg.norm(trendline_x_axis)) * trendline_x_axis
    # avg_tri_norm ist nicht senkrecht zur x-Achse
    # von pcc + avg_tri_norm und zurück auf x-Achse projizieren
    trendline_avg_norm_point = datamean + np.dot(avg_tri_norm, trendline_x_axis) / np.dot(trendline_x_axis,
                                                                                          trendline_x_axis) \
                               * trendline_x_axis
    # y-Achse ist verbindung von pcc+avg_tri_norm mit dem projizierten Punkt
    trendline_y_axis = (datamean + avg_tri_norm) - trendline_avg_norm_point
    trendline_y_axis = (1 / np.linalg.norm(trendline_y_axis)) * trendline_y_axis
    trendline_z_axis = np.cross(trendline_x_axis, trendline_y_axis)

    # Dreiecksmittelpunkte im Sinne der Trendlinie entlang sortieren
    sorted_centerpoints = []

    for i in range(len(triangles)):
        sorted_centerpoints.append(tri_centerpoints(triangles)[sort_tri_id_by_trendline(trendline)[i]])
    sorted_centerpoints = np.asarray(sorted_centerpoints)

    # Start- und Endkante des Patches finden:
    # -> wenn bei großen Randdreiecken nur der Mittelpunkt benutzt wird, wird das Tape zu kurz. Daher werden
    # bei den Randdreiecken noch die Eckpunkte miteinbezogen
    startverts = triangles[sort_tri_id_by_trendline(trendline)[0]] #Comment_DB: startverts is start vertices
    endverts = triangles[sort_tri_id_by_trendline(trendline)[-1]]

    dist_startverts = []
    dist_endverts = []
    for i in range(3):
        dist_startverts.append(distance(startverts[i], datamean))
        dist_endverts.append(distance(endverts[i], datamean))
    startvert_3d = startverts[dist_startverts.index(max(dist_startverts))]
    endvert_3d = endverts[dist_endverts.index(max(dist_endverts))]

    # In Ebene projiziert:
    startvert_proj = project_pointtoplane(startvert_3d, trendline_z_axis, datamean)
    endvert_proj = project_pointtoplane(endvert_3d, trendline_z_axis, datamean)

    # auf Trendline projiziert
    startvert_trendline = project_pointtoline(startvert_3d, datamean + trendline_x_axis,
                                              datamean + 2 * trendline_x_axis)
    endvert_trendline = project_pointtoline(endvert_3d, datamean + trendline_x_axis,
                                            datamean + 2 * trendline_x_axis)

    # Für Start und Endpunkte müssen auch die jeweiligen Projektionen auf die Trendline berechnet werden, um den Abstand
    # zu berechnen (COMMENT_DB: insertion of points into sorted_projection_points list)

    sorted_projection_points = np.insert(sorted_projection_points, 0, startvert_trendline, axis=0)
    sorted_projection_points = np.concatenate((sorted_projection_points, [endvert_trendline]))

    # x-Werte: Abstand zwischen den sorted_projection_points
    x_list = []
    for i in range(len(sorted_projection_points)):
        x = np.linalg.norm((sorted_projection_points[0] - sorted_projection_points[i]))
        x_list.append(x)
    x_list = np.asarray(x_list)

    # Die Dreiecksmittelpunkte werden in x-y-Ebene der Trendline projiziert, um ein 2D Abstandsprofil zu erhalten
    tri_distance_xy_point = [startvert_proj]

    for i in range(len(sorted_centerpoints)):
        tri_distance_xy_point.append(project_pointtoplane(sorted_centerpoints[i], trendline_z_axis,
                                                          datamean))
    tri_distance_xy_point.append(endvert_proj)
    tri_distance_xy_point = np.asarray(tri_distance_xy_point)


    # xy-Distanzplot
    # berechnet den Abstand der auf xy-Ebene projizierten Dreiecksmittelpunkte zur Trendline

    tri_distance_xy = []
    for i in range(len(tri_distance_xy_point)):
        dist = np.linalg.norm(sorted_projection_points[i] - tri_distance_xy_point[i])

        # Vorzeichen ermitteln:
        if ((sorted_projection_points[i] + dist * trendline_y_axis)[0] - tri_distance_xy_point[i][0]) ** 2 < \
                ((sorted_projection_points[i] - dist * trendline_y_axis)[0] - tri_distance_xy_point[i][0]) ** 2:
            tri_distance_xy.append(dist)
        else:
            tri_distance_xy.append(-dist)
    tri_distance_xy = np.asarray(tri_distance_xy)

    # Funktion des xy-Abstands über die Länge der Trendline. Außerdem werden linear Punkte aufgefüllt, um eine
    # äquidistante Schrittgröße zu erreichen
    xy_patch_curve_step_size = equidistant_step_size
    xy_patch_curve = [] #Comment_DB: Smoothed blue curve points in array!
    for i in range(len(tri_distance_xy_point)):
        xy_patch_curve.append([x_list[i], tri_distance_xy[i]])

    for i in range(1, len(tri_distance_xy)):

        x_dist = xy_patch_curve[i][0] - xy_patch_curve[i - 1][0]
        if x_dist > xy_patch_curve_step_size:
            additional_steps = math.floor(x_dist / xy_patch_curve_step_size)

            for j in range(1, additional_steps + 1):

                xy_patch_curve.append([xy_patch_curve[i - 1][0] + j * xy_patch_curve_step_size, \
                                       (xy_patch_curve[i - 1][1] + (xy_patch_curve[i][1] - xy_patch_curve[i - 1][1]) \
                                        / (x_dist) * j * xy_patch_curve_step_size)])

    xy_patch_curve = np.asarray(xy_patch_curve)

    # Entlang der x-Werte sortieren
    xy_patch_curve = xy_patch_curve[xy_patch_curve[:, 0].argsort()]

    # Geglättete y-Werte mit SavitzkyGolay
    y_smooth = smooth_savgol(xy_patch_curve,poly_order,savgol_window_quotient)

    # 2D Knicklinie: Start - und Endpunkte;
    bend_pts_xy = []
    bend_pts_xy.append([xy_patch_curve[0][0], y_smooth[0]]) #Comment_DB: start point 2D (x coord, y coord)

    bend_pts_xy.append([xy_patch_curve[-1][0], y_smooth[-1]]) #Comment_DB: end point 2D (x coord, y coord)
    bend_pts_xy = np.asarray(bend_pts_xy)
    bend_pts_xy = bend_pts_xy[bend_pts_xy[:, 0].argsort()] #Comment_DB: sorted start and endpoints of 2D line that shows bends

    # Einfügen von weiteren Knickpunkten durch Finden von großen Abweichungen zur Kurve:
    set_max_divergence = max_distance #Comment_DB: User input
    insert_pts = True
    while insert_pts: #Comment_DB: this while loop goes all the way till "startparameter extrahieren"
        bend_pts_xy_curve = []
        bend_pts_xy_curve.append([bend_pts_xy[0][0], bend_pts_xy[0][1]]) #Comment_DB: only the first bend point (starting point at edge) appended to bend points curve list

        j = 1 #Comment_DB: at this point, bend_pts_xy curve only has the starting point in it, thus j = 1 is the number of points in the list. j = 1 is also the index of the NEXT point!

        for i in range(1, len(bend_pts_xy)): #Comment_DB: len(bend_pts_xy) is 2 for first iteration
            while bend_pts_xy_curve[-1][0] < bend_pts_xy[i][0]: #Comment_DB: while last x coord VALUE less than ith x coord VALUE in bend_pts_xy (If greater, then that means last point is reached)
                y_add = bend_pts_xy_curve[-1][1] + (bend_pts_xy[i - 1][1] - bend_pts_xy[i][1]) / \
                        (bend_pts_xy[i - 1][0] - bend_pts_xy[i][0]) * (xy_patch_curve[j][0] - xy_patch_curve[j - 1][0]) #Comment_DB: y = b + mx (finds next change in y linearly --> Produces a linear plot until end point at edge!!)
                bend_pts_xy_curve.append([xy_patch_curve[j][0], y_add]) #Comment_DB: append the NEXT point into the list
                j = j + 1 #Comment_DB: NEXT POINT

        bend_pts_xy_curve = np.asarray(bend_pts_xy_curve) #Comment_DB: This is now one linear curve from start to end point. Everything here is dependent on xy_patch_curve. Below will take divergence into consideration

        #Comment_DB: x,y coordinates of sav-gol curve in an array
        xy_savgol_curve = np.column_stack((xy_patch_curve[:, 0], y_smooth))

        #Comment_DB: curve_divergence in terms of y-distance # Größte Abweichung von geglätteter Kurve: (COMMENT_DB: By now all the points in the above (linear) line have been appended)
        curve_divergence_y = []
        for i in range(len(bend_pts_xy_curve)):
            curve_divergence_y.append([bend_pts_xy_curve[i][0], ((bend_pts_xy_curve[i][0]-xy_patch_curve[i][0])**2+(bend_pts_xy_curve[i][1]-y_smooth[i])**2)**0.5]) #Comment_DB: (x-coord vs. change in y-coord) take the x coord and y-distance between linear curve and sav-gol curve and append
        curve_divergence_y = np.asarray(curve_divergence_y)

        max_divergence = max([(v, i) for i, v in enumerate(curve_divergence_y[:, 1])]) #Comment_DB: returns distance, counter (Uses new curve_divergence)

        #Comment_DB: We know at which x-coord of bend_pts_xy_curve the max_divergence happens --> counter i

        bend_pts_xy = np.insert(bend_pts_xy, -1,
                                np.array([bend_pts_xy_curve[max_divergence[1]][0], y_smooth[max_divergence[1]]]), axis=0) #Comment_DB: insert a corner at x coord (counter i) and y coord (counter i) of max divergence
        bend_pts_xy = bend_pts_xy[bend_pts_xy[:, 0].argsort()] #Comment_DB: Bend points sorted in an array

        # no further points, if the chosen maximum distance is not surpassed
        if max_divergence[0] < set_max_divergence: #Comment_DB: This implies that there will be one extra bend, as the above code will have executed already
            insert_pts = False

        # Aktualisieren der Kurvenfunktion (COMMENT_DB: Same as above, except the corner is also in bend_pts_xy, so two linear lines (next iterations more and more linear lines..))
        bend_pts_xy_curve = []
        bend_pts_xy_curve.append([bend_pts_xy[0][0], bend_pts_xy[0][1]])

        j = 1

    ###### Startparameter extrahieren #####


    ###Start_r in 2D & 3D###
    pos_or_neg_adjust_xy = 1

    ###Start_r_atstart in 2D & 3D### COMMENT_DB: NEW DEFINITION
    Start_r_2d_atstart = bend_pts_xy[1] - bend_pts_xy[0]

    Start_r_3d_atstart = Start_r_2d_atstart[0] * trendline_x_axis + \
                 pos_or_neg_adjust_xy * Start_r_2d_atstart[1] * trendline_y_axis
    Start_r_3d_atstart = 1 / np.linalg.norm(Start_r_3d_atstart) * Start_r_3d_atstart

    ## Start_n_3d_atstart
    Start_n_3d_atstart = np.cross(trendline_z_axis, Start_r_3d_atstart)
    Start_n_3d_atstart = 1 / np.linalg.norm(Start_n_3d_atstart) * Start_n_3d_atstart

    # l-list for Chromo
    l_list = []
    for i in range(1, len(bend_pts_xy)):
        l = np.linalg.norm(bend_pts_xy[i] - bend_pts_xy[i - 1])
        l_list.append(l)
    l_list = np.asarray(l_list)

    # beta-list for Chromo
    beta_list = []
    for i in range(1, len(bend_pts_xy) - 1):
        r0 = bend_pts_xy[i] - bend_pts_xy[i - 1]
        r1 = bend_pts_xy[i + 1] - bend_pts_xy[i]
        angle = math.acos(np.dot(r0, r1) / (np.linalg.norm(r0) * np.linalg.norm(r1))) * 180 / math.pi
        steepness_r0 = r0[1] / r0[0]
        steepness_r1 = r1[1] / r1[0]
        if steepness_r1 < steepness_r0:
            angle = -angle
        if np.abs(angle) < 2: #Comment_DB: Small angles approx. to 0 degrees
            angle = 0
        beta_list.append(angle)
    beta_list = np.asarray(beta_list)
    #print("beta_list", beta_list) #Comment_DB: THIS IS THE SAME AS BETA_LIST FOR PREPROCESSED CHROMO IN LoP, EXCEPT HERE IT IS IN DEGREES

    # L_aim of savgol curve:
    L_aim = 0
    for i in range(1,len(xy_patch_curve)):
        p0=np.array([xy_patch_curve[i-1][0], y_smooth[i - 1]])
        p1=np.array([xy_patch_curve[i][0], y_smooth[i]])
        L_aim =L_aim + np.linalg.norm(p1-p0)

    start_parameter = [l_list, L_aim, beta_list, startvert_proj,
                       endvert_proj, Start_r_3d_atstart, Start_n_3d_atstart]
    return start_parameter

def show_startstrip(input_file,poly_order,savgol_window_quotient,max_distance,bestPatch_patternpoints,patch_start,patch_end):
    #Calculation:
    # Read in the file:
    testpatch_vector = mesh.Mesh.from_file(input_file)

    triangles = testpatch_vector.vectors
    stl_normals = testpatch_vector.normals

    #Creating pointcloud:
    patch_pc = patch_pointcloud_weighted_by_area(triangles)

    # Creating trendline:
    trendline = pc_trendline_projection(patch_pc, tri_centerpoints(triangles))

    # Sorted list of triangle points projected on trendline:
    sorted_projection_points = sorted_trendline_projection(trendline, sort_tri_id_by_trendline(trendline))

    # Sorted triangle normals:
    Normlist = tri_normals(sort_tri_id_by_trendline(trendline), triangles, stl_normals)

    # Average normal to determine the plane which we will use to create the 2D stripe:
    tri_surfaces = tri_areas(triangles)
    sorted_tri_surfaces = []
    for i in range(len(tri_surfaces)):
        sorted_tri_surfaces.append(tri_surfaces[sort_tri_id_by_trendline(trendline)[i]])
    total_surface = sum(sorted_tri_surfaces)
    weighted_norms = []
    for i in range(len(sorted_tri_surfaces)):
        weighted_norms.append(sorted_tri_surfaces[i] * Normlist[i] / total_surface)

    avg_tri_norm = sum(weighted_norms)

    # Singular value decomposition of the pointcloud to determine the main axis / directions of the patch
    datamean = patch_pc.mean(axis=0)
    uu, dd, vv = np.linalg.svd(patch_pc - datamean)

    # Definition der Hauptachsen, in denen die Winkel berechnet werden sollen
    trendline_x_axis = vv[0]
    trendline_x_axis = (1 / np.linalg.norm(trendline_x_axis)) * trendline_x_axis
    # avg_tri_norm ist nicht senkrecht zur x-Achse
    # von pcc + avg_tri_norm und zurück auf x-Achse projizieren
    trendline_avg_norm_point = datamean + np.dot(avg_tri_norm, trendline_x_axis) / np.dot(trendline_x_axis,
                                                                                          trendline_x_axis) \
                               * trendline_x_axis
    # y-Achse ist verbindung von pcc+avg_tri_norm mit dem projizierten Punkt
    trendline_y_axis = (datamean + avg_tri_norm) - trendline_avg_norm_point
    trendline_y_axis = (1 / np.linalg.norm(trendline_y_axis)) * trendline_y_axis
    trendline_z_axis = np.cross(trendline_x_axis, trendline_y_axis)

    # Dreiecksmittelpunkte im Sinne der Trendlinie entlang sortieren
    sorted_centerpoints = []

    for i in range(len(triangles)):
        sorted_centerpoints.append(tri_centerpoints(triangles)[sort_tri_id_by_trendline(trendline)[i]])
    sorted_centerpoints = np.asarray(sorted_centerpoints)

    # Start- und Endkante des Patches finden:
    # -> wenn bei großen Randdreiecken nur der Mittelpunkt benutzt wird, wird das Tape zu kurz. Daher werden
    # bei den Randdreiecken noch die Eckpunkte miteinbezogen
    startverts = triangles[sort_tri_id_by_trendline(trendline)[0]]
    endverts = triangles[sort_tri_id_by_trendline(trendline)[-1]]
    dist_startverts = []
    dist_endverts = []
    for i in range(3):
        dist_startverts.append(distance(startverts[i],datamean))
        dist_endverts.append(distance(endverts[i],datamean))
    startvert_3d = startverts[dist_startverts.index(max(dist_startverts))]
    endvert_3d = endverts[dist_endverts.index(max(dist_endverts))]

    # In Hauptebene projiziert:
    startvert_proj = project_pointtoplane(startvert_3d, trendline_z_axis, datamean)
    endvert_proj = project_pointtoplane(endvert_3d, trendline_z_axis, datamean)

    # auf Trendline projiziert
    startvert_trendline = project_pointtoline(startvert_3d, datamean + trendline_x_axis,
                                              datamean + 2 * trendline_x_axis)
    endvert_trendline = project_pointtoline(endvert_3d, datamean + trendline_x_axis,
                                            datamean + 2 * trendline_x_axis)

    #Für Start und Endpunkte müssen auch die jeweiligen Projektionen auf die Trendline berechnet werden, um den Abstand
    #zu berechnen

    sorted_projection_points = np.insert(sorted_projection_points, 0, startvert_trendline, axis=0)
    sorted_projection_points = np.concatenate((sorted_projection_points, [endvert_trendline]))

    # x-Werte: Abstand zwischen den sorted_projection_points
    x_list = []
    for i in range(len(sorted_projection_points)):
        x = np.linalg.norm((sorted_projection_points[0] - sorted_projection_points[i]))
        x_list.append(x)
    x_list = np.asarray(x_list)

    #Die Dreiecksmittelpunkte werden in x-y-Ebene der Trendline projiziert, um ein 2D Abstandsprofil zu erhalten
    tri_distance_xy_point = [startvert_proj]

    for i in range(len(sorted_centerpoints)):
        tri_distance_xy_point.append(project_pointtoplane(sorted_centerpoints[i], trendline_z_axis,
                                                          datamean))
    tri_distance_xy_point.append(endvert_proj)
    tri_distance_xy_point = np.asarray(tri_distance_xy_point)

    # xy-Distanzplot
    # berechnet den Abstand der auf xy-Ebene projizierten Dreiecksmittelpunkte zur Trendline
    tri_distance_xy = []
    for i in range(len(tri_distance_xy_point)):
        dist = np.linalg.norm(sorted_projection_points[i] - tri_distance_xy_point[i])

        # Vorzeichen ermitteln:
        if ((sorted_projection_points[i] + dist * trendline_y_axis)[0] - tri_distance_xy_point[i][0]) ** 2 < \
                ((sorted_projection_points[i] - dist * trendline_y_axis)[0] - tri_distance_xy_point[i][0]) ** 2:
            tri_distance_xy.append(dist)
        else:
            tri_distance_xy.append(-dist)
    tri_distance_xy = np.asarray(tri_distance_xy)

    # Funktion des xy-Abstands über die Länge der Trendline. Außerdem werden linear Punkte aufgefüllt, um eine
    # äquidistante Schrittgröße zu erreichen
    xy_patch_curve_step_size = equidistant_step_size
    xy_patch_curve = []
    for i in range(len(tri_distance_xy_point)):
        xy_patch_curve.append([x_list[i], tri_distance_xy[i]])

    #Die Punkte sicherheitshalber entlang der x_Achse sortieren
    xy_patch_curve=np.asarray(xy_patch_curve)
    xy_patch_curve = xy_patch_curve[xy_patch_curve[:, 0].argsort()]
    #zurück in Liste überführen um Punkte hinzuzufügen
    xy_patch_curve_raw_list = []
    for i in range(len(xy_patch_curve)):
        xy_patch_curve_raw_list.append([xy_patch_curve[i][0],xy_patch_curve[i][1]])


    for i in range(1, len(tri_distance_xy)):

        x_dist = xy_patch_curve_raw_list[i][0] - xy_patch_curve_raw_list[i - 1][0]
        if x_dist > xy_patch_curve_step_size:
            additional_steps = math.floor(x_dist / xy_patch_curve_step_size)


            for j in range(1, additional_steps + 1):


                xy_patch_curve_raw_list.append([xy_patch_curve_raw_list[i - 1][0] + j * xy_patch_curve_step_size, \
                                       (xy_patch_curve_raw_list[i - 1][1] + (xy_patch_curve_raw_list[i][1] - xy_patch_curve_raw_list[i - 1][1]) \
                                        / (x_dist) * j * xy_patch_curve_step_size)])

    xy_patch_curve = np.asarray(xy_patch_curve_raw_list)

    xy_patch_curve = xy_patch_curve[xy_patch_curve[:, 0].argsort()]

    y_smooth = smooth_savgol(xy_patch_curve, poly_order, savgol_window_quotient)

    # 2D Knicklinie: Start - und Endpunkte; lokale Extrema der geglätteten SavGol-Kurve
    bend_pts_xy = []
    bend_pts_xy.append([xy_patch_curve[0][0], y_smooth[0]])
    bend_pts_xy.append([xy_patch_curve[-1][0], y_smooth[-1]])
    bend_pts_xy = np.asarray(bend_pts_xy)
    bend_pts_xy = bend_pts_xy[bend_pts_xy[:, 0].argsort()]

    # Einfügen von weiteren Knickpunkten durch Finden von großen Abweichungen zur Kurve:
    set_max_divergence = max_distance
    insert_pts = True

    while insert_pts:
        bend_pts_xy_curve = []
        bend_pts_xy_curve.append([bend_pts_xy[0][0], bend_pts_xy[0][1]])
        j = 1
        for i in range(1, len(bend_pts_xy)):
            while bend_pts_xy_curve[-1][0] < bend_pts_xy[i][0]:
                y_add = bend_pts_xy_curve[-1][1] + (bend_pts_xy[i - 1][1] - bend_pts_xy[i][1]) / \
                        (bend_pts_xy[i - 1][0] - bend_pts_xy[i][0]) * (xy_patch_curve[j][0] - xy_patch_curve[j - 1][0])
                bend_pts_xy_curve.append([xy_patch_curve[j][0], y_add])
                j = j + 1

        bend_pts_xy_curve = np.asarray(bend_pts_xy_curve)

        # Größte Abweichung von geglätteter Kurve:
        curve_divergence = []
        for i in range(len(bend_pts_xy_curve)):
            curve_divergence.append([bend_pts_xy_curve[i][0], abs(bend_pts_xy_curve[i][1] - y_smooth[i])])
        curve_divergence = np.asarray(curve_divergence)
        max_divergence = (max([(v, i) for i, v in enumerate(curve_divergence[:, 1])]))

        #Comment_DB: x,y coordinates of sav-gol curve in an array

        bend_pts_xy = np.insert(bend_pts_xy, -1,np.array([bend_pts_xy_curve[max_divergence[1]][0], y_smooth[max_divergence[1]]]), axis=0)
        bend_pts_xy = bend_pts_xy[bend_pts_xy[:, 0].argsort()]
        # no further points, if the chosen maximum distance is not surpassed
        if max_divergence[0] < set_max_divergence:
            insert_pts = False

        # Aktualisieren der Kurvenfunktion
        bend_pts_xy_curve = []
        bend_pts_xy_curve.append([bend_pts_xy[0][0], bend_pts_xy[0][1]])

        j = 1
        for i in range(1, len(bend_pts_xy)):
            while bend_pts_xy_curve[-1][0] < bend_pts_xy[i][0]:
                y_add = bend_pts_xy_curve[-1][1] + (bend_pts_xy[i - 1][1] - bend_pts_xy[i][1]) / \
                        (bend_pts_xy[i - 1][0] - bend_pts_xy[i][0]) * (xy_patch_curve[j][0] - xy_patch_curve[j - 1][0])
                bend_pts_xy_curve.append([xy_patch_curve[j][0], y_add])
                j = j + 1

    ###2D-xy-PLOT
    plt.plot(xy_patch_curve[:, 0], xy_patch_curve[:, 1], 'bo', linewidth=2.0,label='ohne Glättung')  # äquidistante Punkte

    plt.plot(xy_patch_curve[:, 0], y_smooth, color='r', linewidth = 3, label ='Savitzky-Golay')  # SavGol-Glättung

    plt.plot(bend_pts_xy[:, 0], bend_pts_xy[:, 1], color='green', linewidth=3.0, label='lineare Angleichung')  # Streckenweise linear (nur Eckpunkte)
    plt.axis([x_list[0] - 50, x_list[-1] + 50, -1 * max(tri_distance_xy)-50, max(tri_distance_xy)+50])

    plt.xlabel('[ mm ]')
    plt.ylabel('[ mm ]')
    plt.legend()

    ###############3D-PLOTTING################

    figure = pyplot.figure() #Comment_DB: 3D plot of objective shape
    axes = mplot3d.Axes3D(figure)
    your_mesh = mesh.Mesh.from_file(input_file)
    patch_visual = mplot3d.art3d.Poly3DCollection(your_mesh.vectors, linewidths=1, alpha=0.5, edgecolor=[1, 1, 1], label = 'Geometriebereich') #Comment_DB: added edgecolor to make the edges visible

    axes.scatter([999999990],[9999999900],[9999999900],linewidths=0.0001, alpha = 0.5, label = "Geometriebereich") #Comment_DB: label in legend

    axes.scatter(datamean[0],datamean[1],datamean[2],c='g')

    # TRENDLINE
    axes.plot(*pc_axes(patch_pc)[0].T,label='Trendlinie', c='red') #Comment_DB: *pc_axes is *args, and .T is np.transpose
    axes.scatter(patch_start[0], patch_start[1], patch_start[2], c="black")
    axes.scatter(patch_end[0],patch_end[1],patch_end[2],c='black')

    # von PCC gemittelter Normalenvektor
    x3, y3, z3 = [datamean[0], datamean[0] + 500 * avg_tri_norm[0]], [datamean[1],
                                                                      datamean[1] + 500 * avg_tri_norm[1]], [
                     datamean[2], datamean[2] + 500 * avg_tri_norm[2]]
    plt.plot(x3,y3,z3,marker='o',c='green')

    patch_meshpoints = [] #Comment_DB: is not used

    for i in range(len(bestPatch_patternpoints) - 2):
        verts = [list(
            zip([bestPatch_patternpoints[i][0], bestPatch_patternpoints[i + 1][0], bestPatch_patternpoints[i + 2][0]], \
                [bestPatch_patternpoints[i][1], bestPatch_patternpoints[i + 1][1], bestPatch_patternpoints[i + 2][1]], \
                [bestPatch_patternpoints[i][2], bestPatch_patternpoints[i + 1][2], bestPatch_patternpoints[i + 2][2]]))] #Comment_DB: DARK BLUE LoP PATCH
        axes.add_collection3d(Poly3DCollection(verts), zs='z') #Comment_DB: INSERT LoP PATCH IN GRAPH
        patch_meshpoints.append(verts) #Comment_DB: is not used
    axes.scatter(bestPatch_patternpoints[:, 0], bestPatch_patternpoints[:, 1], bestPatch_patternpoints[:, 2], c='r')

    face_color = [0.5, 0.5, 1]  # alternative: matplotlib.colors.rgb2hex([0.5, 0.5, 1])
    patch_visual.set_facecolor(face_color)
    axes.legend()
    axes.add_collection3d(patch_visual) #Comment_DB: stl mesh file

    axes.autoscale(enable=False, axis='both')  # you will need this line to change the Z-axis
    axes.set_xbound(-150, 150)
    axes.set_ybound(-50, 250)
    axes.set_zbound(-150, 150)

    pyplot.axis('off')
    pyplot.show(figure)
    return

if __name__ == "__main__":
    show_startstrip(input_file)
