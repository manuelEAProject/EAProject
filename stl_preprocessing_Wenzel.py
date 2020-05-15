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
def tri_normals(id_list, triangle_vectors_of_stl, stl_normals):
    normals=[]
    #We generate our own normals with the id_list. Notice that the triangle order of stl_normals and the vertices
    #(triangles) is not the same

    for i in range(len(id_list)):

        v1= triangle_vectors_of_stl[int(id_list[i])][0] - triangle_vectors_of_stl[int(id_list[i])][1]
        v2= triangle_vectors_of_stl[int(id_list[i])][0] - triangle_vectors_of_stl[int(id_list[i])][2]
        n=np.cross(v1,v2)
        n=(1 / np.linalg.norm(n)) * n
        normals.append(n)

    normals=np.asarray(normals)


    # Normalenvektoren werden immer in positive z-Richtung ausgerichtet

    # the following average stl_normal always point at the outside of the object:
    avg_stl_normal = sum(stl_normals) / len(stl_normals)
    # average of the created normals:
    avg_sorted_normal = sum(normals) / len(normals)

    true_when_stl_and_tri_normal_not_same_direction = avg_sorted_normal[0] * avg_stl_normal[0] < 0
    true_when_z_from_tri_normal_neg = avg_sorted_normal[2] < 0

    if true_when_stl_and_tri_normal_not_same_direction or true_when_z_from_tri_normal_neg:
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
def tri_centerpoints():
    tri_centerpoints=[]
    for i in range(len(triangle_vectors_of_stl)):
        center=np.array([(triangle_vectors_of_stl[i][0][0] + triangle_vectors_of_stl[i][1][0] + triangle_vectors_of_stl[i][2][0]) / 3, (triangle_vectors_of_stl[i][0][1] + triangle_vectors_of_stl[i][1][1] + triangle_vectors_of_stl[i][2][1]) / 3, (triangle_vectors_of_stl[i][0][2] + triangle_vectors_of_stl[i][1][2] + triangle_vectors_of_stl[i][2][2]) / 3])
        tri_centerpoints.append(center)

    tri_centerpoints=np.asarray(tri_centerpoints)
    return tri_centerpoints

def tri_areas():
    #Berechnung der Dreiecksflächen und Speichern in einer Liste
    tri_surface = []
    for i in range(len(tri_centerpoints)):
        tri_surface.append(0.5 * (
            np.linalg.norm((triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][1]) - (triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][2]))))
    tri_surface = np.asarray(tri_surface)
    return tri_surface

#Punktwolke aus den Dreieckspunkten und Dreiecksmittelpunkten
def patch_pointcloud_weighted_by_area():
    #
    centerpoints_weights_area_tri = weights_center_points_by_percentil_area()

    # für jedes Dreieck werden die Mittelpunkte so oft gewertet (in die Punktwolke getan) wie es in centerpoints_weights_area_tri steht
    num_of_triangles = len(triangle_vectors_of_stl)
    pointcloud_weighted=[]

    for i in range(num_of_triangles):
        for j in range(centerpoints_weights_area_tri[i]):
            pointcloud_weighted.append(tri_centerpoints[i])

    pointcloud_weighted=np.asarray(pointcloud_weighted)

    return pointcloud_weighted

def weights_center_points_by_percentil_area():
    # In dieser Funktion wird eine Punktwolke aus den Dreiecksmittelpunkten der stl Datei generiert. Da mit dieser
    # Punktwolke später die Hauptwertzerlegung erfolgt, müssen große Dreiecke gegenüber kleinen Dreiecken stärker
    # gewichtet werden. Dies geschieht indem das Verhältnis der Dreiecksflächen zueinander berechnet wird und die
    # Mittelpunkte großer Dreiecke öfter zur Punktwolke hinzugefügt werden als die Mittelpunkte kleiner Dreiecke.
    ####Gewichtung großer Dreiecke
    # 1) Berechnung der Dreiecksflächen und Speichern in einer Liste


    area_whole_patch = sum(tri_areas)
    # 2) Die Dreiecksflächen werden miteinander verglichen, um später große Dreiecke gegenüber kleinen stärker zu
    # gewichten. Als "kleinste" Referenzfläche wird das 10-er Quantil der Flächen genommen (90% aller anderen Dreiecke
    # sind größer). In centerpoints_weights_area_tri wird für jedes Dreieck berechnet, um welchen Faktor es größer als das
    # Referenzdreicek ist (der Faktor wird aufgerundet). Zur Pointcloud wird von jedem Dreieck der Mittelpunkt so oft
    # hinzugefügt wie hoch der Flächenfaktor aus der centerpoints_weights_area_tri ist (mindestens einmal).
    # Um zu große Rechenzeiten zu vermeiden wird im Vorhinein abgeschätzt wie viele Punkte die Punktwolke
    # haben wird und bei Überschreiten eines Grenzwerts( max_points_in_pc) wird das Programm abgebrochen.
    #
    lower_percentil_area = np.percentile(tri_areas, percentile_pc)
    estimated_number_points_in_pc = math.ceil(area_whole_patch / lower_percentil_area)

    # Abbruchbedingung: in der Pointcloud sollen maximal "max_points_in_pc" berechnet werden.
    if max_points_in_pc < estimated_number_points_in_pc:
        print("ERROR: Please use a .stl-object with reduced resolution ")
        print("Number of triangles: ", len(triangle_vectors_of_stl))
        print("Estimated number of points in pointcloud:", estimated_number_points_in_pc)
        print("Allowed number of points in pointcloud:", max_points_in_pc)
        exit(1)

    # Im Folgenden wird jedes Dreieck mit dem kleinsten Dreieck verglichen und in centerpoints_weights_area_tri festgehalten, wie
    # oft das kleinste Dreieck in das jeweilige Dreieck hineinpasst.
    centerpoints_weights_area_tri = []
    for i in range(len(tri_areas)):
        centerpoints_weights_area_tri.append(math.ceil(tri_areas[i] / lower_percentil_area))

    return centerpoints_weights_area_tri

#Gerade durch Punktwolke mit kleinstem Abstand zu allen Punkten (COMMENT_DB: pc --> point cloud)
def pc_trendline():
    # ->Hauptkomponentenanalyse: https://de.wikipedia.org/wiki/Hauptkomponentenanalyse/Principal Component Analysis(PCA)
    # aka Singulärwertzerlegung / SVD
    # Generate some data that lies along a line
    #SOURCE: https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    #Explanation: https://www.tutorialspoint.com/scipy/scipy_linalg.htm
    maxvals = np.amax(patch_pc_weighted, axis=0)

    # It's a straight line, so we only need 2 points.
    linepts = first_principal_components_pc_weighted[2][0]* np.mgrid[-2*max(maxvals):2*max(maxvals):2j][:, np.newaxis]

    # shift by the mean to get the line in the right place
    linepts += center_point_of_cloud_weighted #COMMENT_DB: linepts = linepts + datamean (shifting it in the same direction as the vv[0] direction vector)

    return linepts

def pc_trendline_projection():
    #SOURCE: https://gamedev.stackexchange.com/questions/72528/how-can-i-project-a-3d-point-onto-a-3d-line
    trendline_projection = [trendline[0]] #Comment_DKu_Wenzel: Warum linepts0 in der Liste?

    for i in range(len(tri_centerpoints)):
        trendline_projection.append(project_pointtoline(tri_centerpoints[i], trendline[0], trendline[1]))

    trendline_projection=np.asarray(trendline_projection)

    return trendline_projection

def first_principal_components_pc_weighted():
    # ->Hauptkomponentenanalyse: https://de.wikipedia.org/wiki/Hauptkomponentenanalyse/Principal Component Analysis(PCA)
    # aka Singulärwertzerlegung / SVD
    # SOURCE: https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    # Explanation: https://www.tutorialspoint.com/scipy/scipy_linalg.htm

    # Do an SVD on the mean-centered data.
    first_principal_components_pc_weighted = np.linalg.svd(patch_pc_weighted - center_point_of_cloud_weighted)
    # Now the first principal component contains [uu, dd, vv] , where vv[0] is the direction
    return first_principal_components_pc_weighted

#Dreiecks-ID Reihenfolge entlang der Trendlinie
def sort_tri_id_by_trendline():
    # Sortiert die projizierten Punkte der Trendlinie entlang und gibt so eine Reihenfolge der
    # Dreiecke entlang der Trendlinienrichtung an. Es wird eine geordnete Liste der
    # Dreiecks-IDs zurückgegeben
    distance_startpoint_trendline_to_projectionpoint = []
    for i in range(1, len(trendline_projection_points)):
        distance = np.array([np.linalg.norm((trendline_projection_points[i] - trendline_projection_points[0])), i - 1])
        distance_startpoint_trendline_to_projectionpoint.append(distance)

    # Put distances in a array for sorting triangles by size and keeping the global id of triangles
    distance_startpoint_trendline_to_projectionpoint = np.asarray(distance_startpoint_trendline_to_projectionpoint)

    sorted_distance_by_size = distance_startpoint_trendline_to_projectionpoint[distance_startpoint_trendline_to_projectionpoint[:, 0].argsort()]
    # Get id and put them in a list (from array
    id_list = list(map(int, sorted_distance_by_size[:,1]))

    return id_list

#Liste mit den sortierten auf die Trendline projizierten Punkte
def sort_trendline_projection(trendline_projection_points, sorted_tri_ids):
    #Die auf die Trendline projizierten Mittelpunkte der Dreiecke werden entlang der Trendline sortiert
    sorted_trendline = []
    for i in range(len(sorted_tri_ids)):
        sorted_trendline.append(trendline_projection_points[sorted_tri_ids[i] + 1])
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
    # Comment_DKu_Wenzel: global defined so they dont have to be calculated multiple times. Stay always the same for same stl file

    #Read in the file:
    patch_vectors_of_stl_input = mesh.Mesh.from_file(input_file) #Comment_DB: stl mesh
    global triangle_vectors_of_stl
    triangle_vectors_of_stl = patch_vectors_of_stl_input.vectors #Comment_DB: triangle edges (wireframe)

    global stl_normals
    stl_normals = patch_vectors_of_stl_input.normals

    global tri_centerpoints # Comment_DKu_Wenzel: basically the unweighted point cloud
    tri_centerpoints= tri_centerpoints()

    global tri_areas
    tri_areas = tri_areas()

    #Creating pointcloud:
    global patch_pc_weighted
    patch_pc_weighted = patch_pointcloud_weighted_by_area() #!!!!!Comment_DB: The blue points!
    # Calculate the mean of the points, i.e. the 'center' of the cloud
    global center_point_of_cloud_weighted
    center_point_of_cloud_weighted = patch_pc_weighted.mean(axis=0)  # Mean of x,y,z-Values
    # Do Principal Component Analysis(PCA) on the mean-centered data. AKA SVD
    # The first principal component contains [uu, dd, vv] , where vv[0] is the direction
    global first_principal_components_pc_weighted
    first_principal_components_pc_weighted = np.linalg.svd(patch_pc_weighted - center_point_of_cloud_weighted)

    # Creating trendline
    global trendline
    trendline = pc_trendline()

    # Creating trendline_projection_points:
    global trendline_projection_points
    trendline_projection_points = pc_trendline_projection()

    global tri_id_sorted
    tri_id_sorted = sort_tri_id_by_trendline()
    #Sorted list of triangle points with trendline_projection_points:
    global sorted_projection_points
    sorted_projection_points = sort_trendline_projection(trendline_projection_points, tri_id_sorted)

    #Sorted triangle normals:
    sorted_triangle_normals = tri_normals(tri_id_sorted, triangle_vectors_of_stl, stl_normals)

    global avg_tri_norm
    #Average normal to determine the plane which we will use to create the 2D stripe:
    sorted_tri_surfaces = []
    for i in range(len(tri_areas)):
        sorted_tri_surfaces.append(tri_areas[tri_id_sorted[i]])

    weighted_norms = []
    for i in range(len(sorted_tri_surfaces)):
        weighted_norms.append( (sorted_tri_surfaces[i]/sum(tri_areas)) * sorted_triangle_normals[i])

    avg_tri_norm=sum(weighted_norms)


    # Singular value decomposition of the pointcloud to determine the main axis / directions of the patch

    uu, dd, vv = np.linalg.svd(patch_pc_weighted - center_point_of_cloud_weighted)

    # Definition der Hauptachsen, in denen die Winkel berechnet werden sollen
    trendline_x_axis = vv[0]
    trendline_x_axis = (1 / np.linalg.norm(trendline_x_axis)) * trendline_x_axis
    # avg_tri_norm ist nicht senkrecht zur x-Achse
    # von pcc + avg_tri_norm und zurück auf x-Achse projizieren
    trendline_avg_norm_point = center_point_of_cloud_weighted + np.dot(avg_tri_norm, trendline_x_axis) / np.dot(trendline_x_axis,
                                                                                                                trendline_x_axis) \
                               * trendline_x_axis
    # y-Achse ist verbindung von pcc+avg_tri_norm mit dem projizierten Punkt
    trendline_y_axis = (center_point_of_cloud_weighted + avg_tri_norm) - trendline_avg_norm_point
    trendline_y_axis = (1 / np.linalg.norm(trendline_y_axis)) * trendline_y_axis
    trendline_z_axis = np.cross(trendline_x_axis, trendline_y_axis)

    # Dreiecksmittelpunkte im Sinne der Trendlinie entlang sortieren
    sorted_centerpoints = []

    for i in range(len(triangle_vectors_of_stl)):
        sorted_centerpoints.append(tri_centerpoints[tri_id_sorted[i]])
    sorted_centerpoints = np.asarray(sorted_centerpoints)

    # Start- und Endkante des Patches finden:
    # -> wenn bei großen Randdreiecken nur der Mittelpunkt benutzt wird, wird das Tape zu kurz. Daher werden
    # bei den Randdreiecken noch die Eckpunkte miteinbezogen
    startverts = triangle_vectors_of_stl[tri_id_sorted[0]] #Comment_DB: startverts is start vertices
    endverts = triangle_vectors_of_stl[tri_id_sorted[-1]]

    dist_startverts = []
    dist_endverts = []
    for i in range(3):
        dist_startverts.append(distance(startverts[i], center_point_of_cloud_weighted))
        dist_endverts.append(distance(endverts[i], center_point_of_cloud_weighted))
    startvert_3d = startverts[dist_startverts.index(max(dist_startverts))]
    endvert_3d = endverts[dist_endverts.index(max(dist_endverts))]

    # In Ebene projiziert:
    startvert_proj = project_pointtoplane(startvert_3d, trendline_z_axis, center_point_of_cloud_weighted)
    endvert_proj = project_pointtoplane(endvert_3d, trendline_z_axis, center_point_of_cloud_weighted)

    # auf Trendline projiziert
    startvert_trendline = project_pointtoline(startvert_3d, center_point_of_cloud_weighted + trendline_x_axis,
                                              center_point_of_cloud_weighted + 2 * trendline_x_axis)
    endvert_trendline = project_pointtoline(endvert_3d, center_point_of_cloud_weighted + trendline_x_axis,
                                            center_point_of_cloud_weighted + 2 * trendline_x_axis)

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
                                                          center_point_of_cloud_weighted))
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

    for i in range(len(triangle_vectors_of_stl)):
        sorted_centerpoints.append(tri_centerpoints[tri_id_sorted[i]])
    sorted_centerpoints = np.asarray(sorted_centerpoints)

    # Start- und Endkante des Patches finden:
    # -> wenn bei großen Randdreiecken nur der Mittelpunkt benutzt wird, wird das Tape zu kurz. Daher werden
    # bei den Randdreiecken noch die Eckpunkte miteinbezogen
    startverts = triangle_vectors_of_stl[tri_id_sorted[0]]
    endverts = triangle_vectors_of_stl[tri_id_sorted[-1]]
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

    patch_visual = mplot3d.art3d.Poly3DCollection(triangle_vectors_of_stl, linewidths=1, alpha=0.5, edgecolor=[1, 1, 1], label ='Geometriebereich') #Comment_DB: added edgecolor to make the edges visible

    axes.scatter([999999990],[9999999900],[9999999900],linewidths=0.0001, alpha = 0.5, label = "Geometriebereich") #Comment_DB: label in legend

    axes.scatter(datamean[0],datamean[1],datamean[2],c='g')

    pc_axes=np.asarray(trendline)

    # TRENDLINE
    axes.plot(*pc_axes.T,label='Trendlinie', c='red') #Comment_DB: *pc_axes is *args, and .T is np.transpose
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
