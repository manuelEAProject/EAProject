import numpy as np
import math
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from scipy.signal import savgol_filter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from pyquaternion import Quaternion

equidistant_step_size = 3
percentile_pc = 5
max_points_in_pc = 10000

#Berechnet den Abstand von 2 Punkten
def calc_distance_between_two_points(p1, p2):
    distance = np.linalg.norm(p2-p1)
    return distance
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
def calc_start_end_point_3D_from_stl_triangl_vector(center_point_of_cloud_weighted, triangle_ID_sorted_on_trendline, triangle_vectors_of_stl):
    # -> wenn bei großen Randdreiecken nur der Mittelpunkt benutzt wird, wird das Tape zu kurz. Daher werden
    # bei den Randdreiecken noch die Eckpunkte miteinbezogen
    startverts = triangle_vectors_of_stl[triangle_ID_sorted_on_trendline[0]]  # Comment_DB: startverts is start vertices
    endverts = triangle_vectors_of_stl[triangle_ID_sorted_on_trendline[-1]]
    dist_startverts = []
    dist_endverts = []
    for i in range(3):
        dist_startverts.append(calc_distance_between_two_points(startverts[i], center_point_of_cloud_weighted))
        dist_endverts.append(calc_distance_between_two_points(endverts[i], center_point_of_cloud_weighted))
    startvert_3d = startverts[dist_startverts.index(max(dist_startverts))]
    endvert_3d = endverts[dist_endverts.index(max(dist_endverts))]
    return endvert_3d, startvert_3d


####MAIN
def startparam(input_file,poly_order,savgol_window_quotient,max_distance):
    print(
        "Preprocessing läuft... Bei hochaufgelösten stl-Dateien und großen Flächengrößenunterschieden kann es zu "
        "längeren Berechnungszeiten kommen")
    # Comment_DKu_Wenzel: global defined so they dont have to be calculated multiple times. Stay always the same for same stl file

    #Read in the file:
    patch_vectors_of_stl_input = mesh.Mesh.from_file(input_file) #Comment_DB: stl mesh

    stl_normals = (patch_vectors_of_stl_input.normals)

    global triangle_vectors_of_stl
    triangle_vectors_of_stl = patch_vectors_of_stl_input.vectors #Comment_DB: triangle edges (wireframe)
    # Number of triangles, better then using always different len(list)


    global num_of_triangles
    num_of_triangles = len(triangle_vectors_of_stl)

    tri_normals = calc_tri_normals_from_stl(stl_normals,triangle_vectors_of_stl)
    tri_areas = calc_tri_areas(triangle_vectors_of_stl)

    global tri_centerpoints
    tri_centerpoints = calc_tri_centerpoints(triangle_vectors_of_stl)  # Comment_DKu_Wenzel: basically the unweighted point cloud

    global tri_corner_points
    # Die Dreiecksmittelpunkte werden in Ebene der Trendline projiziert
    tri_corner_points = calc_tri_corner_points(triangle_vectors_of_stl)

    global avg_tri_normal_weighted
    avg_tri_normal_weighted = calc_avg_tri_norm_weighted_by_area(tri_areas, tri_normals)

    #Creating pointcloud:
    point_cloud_tri_centerpoints_weighted = calc_patch_pointcloud_weighted_by_area(tri_areas, tri_centerpoints) #!!!!!Comment_DB: The blue points!
    # Calculate the mean of the points, i.e. the 'center' of the cloud
    global center_point_of_cloud_weighted
    center_point_of_cloud_weighted = point_cloud_tri_centerpoints_weighted.mean(axis=0)  # Mean of x,y,z-Values

    global trendline_global_KOS
    trendline_x_axis, trendline_y_axis, trendline_z_axis = calc_trendline_axis_with_svd(point_cloud_tri_centerpoints_weighted, center_point_of_cloud_weighted)



    # Wenn trendline x-Achse in negative x-Richtung definiert wurde von svd
    if trendline_x_axis[0]<0: # Rotation von 180° um y-Achse
        trendline_x_axis= -(trendline_x_axis)
        trendline_z_axis= -(trendline_z_axis)

    trendline_global_KOS = np.vstack((trendline_x_axis, trendline_y_axis, trendline_z_axis))


    # Creating trendline
    global trendline
    trendline = calc_trendline(point_cloud_tri_centerpoints_weighted, center_point_of_cloud_weighted, trendline_x_axis)

    # Creating trendline_projection_points:
    trendline_projection_points_tri_centerpoints = project_tri_centerpoints_to_trendline(tri_centerpoints)
    triangle_ID_sorted_on_trendline = sort_tri_id_by_trendline(trendline_projection_points_tri_centerpoints,trendline)

    #Sorted list of triangle points with trendline_projection_points:
    global sorted_projection_points_tri_centerpoints, sorted_centerpoints
    sorted_projection_points_tri_centerpoints = sort_list_by_tri_id(trendline_projection_points_tri_centerpoints,triangle_ID_sorted_on_trendline)
    sorted_centerpoints = sort_list_by_tri_id(tri_centerpoints,triangle_ID_sorted_on_trendline)

    # Start- und Endkante des Patches finden:
    endvert_3d, startvert_3d = calc_start_end_point_3D_from_stl_triangl_vector(center_point_of_cloud_weighted, triangle_ID_sorted_on_trendline,
                                                                               triangle_vectors_of_stl)


    # In Ebene projiziert:
    global startpoint_project_to_trendline_plan
    startpoint_project_to_trendline_plan = project_pointtoplane(startvert_3d, trendline_z_axis, center_point_of_cloud_weighted)
    endpoint_project_to_trendline_plan = project_pointtoplane(endvert_3d, trendline_z_axis, center_point_of_cloud_weighted)

    # auf Trendline projiziert
    startpoint_project_to_trendline = project_pointtoline(startvert_3d, trendline[0], trendline[1])
    endpoint_project_to_trendline = project_pointtoline(endvert_3d, trendline[0], trendline[1])

    # Für Start und Endpunkte müssen auch die jeweiligen Projektionen auf die Trendline berechnet werden, um den Abstand
    # zu berechnen (COMMENT_DB: insertion of points into sorted_projection_points list)
    sorted_projection_points_tri_centerpoints = np.insert(sorted_projection_points_tri_centerpoints, 0, startpoint_project_to_trendline, axis=0)
    sorted_projection_points_tri_centerpoints = np.concatenate((sorted_projection_points_tri_centerpoints, [endpoint_project_to_trendline]))

    tri_centerpoints_projected_to_trendline_plane = project_tri_centerpoints_to_trendline_plane(endpoint_project_to_trendline_plan,sorted_centerpoints, startpoint_project_to_trendline_plan, trendline_z_axis)


    #tri_corner_and_center_points = np.concatenate((tri_centerpoints, tri_corner_points))
    global tri_corner_points_projected_to_trendline
    tri_corner_points_projected_to_trendline = project_tri_centerpoints_to_trendline(tri_corner_points)

    global tri_corner_points_projected_to_trendline_plane
    tri_corner_points_projected_to_trendline_plane = project_tri_corner_points_to_trendline_plane(tri_corner_points,
                                                                                                  trendline_z_axis)
    global y_list, x_list

    global new_bending_direction_points_global_KOS



    #Interpolation
    grid_lin_interpolation, x_list_interpol_trendline, y_list_interpol_trendline, new_bending_direction_points_global_KOS = interpolate_start_geometrie(max_distance,poly_order, savgol_window_quotient)

    draw_start_point = True
    if draw_start_point:
        startpoint_project_to_trendline_plan = new_bending_direction_points_global_KOS[0][:]
        endpoint_project_to_trendline_plan = new_bending_direction_points_global_KOS[-1][:]


        # Global trendline y-Achse ist lokale z-Achse. -> z_global = - y_lokal für richtige Orientierung des KOS

        trendline_x_axis = trendline_new_direction_global_KOS[0][:]
        trendline_z_axis = -trendline_new_direction_global_KOS[1][:]
        trendline_y_axis = trendline_new_direction_global_KOS[2][:]


    #Start und Endwerte meistens NaN!
    i=0
    while y_list_interpol_trendline[i] != y_list_interpol_trendline[i] : i+=1
    j = -1
    while y_list_interpol_trendline[j] != y_list_interpol_trendline[j]:  j-=1

    x_list = x_list_interpol_trendline[i:j]
    y_list = y_list_interpol_trendline[i:j]

    #x_list = calc_x_values_from_projectetion_on_trendline(tri_corner_points_projected_to_trendline)

    #y_list = calc_y_values_from_projection_points(tri_corner_points_projected_to_trendline, trendline_y_axis, tri_corner_points_projected_to_trendline_plane)

    # Comment_DKu_Wenzel: An dieser Stelle werden die projezierten Punkte vom globalen KOS in ein lokales KOS umgewandelt
    # x-Werte: Abstand zwischen den sorted_projection_points

    #x_list = calc_x_values_from_projectetion_on_trendline(sorted_projection_points_tri_centerpoints)

    # xy-Distanzplot
    # berechnet den Abstand der auf xy-Ebene projizierten Dreiecksmittelpunkte zur Trendline

    #y_list = calc_y_values_from_projection_points(sorted_projection_points_tri_centerpoints, trendline_y_axis, tri_centerpoints_projected_to_trendline_plane)
    # Funktion des xy-Abstands über die Länge der Trendline. Außerdem werden linear Punkte aufgefüllt, um eine
    # nahezu äquidistante Schrittgröße zu erreichen

    global x_y_points_filled_up
    x_y_points_filled_up = calc_x_z_points_filled_up(x_list, y_list)

    # Comment DKu_Wenzel: y_List aus Inpolierter Fläche besser als Anfangswert!! Kein Savgol Filter mehr??!!
    # Geglättete y-Werte mit SavitzkyGolay
    global y_smooth
    y_smooth = smooth_savgol(x_y_points_filled_up, poly_order, savgol_window_quotient)


    # 2D Knicklinie: Start - und Endpunkte;
    global bend_pts_xy
    bend_pts_xy = calc_bend_pts(max_distance, x_y_points_filled_up, y_list)

    ###### Startparameter extrahieren #####
    ###Start_r_atstart in 2D & 3D### COMMENT_DB: NEW DEFINITION
    Start_r_2d_atstart = bend_pts_xy[1] - bend_pts_xy[0]

    # todo: Start Lösungsrichtung muss hier angepasst werden: trendline x,y,z_axis
    Start_r_3d_atstart = Start_r_2d_atstart[0] * trendline_x_axis + \
                         Start_r_2d_atstart[1] * trendline_y_axis
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

    # L_aim of savgol curve: #Comment_DKu_Wenzel: Kein Savgol mehr! Nimmt kurve von Interpolierter Fläche!!
    L_aim = 0
    for i in range(1, len(x_y_points_filled_up)):
        p0=np.array([x_y_points_filled_up[i - 1][0], x_y_points_filled_up[i - 1][1]])
        p1=np.array([x_y_points_filled_up[i][0], x_y_points_filled_up[i][1]])
        L_aim =L_aim + np.linalg.norm(p1-p0)



    start_parameter = [l_list, L_aim, beta_list, startpoint_project_to_trendline_plan,
                       endpoint_project_to_trendline_plan, Start_r_3d_atstart, Start_n_3d_atstart]
    return start_parameter

def calc_tri_normals_from_stl(stl_normals,triangle_vectors_of_stl):
    normals=[]
    #We generate our own normals with the id_list. Notice that the triangle order of stl_normals and the vertices
    #(triangles) is not the same

    #Comment_DKu_Wenzel: We are just taking the length of tri_id_sorted. So we iterate normal and not the different id
    for i in range(num_of_triangles):

        v1= triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][1]
        v2= triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][2]
        n=np.cross(v1,v2)
        n=(1 / np.linalg.norm(n)) * n
        normals.append(n)

    normals=np.asarray(normals)


    # Normalenvektoren werden immer in positive z-Richtung ausgerichtet

    # the following average stl_normal always point at the outside of the object:
    avg_stl_normal = sum(stl_normals) / num_of_triangles
    # average of the created normals:
    avg_sorted_normal = sum(normals) / num_of_triangles


    true_when_stl_and_tri_normal_not_same_direction = avg_sorted_normal[0] * avg_stl_normal[0] < 0
    true_when_z_from_tri_normal_neg = avg_sorted_normal[2] < 0

    # Comment_DKu_Wenzel: Die if Bedingung ist seltsam
    if true_when_stl_and_tri_normal_not_same_direction or true_when_z_from_tri_normal_neg:
        normals=np.negative(normals)

    return normals
def calc_tri_centerpoints(triangle_vectors_of_stl):
    tri_centerpoints=[]
    for i in range(num_of_triangles):
        center=np.array([(triangle_vectors_of_stl[i][0][0] + triangle_vectors_of_stl[i][1][0] + triangle_vectors_of_stl[i][2][0]) / 3, (triangle_vectors_of_stl[i][0][1] + triangle_vectors_of_stl[i][1][1] + triangle_vectors_of_stl[i][2][1]) / 3, (triangle_vectors_of_stl[i][0][2] + triangle_vectors_of_stl[i][1][2] + triangle_vectors_of_stl[i][2][2]) / 3])
        tri_centerpoints.append(center)

    tri_centerpoints=np.asarray(tri_centerpoints)
    return tri_centerpoints
def calc_tri_areas(triangle_vectors_of_stl):
    #Berechnung der Dreiecksflächen und Speichern in einer Liste
    tri_surface_area = []
    for i in range(num_of_triangles):
        tri_surface_area.append(0.5 * (
            np.linalg.norm((triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][1]) - (triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][2]))))
    tri_surface_area = np.asarray(tri_surface_area)
    return tri_surface_area
def calc_avg_tri_norm_weighted_by_area(tri_areas,triangle_normals):
    weighted_norms = []
    for i in range(num_of_triangles):
        weighted_norms.append((tri_areas[i] / sum(tri_areas)) * triangle_normals[i])
    avg_tri_normal_weighted = sum(weighted_norms)

    return avg_tri_normal_weighted
def calc_patch_pointcloud_weighted_by_area(tri_areas, tri_centerpoints):
    #
    centerpoints_weights_area_tri = calc_weights_for_center_points_by_percentil_area(tri_areas)

    # für jedes Dreieck werden die Mittelpunkte so oft gewertet (in die Punktwolke getan) wie es in centerpoints_weights_area_tri steht
    pointcloud_weighted=[]

    for i in range(num_of_triangles):
        for j in range(centerpoints_weights_area_tri[i]):
            pointcloud_weighted.append(tri_centerpoints[i])

    pointcloud_weighted=np.asarray(pointcloud_weighted)

    return pointcloud_weighted
def calc_weights_for_center_points_by_percentil_area(tri_areas):
    # In dieser Funktion wird eine Punktwolke aus den Dreiecksmittelpunkten der stl Datei generiert. Da mit dieser
    # Punktwolke später die Hauptwertzerlegung erfolgt, müssen große Dreiecke gegenüber kleinen Dreiecken stärker
    # gewichtet werden. Dies geschieht indem das Verhältnis der Dreiecksflächen zueinander berechnet wird und die
    # Mittelpunkte großer Dreiecke öfter zur Punktwolke hinzugefügt werden als die Mittelpunkte kleiner Dreiecke.

    ###Gewichtung großer Dreiecke
    # 1) Berechnung der Dreiecksflächen und Speichern in einer Liste
    area_whole_patch = sum(tri_areas)
    # 2) Die Dreiecksflächen werden miteinander verglichen, um später große Dreiecke gegenüber kleinen stärker zu
    # gewichten. Als "kleinste" Referenzfläche wird das 10-er Quantil der Flächen genommen (90% aller anderen Dreiecke
    # sind größer). In centerpoints_weights_area_tri wird für jedes Dreieck berechnet, um welchen Faktor es größer als das
    # Referenzdreicek ist (der Faktor wird aufgerundet). Zur Pointcloud wird von jedem Dreieck der Mittelpunkt so oft
    # hinzugefügt wie hoch der Flächenfaktor aus der centerpoints_weights_area_tri ist (mindestens einmal).
    # Um zu große Rechenzeiten zu vermeiden wird im Vorhinein abgeschätzt wie viele Punkte die Punktwolke
    # haben wird und bei Überschreiten eines Grenzwerts( max_points_in_pc) wird das Programm abgebrochen.
    lower_percentil_area = np.percentile(tri_areas, percentile_pc)
    estimated_number_points_in_pc = math.ceil(area_whole_patch / lower_percentil_area)

    # Abbruchbedingung: in der Pointcloud sollen maximal "max_points_in_pc" berechnet werden.
    if max_points_in_pc < estimated_number_points_in_pc:
        print("ERROR: Please use a .stl-object with reduced resolution ")
        print("Number of triangles: ", num_of_triangles)
        print("Estimated number of points in pointcloud:", estimated_number_points_in_pc)
        print("Allowed number of points in pointcloud:", max_points_in_pc)
        exit(1)

    # Im Folgenden wird jedes Dreieck mit dem kleinsten Dreieck verglichen und in centerpoints_weights_area_tri festgehalten, wie
    # oft das kleinste Dreieck in das jeweilige Dreieck hineinpasst.
    centerpoints_weights_area_tri = []
    for i in range(num_of_triangles):
        centerpoints_weights_area_tri.append(math.ceil(tri_areas[i] / lower_percentil_area))

    return centerpoints_weights_area_tri
def calc_trendline(patch_pc_weighted, center_point_of_cloud_weighted, trendline_x_axis):
    # Gerade durch Punktwolke mit kleinstem Abstand zu allen Punkten (COMMENT_DB: pc --> point cloud)
    # ->Hauptkomponentenanalyse: https://de.wikipedia.org/wiki/Hauptkomponentenanalyse/Principal Component Analysis(PCA)
    # aka Singulärwertzerlegung / SVD
    # Generate some data that lies along a line
    #SOURCE: https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    #Explanation: https://www.tutorialspoint.com/scipy/scipy_linalg.htm
    maxvals = np.amax(patch_pc_weighted, axis=0)

    # It's a straight line, so we only need 2 points.
    linepts = trendline_x_axis* np.mgrid[-2*max(maxvals):2*max(maxvals):2j][:, np.newaxis]

    # shift by the mean to get the line in the right place
    linepts += center_point_of_cloud_weighted #COMMENT_DB: linepts = linepts + center_point_of_cloud_weighted (shifting it in the same direction as the vv[0] direction vector)

    return linepts
def project_tri_centerpoints_to_trendline(tri_centerpoints):
    #SOURCE: https://gamedev.stackexchange.com/questions/72528/how-can-i-project-a-3d-point-onto-a-3d-line
    trendline_projection = []
    for i in range(len(tri_centerpoints)):
        trendline_projection.append(project_pointtoline(tri_centerpoints[i], trendline[0], trendline[1]))
    trendline_projection=np.asarray(trendline_projection)

    return trendline_projection
def calc_trendline_axis_with_svd(patch_pc_weighted, center_point_of_cloud_weighted):
    # Do Principal Component Analysis(PCA) on the mean-centered data. AKA SVD
    # The first principal component contains [uu, dd, vv] , where vv[0] is the direction
    first_principal_components_pc_weighted = np.linalg.svd(patch_pc_weighted - center_point_of_cloud_weighted)
    # Definition der Hauptachsen
    trendline_x_axis = first_principal_components_pc_weighted[2][0] # first_principal_components_pc_weighted[2][0]: is direction of trendline
    trendline_x_axis = (1 / np.linalg.norm(trendline_x_axis)) * trendline_x_axis
    # avg_tri_norm ist nicht senkrecht zur x-Achse
    # von pcc + avg_tri_norm und zurück auf x-Achse projizieren
    trendline_avg_norm_point = center_point_of_cloud_weighted + np.dot(avg_tri_normal_weighted,
                                                                       trendline_x_axis) / np.dot(trendline_x_axis,
                                                                                                  trendline_x_axis) * trendline_x_axis
    # y-Achse ist verbindung von pcc+avg_tri_norm mit dem projizierten Punkt
    trendline_y_axis = (center_point_of_cloud_weighted + avg_tri_normal_weighted) - trendline_avg_norm_point
    trendline_y_axis = (1 / np.linalg.norm(trendline_y_axis)) * trendline_y_axis
    trendline_z_axis = np.cross(trendline_x_axis, trendline_y_axis)
    return trendline_x_axis, trendline_y_axis, trendline_z_axis
def sort_tri_id_by_trendline(trendline_projection_points_tri_centerpoints,trendline):
    # Sortiert die projizierten Punkte der Trendlinie entlang und gibt so eine Reihenfolge der
    # Dreiecke entlang der Trendlinienrichtung an. Es wird eine geordnete Liste der
    # Dreiecks-IDs zurückgegeben
    distance_startpoint_trendline_to_projectionpoint = []
    for i in range(num_of_triangles):
        distance = np.array([np.linalg.norm((trendline_projection_points_tri_centerpoints[i] - trendline[0])), i])
        distance_startpoint_trendline_to_projectionpoint.append(distance)

    # Put distances in a array for sorting triangles by size and keeping the global id of triangles
    distance_startpoint_trendline_to_projectionpoint = np.asarray(distance_startpoint_trendline_to_projectionpoint)

    sorted_distance_by_size = distance_startpoint_trendline_to_projectionpoint[distance_startpoint_trendline_to_projectionpoint[:, 0].argsort()]
    # Get id and put them in a list (from array
    id_list = list(map(int, sorted_distance_by_size[:,1]))

    return id_list
def sort_list_by_tri_id(list_to_sort,triangle_ID_sorted_on_trendline):
    # Liste mit den sortierten auf die Trendline projizierten Punkte
    #Die auf die Trendline projizierten Mittelpunkte der Dreiecke werden entlang der Trendline sortiert
    sorted_trendline = []
    for i in range(num_of_triangles):
        sorted_trendline.append(list_to_sort[triangle_ID_sorted_on_trendline[i]])
    sorted_trendline=np.asarray(sorted_trendline)
    return sorted_trendline
def smooth_savgol(equidistant_data_set,polynom_order,savgol_window_quotient):
    savgol_window = int(len(equidistant_data_set) / savgol_window_quotient)
    #polynom_order = 3 # Comment_DKu_Wenzel: Warum?
    # window must be impair:
    if savgol_window % 2 == 0:
        savgol_window += 1
    return savgol_filter(equidistant_data_set[:,1], savgol_window, polynom_order) # (data,window size, polynomial order)
def find_nearest(array, value):
# find the index of the closest value in an array to a given value
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
def calc_bend_pts(max_distance, x_y_points_filled_up, y_smooth):

    bend_pts_xy = []
    bend_pts_xy.append([x_y_points_filled_up[0][0], y_smooth[0]])  # Comment_DB: start point 2D (x coord, y coord)
    bend_pts_xy.append([x_y_points_filled_up[-1][0], y_smooth[-1]])  # Comment_DB: end point 2D (x coord, y coord)
    bend_pts_xy = np.asarray(bend_pts_xy)

    # Einfügen von weiteren Knickpunkten durch Finden von großen Abweichungen zur Kurve:
    insert_pts = True
    while insert_pts:
        points_on_line_between_bends_filled_up = []
        points_on_line_between_bends_filled_up.append([bend_pts_xy[0][0], bend_pts_xy[0][1]])  # Comment_DB: only the first bend point (starting point at edge) appended to bend points curve list

        points_on_line_between_bends_filled_up = calc_points_on_line_between_bends_filled_up(bend_pts_xy,
                                                                                              points_on_line_between_bends_filled_up, x_y_points_filled_up)

        max_divergence = calc_point_of_max_divergence_between_smooth_and_lin_curve(
            points_on_line_between_bends_filled_up, x_y_points_filled_up, y_smooth)

        # Comment_DB: We know at which x-coord of points_on_line_between_bends_filled_up the max_divergence happens --> counter i

        bend_pts_xy = np.insert(bend_pts_xy, -1,
                                np.array([points_on_line_between_bends_filled_up[max_divergence[1]][0], y_smooth[max_divergence[1]]]),
                                axis=0)  # Comment_DB: insert a corner at x coord (counter i) and y coord (counter i) of max divergence
        bend_pts_xy = bend_pts_xy[bend_pts_xy[:, 0].argsort()]  # Comment_DB: Bend points sorted in an array

        # no further points, if the chosen maximum distance is not surpassed
        if max_divergence[0] < max_distance:  # Comment_DB: This implies that there will be one extra bend, as the above code will have executed already, max_distance: User Input
            insert_pts = False

    return bend_pts_xy
def calc_point_of_max_divergence_between_smooth_and_lin_curve(points_on_line_between_bends_filled_up, x_y_points_filled_up, y_smooth):
    # Comment_DB: curve_divergence in terms of y-distance # Größte Abweichung von geglätteter Kurve: (COMMENT_DB: By now all the points in the above (linear) line have been appended)
    curve_divergence_y = []
    for i in range(len(points_on_line_between_bends_filled_up)):
        curve_divergence_y.append([points_on_line_between_bends_filled_up[i][0], (
                (points_on_line_between_bends_filled_up[i][0] - x_y_points_filled_up[i][0]) ** 2 + (
                points_on_line_between_bends_filled_up[i][1] - y_smooth[i]) ** 2) ** 0.5])  # Comment_DB: (x-coord vs. change in y-coord) take the x coord and y-distance between linear curve and sav-gol curve and append
    curve_divergence_y = np.asarray(curve_divergence_y)
    max_divergence = max([(v, i) for i, v in enumerate(curve_divergence_y[:, 1])])  # Comment_DB: returns distance, counter (Uses new curve_divergence)
    return max_divergence
def calc_points_on_line_between_bends_filled_up(bend_pts_xy, bend_pts_xy_curve, x_y_points_filled_up):
    j = 1  # Comment_DB: at this point, bend_pts_xy curve only has the starting point in it, thus j = 1 is the number of points in the list. j = 1 is also the index of the NEXT point!
    for i in range(1, len(bend_pts_xy)):  # Comment_DB: len(bend_pts_xy) is 2 for first iteration

        slope_between_bends = (bend_pts_xy[i - 1][1] - bend_pts_xy[i][1]) / (bend_pts_xy[i - 1][0] - bend_pts_xy[i][0])
        while bend_pts_xy_curve[-1][0] < bend_pts_xy[i][0]:  # Comment_DB: while last x coord VALUE less than ith x coord VALUE in bend_pts_xy (If greater, then that means last point is reached)
            y_add = bend_pts_xy_curve[-1][1] + slope_between_bends * (
                        x_y_points_filled_up[j][0] - x_y_points_filled_up[j - 1][0])  # Comment_DB: y = b + mx (finds next change in y linearly --> Produces a linear plot until end point at edge!!)
            bend_pts_xy_curve.append([x_y_points_filled_up[j][0], y_add])  # Comment_DB: append the NEXT point into the list
            j = j + 1  # Comment_DB: NEXT POINT
    bend_pts_xy_curve = np.asarray(bend_pts_xy_curve)  # Comment_DB: This is now one linear curve from start to end point. Everything here is dependent on xy_patch_curve. Below will take divergence into consideration
    return bend_pts_xy_curve
def calc_x_z_points_filled_up(x_list, y_list):
    xy_patch_curve = []  # Comment_DB: Smoothed blue curve points in array!
    for i in range(len(y_list)):
        xy_patch_curve.append([x_list[i], y_list[i]])

    for i in range(1, len(y_list)):

        x_dist = x_list[i] - x_list[i - 1]
        if x_dist > equidistant_step_size:
            additional_steps = math.floor(x_dist / equidistant_step_size)

            for j in range(1, additional_steps + 1):
                xy_patch_curve.append([x_list[i - 1] + j * equidistant_step_size,
                                       (y_list[i - 1] + (y_list[i] - y_list[i - 1]) / (
                                           x_dist) * j * equidistant_step_size)])

    xy_patch_curve = np.asarray(xy_patch_curve)
    # Entlang der x-Werte sortieren
    xy_patch_curve = xy_patch_curve[xy_patch_curve[:, 0].argsort()]

    return xy_patch_curve
def calc_y_values_from_projection_points(sorted_projection_points_tri_centerpoints, trendline_y_axis,tri_centerpoints_projected_to_trendline_plane):
    y_list = []
    for i in range(len(tri_centerpoints_projected_to_trendline_plane)):
        dist = np.linalg.norm(
            sorted_projection_points_tri_centerpoints[i] - tri_centerpoints_projected_to_trendline_plane[i])

        # Vorzeichen ermitteln:
        if ((sorted_projection_points_tri_centerpoints[i] + dist * trendline_y_axis)[0] -
            tri_centerpoints_projected_to_trendline_plane[i][0]) ** 2 < \
                ((sorted_projection_points_tri_centerpoints[i] - dist * trendline_y_axis)[0] -
                 tri_centerpoints_projected_to_trendline_plane[i][0]) ** 2:
            y_list.append(dist)
        else:
            y_list.append(-dist)
    y_list = np.asarray(y_list)

    return y_list
def project_tri_centerpoints_to_trendline_plane(endpoint_project_to_trendline_plan,sorted_centerpoints, startpoint_project_to_trendline_plan, trendline_z_axis):
    tri_centerpoints_projected_to_trendline_plane = [startpoint_project_to_trendline_plan]
    for i in range(num_of_triangles):
        tri_centerpoints_projected_to_trendline_plane.append(
            project_pointtoplane(sorted_centerpoints[i], trendline_z_axis,
                                 center_point_of_cloud_weighted))
    tri_centerpoints_projected_to_trendline_plane.append(endpoint_project_to_trendline_plan)
    tri_centerpoints_projected_to_trendline_plane = np.asarray(tri_centerpoints_projected_to_trendline_plane)
    return tri_centerpoints_projected_to_trendline_plane
def calc_x_values_from_projectetion_on_trendline(sorted_projection_points_tri_centerpoints):
    x_list = []
    for i in range(len(sorted_projection_points_tri_centerpoints)):
        x = np.linalg.norm(
            (sorted_projection_points_tri_centerpoints[0] - sorted_projection_points_tri_centerpoints[i]))
        x_list.append(x)
    x_list = np.asarray(x_list)
    return x_list

def calc_tri_corner_points(triangle_vectors_of_stl):
    tri_corner_points=[]
    for i in range(num_of_triangles):
        triangel = triangle_vectors_of_stl[i]
        for j in range(3):
            tri_corner_points.append(triangel[j])

    # Doppelte Ecken gelöscht
    tri_corner_points = np.unique(tri_corner_points, axis=0)

    return tri_corner_points
def project_tri_corner_points_to_trendline_plane(tri_corner_points, trendline_z_axis):
    tri_cornerpoints_projected_to_trendline_plane = []
    for i in range(len(tri_corner_points)):
        tri_cornerpoints_projected_to_trendline_plane.append(
            project_pointtoplane(tri_corner_points[i], trendline_z_axis,
                                 center_point_of_cloud_weighted))

    tri_cornerpoints_projected_to_trendline_plane = np.asarray(tri_cornerpoints_projected_to_trendline_plane)
    return tri_cornerpoints_projected_to_trendline_plane


def interpolate_start_geometrie(max_distance,poly_order, savgol_window_quotient, grid_resolution=2500j):
    grid_ressolution_int = int(abs(grid_resolution))

    grid_x, grid_y, max_x, max_y, min_x, min_y, z_grid_values_linear = interpolate_geometrie(grid_resolution)

    show_interpolation_and_draw_start_end_points(max_x, max_y, min_x, min_y, z_grid_values_linear)

    # todo: Hier wird bisher von der trendlinie die Biegelienie entnommen. Können wir eine "Schräge_Ebene auswählen?"
    ### Comment_DKu_Wenzel: Vorsicht!! Hier gibt die z-Achse das Höhenprofil an. In start_parameter ist es y.

    # Straight line through y=0
    y_0_grid_point_index = np.asarray(np.round(max_y / (max_y - min_y) * grid_ressolution_int), dtype=np.int32)
    x_0_grid_point_index = np.asarray(np.round(max_x / (max_x - min_x) * grid_ressolution_int), dtype=np.int32)
    y_value_at_gridpoint = grid_y[0,-y_0_grid_point_index] # Should be around 0. Most of the time not exact 0
    z_values_at_y_0 = z_grid_values_linear[:, -y_0_grid_point_index]

    # Tilted line, line given by two points xdata and ydata
    new_bending_direction_points_tilted_KOS, new_bending_direction_points_tilted_KOS_left, new_bending_direction_points_tilted_KOS_right, new_bending_direction_points_global_KOS, x_values, x_values_trim, y_values_trim, bend_pts_xz_global = calc_tilted_bending_points(
        grid_ressolution_int, grid_x, max_y, min_y ,max_x, min_x, y_0_grid_point_index, x_0_grid_point_index, z_grid_values_linear,xdata,ydata,max_distance,poly_order, savgol_window_quotient)




    ###2D-xy-PLOT
    pyplot.figure()

    plt.subplot(221)
    plt.imshow(z_grid_values_linear.T, extent=(min_x, max_x, min_y, max_y), origin='lower')
    plt.plot(x_values_trim, y_values_trim, 'bo', linewidth=2.0, label='Schnitt')
    plt.title('Linear')

    plt.subplot(222)
    plt.plot(new_bending_direction_points_tilted_KOS[:, 0], new_bending_direction_points_tilted_KOS[:, 2], 'bo', linewidth=2.0,
             label='z_values_new_bending_direction')
    plt.plot(bend_pts_xz_local[:, 0], bend_pts_xz_local[:, 1], color='green', linewidth=3.0,
             label='lineare Angleichung')  # Streckenweise linear (nur Eckpunkte)

    plt.subplot(223)
    plt.title('left')
    plt.plot(new_bending_direction_points_tilted_KOS_left[:, 0], new_bending_direction_points_tilted_KOS_left[:, 2], 'bo', linewidth=2.0, label='z_values_at_y_0')
    plt.plot(bend_pts_xz_local_left[:, 0], bend_pts_xz_local_left[:, 1], color='green', linewidth=3.0,
             label='lineare Angleichung')  # Streckenweise linear (nur Eckpunkte)

    plt.subplot(224)
    plt.title('right')
    plt.plot(new_bending_direction_points_tilted_KOS_right[:, 0], new_bending_direction_points_tilted_KOS_right[:, 2], 'bo', linewidth=2.0, label='z_values_at_y_0')
    plt.plot(bend_pts_xy_local_right[:, 0], bend_pts_xy_local_right[:, 1], color='green', linewidth=3.0,
             label='lineare Angleichung')  # Streckenweise linear (nur Eckpunkte)
    plt.show()



    return z_grid_values_linear, new_bending_direction_points_tilted_KOS[:, 0], new_bending_direction_points_tilted_KOS[:, 2], new_bending_direction_points_global_KOS

def calc_tilted_bending_points(grid_ressolution_int, grid_x, max_y, min_y, max_x, min_x, y_0_grid_point_index, x_0_grid_point_index, z_grid_values_linear, xdata, ydata,max_distance,poly_order, savgol_window_quotient):
    ###### Schräge Gerade mit y(x)
    # Wir brauchen 2 Geraden:   • Einmal die mathematische beschreibung mit Coordinaten für Plot.
    #                           • Einmal die mit Indizies für die Extraction aus dem Grid

    # Schrittweite
    dy = (max_y - min_y) / grid_ressolution_int
    dx = (max_x - min_x) / grid_ressolution_int



    # Steigung(x_slope) und y-Achsenabschnitt(y_intercept) mit Leastsquare, y = x_slope*x + y_intercept
    A = np.vstack([xdata, np.ones(len(xdata))]).T
    x_slope, y_intercept = np.linalg.lstsq(A, ydata,rcond=None)[0]

    x_values = grid_x[:, 0]


    # x_y_z_Axis of new KOS
    trendline_new_direction_current_KOS, x_trendline_new_direction, z_trendline_new_direction = calc_local_trendline_KOS(x_slope)

    end_point_drawn, start_point_drawn, x_start_index, x_end_index= calculate_3D_Start_End_from_xdata_ydata(dx, dy, grid_ressolution_int,
                                                                                 x_0_grid_point_index, xdata,
                                                                                 y_0_grid_point_index, ydata,
                                                                                 z_grid_values_linear)


    # Eckpunkte
    width =5
    delta_length_start_bend = 0
    delta_length_end_bend = 0

    end_point_drawn_left, start_point_drawn_left, x_end_index_left, x_start_index_left = calc_start_end_point_side(True,delta_length_end_bend,delta_length_start_bend,
                                                                                                                   dx, dy,end_point_drawn,grid_ressolution_int,
                                                                                                                   start_point_drawn, width, x_0_grid_point_index, x_trendline_new_direction,
                                                                                                                   y_0_grid_point_index, z_grid_values_linear,z_trendline_new_direction)

    end_point_drawn_right, start_point_drawn_right, x_end_index_right, x_start_index_right = calc_start_end_point_side(False,delta_length_end_bend,delta_length_start_bend,
                                                                                                                   dx,dy,end_point_drawn,grid_ressolution_int, start_point_drawn,
                                                                                                                   width, x_0_grid_point_index,x_trendline_new_direction, y_0_grid_point_index,
                                                                                                                   z_grid_values_linear, z_trendline_new_direction)



    # Calculate points on line  # Comment_DKu_Wenzel: Rename this function
    new_bending_direction_points_left_current_KOS, y_intercept_left, x_values_indizes_trim_left, x_values_trim_left, y_values_indizes_trim_left, y_values_trim_left = calc_new_direction_bending_points(
        start_point_drawn_left, end_point_drawn_left,dy, grid_ressolution_int, x_slope, x_values, y_0_grid_point_index, z_grid_values_linear, x_start_index_left, x_end_index_left)
    ### Rotating 3D-Points to global and tilted KOS
    global new_bending_direction_points_global_KOS_left, start_end_point_drawn_global_KOS_left
    new_bending_direction_points_tilted_KOS_left, new_bending_direction_points_global_KOS_left, start_end_point_drawn_global_KOS_left, trendline_new_direction_global_KOS_left = new_bending_points_in_global_and_tilted_KOS(
        start_point_drawn_left, end_point_drawn_left, y_intercept_left, new_bending_direction_points_left_current_KOS,
        trendline_new_direction_current_KOS)
    global bend_pts_xz_local_left, bend_pts_xz_global_left
    bend_pts_xz_local_left, bend_pts_xz_global_left = calc_local_and_global_bendpoints(max_distance,
                                                                                       new_bending_direction_points_tilted_KOS_left,
                                                                                       poly_order,
                                                                                       savgol_window_quotient,
                                                                                       trendline_new_direction_current_KOS,
                                                                                       y_intercept_left)



    new_bending_direction_points_right_current_KOS, y_intercept_right, x_values_indizes_trim_right, x_values_trim_right, y_values_indizes_trim_right, y_values_trim_right = calc_new_direction_bending_points(
        end_point_drawn_right, start_point_drawn_right, dy, grid_ressolution_int, x_slope, x_values,
        y_0_grid_point_index, z_grid_values_linear, x_start_index_right, x_end_index_right)
    global new_bending_direction_points_global_KOS_right, start_end_point_drawn_global_KOS_right
    new_bending_direction_points_tilted_KOS_right, new_bending_direction_points_global_KOS_right, start_end_point_drawn_global_KOS_right, trendline_new_direction_global_KOS_right = new_bending_points_in_global_and_tilted_KOS(
        start_point_drawn_right, end_point_drawn_right, y_intercept_right, new_bending_direction_points_right_current_KOS,
        trendline_new_direction_current_KOS)
    global bend_pts_xy_local_right, bend_pts_xy_global_right
    bend_pts_xy_local_right, bend_pts_xy_global_right = calc_local_and_global_bendpoints(max_distance,
                                                                                       new_bending_direction_points_tilted_KOS_right,
                                                                                       poly_order,
                                                                                       savgol_window_quotient,
                                                                                       trendline_new_direction_current_KOS,
                                                                                       y_intercept_right)



    new_bending_direction_points_current_KOS, y_intercept, x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim = calc_new_direction_bending_points(
        start_point_drawn, end_point_drawn, dy, grid_ressolution_int, x_slope, x_values, y_0_grid_point_index,
        z_grid_values_linear, x_start_index, x_end_index)

    global start_end_point_drawn_global_KOS, trendline_new_direction_global_KOS
    new_bending_direction_points_tilted_KOS, new_bending_direction_points_global_KOS, start_end_point_drawn_global_KOS, trendline_new_direction_global_KOS = new_bending_points_in_global_and_tilted_KOS(
        start_point_drawn, end_point_drawn, y_intercept, new_bending_direction_points_current_KOS,
        trendline_new_direction_current_KOS)
    global bend_pts_xz_local, bend_pts_xy_global
    bend_pts_xz_local, bend_pts_xy_global = calc_local_and_global_bendpoints(max_distance,
                                                                             new_bending_direction_points_tilted_KOS,
                                                                             poly_order,
                                                                             savgol_window_quotient,
                                                                             trendline_new_direction_current_KOS,
                                                                             y_intercept)





    return new_bending_direction_points_tilted_KOS, new_bending_direction_points_tilted_KOS_left,new_bending_direction_points_tilted_KOS_right, new_bending_direction_points_global_KOS, x_values, x_values_trim, y_values_trim, bend_pts_xy_global


def calc_local_trendline_KOS(x_slope):
    x_trendline_new_direction = np.asarray((1, x_slope, 0), dtype=np.float32)
    x_trendline_new_direction = 1 / np.linalg.norm(x_trendline_new_direction) * x_trendline_new_direction
    y_trendline_new_direction = np.asarray((-x_slope, 1, 0), dtype=np.float32)
    y_trendline_new_direction = 1 / np.linalg.norm(y_trendline_new_direction) * y_trendline_new_direction
    z_trendline_new_direction = np.asarray((0, 0, 1), dtype=np.float32)
    trendline_new_direction_current_KOS = np.vstack(
        (x_trendline_new_direction, y_trendline_new_direction, z_trendline_new_direction))
    return trendline_new_direction_current_KOS, x_trendline_new_direction, z_trendline_new_direction


def calc_start_end_point_side(calc_left_side, delta_length_end_bend, delta_length_start_bend, dx, dy, end_point_drawn,
                              grid_ressolution_int, start_point_drawn, width, x_0_grid_point_index,
                              x_trendline_new_direction, y_0_grid_point_index, z_grid_values_linear,
                              z_trendline_new_direction):

    if calc_left_side:
        Start_point_left = start_point_drawn - np.cross(x_trendline_new_direction,
                                                    z_trendline_new_direction) * width / 2 + delta_length_start_bend * x_trendline_new_direction
        End_point_left = end_point_drawn - np.cross(x_trendline_new_direction,
                                                z_trendline_new_direction) * width / 2 + delta_length_end_bend * x_trendline_new_direction

    else:
        Start_point_left = start_point_drawn + np.cross(x_trendline_new_direction,
                                                        z_trendline_new_direction) * width / 2 + delta_length_start_bend * x_trendline_new_direction
        End_point_left = end_point_drawn + np.cross(x_trendline_new_direction,
                                                    z_trendline_new_direction) * width / 2 + delta_length_end_bend * x_trendline_new_direction

    x_data_left_start_end = (Start_point_left[0], End_point_left[0])
    y_data_left_start_end = (Start_point_left[1], End_point_left[1])
    end_point_drawn_left, start_point_drawn_left, x_start_index_left, x_end_index_left = calculate_3D_Start_End_from_xdata_ydata(
        dx, dy, grid_ressolution_int,
        x_0_grid_point_index, x_data_left_start_end,
        y_0_grid_point_index, y_data_left_start_end,
        z_grid_values_linear)
    return end_point_drawn_left, start_point_drawn_left, x_end_index_left, x_start_index_left

    # todo: left remove
def calc_local_and_global_bendpoints(max_distance, new_bending_direction_points_tilted_KOS, poly_order,
                                     savgol_window_quotient, trendline_new_direction_current_KOS, y_intercept):

    x_z_points_filled_up = calc_x_z_points_filled_up(new_bending_direction_points_tilted_KOS[:, 0],
                                                     new_bending_direction_points_tilted_KOS[:, 2])
    #y_smooth_left = smooth_savgol(x_z_points_filled_up, poly_order, savgol_window_quotient)
    bend_pts_xz_local = calc_bend_pts(max_distance, x_z_points_filled_up,
                                           x_z_points_filled_up[:,1])
    # in tilted KOS y=0. Insert this to bendpoints to get (x,y,z)
    bend_pts_xyz_left = np.insert(bend_pts_xz_local, 1, np.zeros(len(bend_pts_xz_local)), axis=1)
    bend_pts_xyz_global = new_bending_points_tilted_to_global(y_intercept, bend_pts_xyz_left,
                                                                  trendline_new_direction_current_KOS)
    return bend_pts_xz_local, bend_pts_xyz_global

def calculate_3D_Start_End_from_xdata_ydata(dx, dy, grid_ressolution_int, x_0_grid_point_index, xdata,
                                            y_0_grid_point_index, ydata, z_grid_values_linear):
    # Start und Endpunkt der Eingezeichnet wurde

    y_end_index = np.asarray(np.round(np.add(np.divide(ydata[1], dy), (grid_ressolution_int - y_0_grid_point_index))),
                             dtype=np.int32)
    x_end_index = np.asarray(np.round(np.add(np.divide(xdata[1], dx), (grid_ressolution_int - x_0_grid_point_index))),
                             dtype=np.int32)
    y_start_index = np.asarray(np.round(np.add(np.divide(ydata[0], dy), (grid_ressolution_int - y_0_grid_point_index))),
                               dtype=np.int32)
    x_start_index = np.asarray(np.round(np.add(np.divide(xdata[0], dx), (grid_ressolution_int - x_0_grid_point_index))),
                               dtype=np.int32)

    # Left and right starting point can be outside the grit, when they are outside, they get a default value. The z-Data would be needed for the plot.
    if x_start_index < 0 or x_start_index >= grid_ressolution_int or \
            y_start_index < 0 or y_start_index >= grid_ressolution_int:
        z_start_data = z_grid_values_linear[0, 0]
    else: z_start_data = z_grid_values_linear[x_start_index, y_start_index]

    if x_end_index < 0 or x_end_index >= grid_ressolution_int or \
            y_end_index < 0 or y_end_index >= grid_ressolution_int:
        z_end_data = z_grid_values_linear[grid_ressolution_int - 1, grid_ressolution_int - 1]
    else: z_end_data = z_grid_values_linear[x_end_index, y_end_index]


    start_point_drawn = (np.vstack((xdata[0], ydata[0], z_start_data)).T)[0][:]
    end_point_drawn = (np.vstack((xdata[1], ydata[1], z_end_data)).T)[0][:]
    return end_point_drawn, start_point_drawn, x_start_index, x_end_index

def calc_new_direction_bending_points(Start_point, End_point, dy, grid_ressolution_int, x_slope, x_values, y_0_grid_point_index,
                                      z_grid_values_linear, x_start_index, x_end_index):


    # Neuer y-Achsen Abschnitt für Eckpunkte: y = x_slope*x + y_intercept
    y_intercept = Start_point[1] - Start_point[0] * x_slope
    # x-y-Achsenabschnitt in Coordinaten
    # y = x_slope*x + y_intercept

    y_values = np.add(np.multiply(x_values, x_slope), y_intercept)
    # x-y-Achsenabschnitt in Indizies
    x_values_indizes = np.asarray(list(range(grid_ressolution_int)), dtype=np.int32)
    y_values_indizes = np.add(np.divide(y_values, dy), (grid_ressolution_int - y_0_grid_point_index))
    y_values_indizes = np.asarray(np.round(y_values_indizes), dtype=np.int32)


    # Trim x and y valules to start/end point and if nec
    x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim = trim_x_y_values_to_start_end_point(
            x_start_index, x_end_index, x_values, x_values_indizes, y_values, y_values_indizes)

    if min(y_values_indizes) < 0 or  max(y_values_indizes) > grid_ressolution_int-1:
      x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim = trim_x_y_values_to_geometry(
                grid_ressolution_int, x_slope, x_values_trim, x_values_indizes_trim, y_values_trim, y_values_indizes_trim)

    # z Values from grid
    z_values_new_bending_direction = []
    for i in range(len(x_values_indizes_trim)):
        z_values_new_bending_direction.append(
            z_grid_values_linear[x_values_indizes_trim[i], y_values_indizes_trim[i]])
    z_values_new_bending_direction = np.asarray(z_values_new_bending_direction, dtype=np.float32)
    # x, y and z stacked together to 3D-Points
    new_bending_direction_points_current_KOS = np.vstack(
        (x_values_trim, y_values_trim, z_values_new_bending_direction)).T
    return new_bending_direction_points_current_KOS, y_intercept, x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim

def new_bending_points_in_global_and_tilted_KOS(start_point_drawn, end_point_drawn, y_intercept,
                                                new_bending_direction_points_current_KOS,
                                                trendline_new_direction_current_KOS):
    # new zero
    new_zero = np.asarray((0, y_intercept, 0), dtype=np.float32)

    start_end_point_current_KOS = np.vstack((start_point_drawn, end_point_drawn))

    start_end_point_drawn_global_KOS = translation_Points_from_OLD_to_trendline_KOS(
        start_end_point_current_KOS,
        trendline_global_KOS,
        center_point_of_cloud_weighted, False, True)
    trendline_new_direction_global_KOS = translation_Points_from_OLD_to_trendline_KOS(
        trendline_new_direction_current_KOS,
        trendline_global_KOS,
        (center_point_of_cloud_weighted-center_point_of_cloud_weighted), False, True)

    new_bending_direction_points_global_KOS = translation_Points_from_OLD_to_trendline_KOS(
        new_bending_direction_points_current_KOS,
        trendline_global_KOS,
        center_point_of_cloud_weighted, False, True)

    new_bending_direction_points_tilted_KOS = translation_Points_from_OLD_to_trendline_KOS(
        new_bending_direction_points_current_KOS, trendline_new_direction_current_KOS, new_zero)

    return new_bending_direction_points_tilted_KOS, new_bending_direction_points_global_KOS, start_end_point_drawn_global_KOS, trendline_new_direction_global_KOS

def new_bending_points_tilted_to_global(y_intercept,
                                                new_bending_direction_points_current_KOS,
                                                trendline_new_direction_current_KOS):
    # new zero
    new_zero = np.asarray((0, y_intercept, 0), dtype=np.float32)

    new_bending_direction_points_trendline_KOS = translation_Points_from_OLD_to_trendline_KOS(
        new_bending_direction_points_current_KOS, trendline_new_direction_current_KOS, new_zero, False, True)

    new_bending_direction_points_global_KOS = translation_Points_from_OLD_to_trendline_KOS(
        new_bending_direction_points_trendline_KOS,
        trendline_global_KOS,
        center_point_of_cloud_weighted, False, True)

    return new_bending_direction_points_global_KOS

def trim_x_y_values_to_geometry(grid_ressolution_int, x_slope, x_values, x_values_indizes, y_values, y_values_indizes):
    # -1 default wert
    y_start_index = -1
    y_end_index = -1
    for k in range(len(x_values)):

        if x_slope >= 0:
            if (y_values_indizes[k] >= 0) & (y_start_index < 0): y_start_index = k
            if (y_values_indizes[k] > grid_ressolution_int) & (y_end_index < 0): y_end_index = k - 1

        if x_slope < 0:
            if (y_values_indizes[k] < grid_ressolution_int) & (y_start_index < 0): y_start_index = k
            if (y_values_indizes[k] <= 0) & (y_end_index < 0): y_end_index = k - 1
    if y_end_index <= 0: y_end_index = grid_ressolution_int - 1
    # Trim indizes and coordinates to grid size
    x_values_indizes_trim = x_values_indizes[y_start_index:y_end_index]
    y_values_indizes_trim = y_values_indizes[y_start_index:y_end_index]
    x_values_trim = x_values[y_start_index:y_end_index]
    y_values_trim = y_values[y_start_index:y_end_index]
    return x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim

def trim_x_y_values_to_start_end_point(x_start_index, x_end_index, x_values, x_values_indizes, y_values,
                                       y_values_indizes):
   # Trim indizes and coordinates to grid size
    x_values_indizes_trim = x_values_indizes[x_start_index:x_end_index]
    y_values_indizes_trim = y_values_indizes[x_start_index:x_end_index]
    x_values_trim = x_values[x_start_index:x_end_index]
    y_values_trim = y_values[x_start_index:x_end_index]
    return x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim

# Functionen in Interpolate start_geometrie
def show_interpolation_and_draw_start_end_points(max_x, max_y, min_x, min_y, z_grid_values_linear):
    figure = pyplot.figure()  # Comment_DB: create a new figure
    plt.imshow(z_grid_values_linear.T, extent=(min_x, max_x, min_y, max_y), origin='lower')
    # plt.plot(x_values_trim, y_values, 'bo', linewidth=2.0, label='Schnitt')
    plt.title('Linear')
    global xdata, ydata
    xdata = []
    ydata = []

    def onclick(event):
        print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              (event.button, event.x, event.y, event.xdata, event.ydata))
        plt.plot(event.xdata, event.ydata, ',')
        figure.canvas.draw()
        global xdata, ydata
        xdata.append(event.xdata)
        ydata.append(event.ydata)

    figure.canvas.mpl_connect('button_press_event', onclick)
    global continue_bool
    continue_bool = False

    def handle_close(event):
        print('Closed Figure!')
        global continue_bool
        continue_bool = True

    figure.canvas.mpl_connect('close_event', handle_close)
    plt.show()
    while continue_bool is False:
        pyplot.pause(2)
def interpolate_geometrie(grid_resolution):
    # Comment_DKu_Wenzel: Interpolation mit Centerpoints teils ungenauer
    tri_centerpoints_rotatet_and_translated = translation_Points_from_OLD_to_trendline_KOS(tri_centerpoints,
                                                                                           trendline_global_KOS,
                                                                                           center_point_of_cloud_weighted)
    tri_corner__points_rotatet_and_translated = translation_Points_from_OLD_to_trendline_KOS(tri_corner_points,
                                                                                             trendline_global_KOS,
                                                                                             center_point_of_cloud_weighted,
                                                                                              True)

    global tri_corner__points_rotatet_and_translated_reverse
    tri_corner__points_rotatet_and_translated_reverse = translation_Points_from_OLD_to_trendline_KOS(tri_corner__points_rotatet_and_translated,
                                                                                             trendline_global_KOS,
                                                                                             center_point_of_cloud_weighted,
                                                                                              False, True)


    interpolate_with_corner_and_centerpoints = False
    if interpolate_with_corner_and_centerpoints:
        points = np.concatenate((tri_centerpoints_rotatet_and_translated, tri_corner__points_rotatet_and_translated))
    else:
        points = tri_corner__points_rotatet_and_translated

    corner_points_X_Y_TrendlineKOS = points[:, 0:2]
    corner_points_Z_TrendlineKOS = points[:, 2]
    max_x = max(corner_points_X_Y_TrendlineKOS[:, 0])
    min_x = min(corner_points_X_Y_TrendlineKOS[:, 0])
    max_y = max(corner_points_X_Y_TrendlineKOS[:, 1])
    min_y = min(corner_points_X_Y_TrendlineKOS[:, 1])
    grid_x, grid_y = np.mgrid[min_x:max_x:grid_resolution, min_y:max_y:grid_resolution]
    corner_points_X_Y_TrendlineKOS = np.asarray(corner_points_X_Y_TrendlineKOS, dtype=np.float32)
    corner_points_Z_TrendlineKOS = np.asarray(corner_points_Z_TrendlineKOS, dtype=np.float32)
    # Interpolating  the Surface Geometry
    z_grid_values_linear = griddata(corner_points_X_Y_TrendlineKOS, corner_points_Z_TrendlineKOS, (grid_x, grid_y),
                                    method='linear')
    return grid_x, grid_y, max_x, max_y, min_x, min_y, z_grid_values_linear
def translation_Points_from_OLD_to_trendline_KOS(points_in_old_KOS, new_trendline_axis_in_old_KOS, new_zero_point_in_old_KOS, add_start_point_to_pointset = False, reverse = False):
    # Gesamtidee: Erst wird trendx zu (1,0,0) rotiert, anschließend werden die Punkte in den center point weighted(cpw) verschoben

    # Basic Coordinate System
    x_axis = np.asarray((1, 0, 0), dtype=np.float32)
    y_axis = np.asarray((0, 1, 0), dtype=np.float32)
    z_axis = np.asarray((0, 0, 1), dtype=np.float32)

    # Rotationswinkel
    anglez, angley = calc_angle_coordinate_rotation_x_trendline(x_axis, z_axis, new_trendline_axis_in_old_KOS)

    if reverse:
        anglez=-anglez
        angley=-angley


    points_in_trendline_KOS = []

    for i in range(len(points_in_old_KOS[:, 0])):
        tri_corner__points_rotatet_i = rotate_point_around_z_and_y_axis_with_given_angle(anglez, angley, z_axis, y_axis,
                                                                                         points_in_old_KOS[i],new_trendline_axis_in_old_KOS,reverse)
        points_in_trendline_KOS.append(tri_corner__points_rotatet_i)

    # Comment_DKu_Wenzel: Startpunkt ist von einem Eckpunkt auf die Trendlinien-x-z-Ebene-Projeziert
    if add_start_point_to_pointset :
        startpoint_project_to_trendline_plan_rotated = rotate_point_around_z_and_y_axis_with_given_angle(anglez, angley,
                                                                                                     z_axis, y_axis,
                                                                                                     startpoint_project_to_trendline_plan,new_trendline_axis_in_old_KOS,reverse)
        points_in_trendline_KOS.append(startpoint_project_to_trendline_plan_rotated)

    points_in_trendline_KOS = np.asarray(points_in_trendline_KOS, dtype=np.float32)

    # rotation of translation vector/new zero
    new_zero_point_in_old_KOS_rotated = rotate_point_around_z_and_y_axis_with_given_angle(anglez, angley, z_axis, y_axis,new_zero_point_in_old_KOS,new_trendline_axis_in_old_KOS,reverse)

    # At the reverse case we have to shift in the old unrotated KOS
    if reverse: new_zero_point_in_old_KOS_rotated = -(new_zero_point_in_old_KOS)


    #  shifting rotated points to new zero
    points_in_trendline_KOS[:, 0] = np.subtract(points_in_trendline_KOS[:, 0], new_zero_point_in_old_KOS_rotated[0])
    points_in_trendline_KOS[:, 1] = np.subtract(points_in_trendline_KOS[:, 1], new_zero_point_in_old_KOS_rotated[1])
    points_in_trendline_KOS[:, 2] = np.subtract(points_in_trendline_KOS[:, 2], new_zero_point_in_old_KOS_rotated[2])

    return points_in_trendline_KOS
def rotate_point_around_z_and_y_axis_with_given_angle(angle1, angle2, axis1, axis2, point_to_rotate, new_trendline_axis_in_old_KOS,reverse):
    if reverse:
        point_rotated_around_axis2 = Quaternion(axis=axis2, angle=angle2).rotate(point_to_rotate)
        point_rotated_around_1_and_2 = Quaternion(axis=axis1, angle=-angle1).rotate(point_rotated_around_axis2)
    else:
        rotated_point_around_axis1 = Quaternion(axis=axis1, angle=-angle1).rotate(
            point_to_rotate)
        point_rotated_around_1_and_2 = Quaternion(axis=axis2, angle=angle2).rotate(rotated_point_around_axis1)

    # If case
        # The first trendline created is defined x,z,y. If the normal of the plain goes in to negativ direction, rotate around x = 180°
        # For all following lokal KOS this statement should be false
    #if new_trendline_axis_in_old_KOS[1][2] < 0: # todo: demonstrator tool
    #    point_rotated_around_1_and_2 = Quaternion(axis=(1,0,0), angle=np.pi).rotate(point_rotated_around_1_and_2)
    return point_rotated_around_1_and_2
def calc_angle_coordinate_rotation_x_trendline(x_axis, z_axis, new_trendline_axis_in_old_KOS):
    new_x_trendline_projected_to_x_y = project_pointtoplane((new_trendline_axis_in_old_KOS[0][:]), z_axis, np.zeros(3))
    new_x_trendline_projected_to_x_y = 1 / np.linalg.norm(
        new_x_trendline_projected_to_x_y) * new_x_trendline_projected_to_x_y

    anglez = math.acos(np.dot(x_axis, new_x_trendline_projected_to_x_y))
    # Wenn y negativ, in die x_Rotation in die andere Richtung korrigieren
    if new_trendline_axis_in_old_KOS[0][1] <= -0.01: anglez = -anglez


    rotated_x_trend_around_z = Quaternion(axis=z_axis, angle=-anglez).rotate(new_trendline_axis_in_old_KOS[0][:])
    rotated_x_trend_around_z = 1 / np.linalg.norm(rotated_x_trend_around_z) * rotated_x_trend_around_z

    angley = math.acos(np.dot(x_axis, rotated_x_trend_around_z) / (
            np.linalg.norm(x_axis) * np.linalg.norm(rotated_x_trend_around_z)))

    # Wenn z negativ, in die x_Rotation in die andere Richtung korrigieren
    if rotated_x_trend_around_z[2] <= -0.01: angley = -angley

    return anglez, angley


def show_startstrip(bestPatch_patternpoints,patch_start,patch_end):
    ###2D-xy-PLOT
    plt.plot(x_y_points_filled_up[:, 0], x_y_points_filled_up[:, 1], 'bo', linewidth=2.0, label='ohne Glättung')  # äquidistante Punkte

    plt.plot(x_y_points_filled_up[:, 0], y_smooth, color='r', linewidth = 3, label ='Savitzky-Golay')  # SavGol-Glättung

    plt.plot(bend_pts_xz_local[:, 0], bend_pts_xz_local[:, 1], color='green', linewidth=3.0, label='lineare Angleichung')  # Streckenweise linear (nur Eckpunkte)
    plt.axis([x_list[0] - 50, x_list[-1] + 50, -1 * max(y_list)-50, max(y_list)+50])

    plt.xlabel('[ mm ]')
    plt.ylabel('[ mm ]')
    plt.legend()

    ###############3D-PLOTTING################ new_bending_direction_points_previous_KOS

    figure = pyplot.figure() #Comment_DB: 3D plot of objective shape
    axes = mplot3d.Axes3D(figure)

    patch_visual = mplot3d.art3d.Poly3DCollection(triangle_vectors_of_stl, linewidths=1, alpha=0.5, edgecolor=[1, 1, 1], label ='Geometriebereich') #Comment_DB: added edgecolor to make the edges visible

    axes.scatter([999999990],[9999999900],[9999999900],linewidths=0.0001, alpha = 0.5, label = "Geometriebereich") #Comment_DB: label in legend

    axes.scatter(center_point_of_cloud_weighted[0],center_point_of_cloud_weighted[1],center_point_of_cloud_weighted[2],c='g')
    axes.scatter(bend_pts_xy_global[:,0], bend_pts_xy_global[:,1],
                 bend_pts_xy_global[:,2], c='g')
    axes.scatter(bend_pts_xz_global_left[:, 0], bend_pts_xz_global_left[:, 1],
                 bend_pts_xz_global_left[:, 2], c='g')
    axes.scatter(bend_pts_xy_global_right[:, 0], bend_pts_xy_global_right[:, 1],
                 bend_pts_xy_global_right[:, 2], c='g')

    axes.plot([bend_pts_xy_global_right[1, 0], bend_pts_xz_global_left[1, 0]], [bend_pts_xy_global_right[1, 1], bend_pts_xz_global_left[1, 1]], [bend_pts_xy_global_right[1, 2], bend_pts_xz_global_left[1, 2]], label='Trendlinie', c='red')  # Comment_DB: *pc_axes is *args, and .T is np.transpose
    pc_axes=np.asarray(trendline)
    
    # TRENDLINE
    axes.plot(*pc_axes.T,label='Trendlinie', c='red') #Comment_DB: *pc_axes is *args, and .T is np.transpose
    axes.scatter(patch_start[0], patch_start[1], patch_start[2], c="black")
    axes.scatter(patch_end[0],patch_end[1],patch_end[2],c='black')

    axes.scatter(start_end_point_drawn_global_KOS[:,0], start_end_point_drawn_global_KOS[:, 1],
                 start_end_point_drawn_global_KOS[:, 2], c='black')

    axes.scatter(start_end_point_drawn_global_KOS_left[:, 0], start_end_point_drawn_global_KOS_left[:, 1],
                 start_end_point_drawn_global_KOS_left[:, 2], c='black')

    axes.scatter(start_end_point_drawn_global_KOS_right[:, 0], start_end_point_drawn_global_KOS_right[:, 1],
                 start_end_point_drawn_global_KOS_right[:, 2], c='black')


    # von PCC gemittelter Normalenvektor
    x3, y3, z3 = [center_point_of_cloud_weighted[0], center_point_of_cloud_weighted[0] + 500 * avg_tri_normal_weighted[0]], \
                 [center_point_of_cloud_weighted[1],center_point_of_cloud_weighted[1] + 500 * avg_tri_normal_weighted[1]], \
                 [center_point_of_cloud_weighted[2], center_point_of_cloud_weighted[2] + 500 * avg_tri_normal_weighted[2]]

    plt.plot(x3,y3,z3,marker='o',c='green')
    """
    x333, y333, z333 = [center_point_of_cloud_weighted[0],500 *trendline_x_axis[0]],[center_point_of_cloud_weighted[1],500 *trendline_x_axis[1]],[center_point_of_cloud_weighted[2],500 *trendline_x_axis[2]]

    x332, y332, z332 = [center_point_of_cloud_weighted[0],500 *trendline_z_axis[0]],[center_point_of_cloud_weighted[1],500 *trendline_z_axis[1]],[center_point_of_cloud_weighted[2],500 *trendline_z_axis[2]]

    x331, y331, z331 = [center_point_of_cloud_weighted[0],500 *trendline_y_axis[0]],[center_point_of_cloud_weighted[1],500 *trendline_y_axis[1]],[center_point_of_cloud_weighted[2],500 *trendline_y_axis[2]]

    plt.plot(x333, y333, z333, marker='o', c='blue')
    plt.plot(x332, y332, z332, marker='o', c='blue')
    plt.plot(x331, y331, z331, marker='o', c='blue')
    """

    # Plot new tildet bending Points
    plt.plot(list(tri_corner__points_rotatet_and_translated_reverse[:, 0]), list(tri_corner__points_rotatet_and_translated_reverse[:, 1]),
             list(tri_corner__points_rotatet_and_translated_reverse[:, 2]), marker='o', c='blue')

    """
    plt.plot(list(new_bending_direction_points_global_KOS[:, 0]), list(new_bending_direction_points_global_KOS[:, 1]), list(new_bending_direction_points_global_KOS[:, 2]), marker='o', c='blue')
    plt.plot(list(new_bending_direction_points_global_KOS_left[:, 0]), list(new_bending_direction_points_global_KOS_left[:, 1]),
             list(new_bending_direction_points_global_KOS_left[:, 2]), marker='o', c='blue')
    plt.plot(list(new_bending_direction_points_global_KOS_right[:, 0]),
             list(new_bending_direction_points_global_KOS_right[:, 1]),
             list(new_bending_direction_points_global_KOS_right[:, 2]), marker='o', c='blue')
    #axes.scatter(tri_corner_points_projected_to_trendline_plane[:, 0], tri_corner_points_projected_to_trendline_plane[:, 1], tri_corner_points_projected_to_trendline_plane[:, 2], c='r')
    #axes.scatter(tri_corner_points_projected_to_trendline[:, 0], tri_corner_points_projected_to_trendline[:, 1], tri_corner_points_projected_to_trendline[:, 2], c='g')
    """
    
    for i in range(len(bestPatch_patternpoints) - 2):
        verts = [list(
            zip([bestPatch_patternpoints[i][0], bestPatch_patternpoints[i + 1][0], bestPatch_patternpoints[i + 2][0]], \
                [bestPatch_patternpoints[i][1], bestPatch_patternpoints[i + 1][1], bestPatch_patternpoints[i + 2][1]], \
                [bestPatch_patternpoints[i][2], bestPatch_patternpoints[i + 1][2],
                 bestPatch_patternpoints[i + 2][2]]))]  # Comment_DB: DARK BLUE LoP PATCH
        axes.add_collection3d(Poly3DCollection(verts), zs='z')  # Comment_DB: INSERT LoP PATCH IN GRAPH
        # patch_meshpoints.append(verts) #Comment_DB: is not used
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
