import numpy as np
import math
import scipy
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
max_points_in_pc = 72000



####MAIN
def startparam(input_file, max_distance, width_for_edge_detection, grid_resolution,interpolate_with_corner_and_centerpoints):
    print(
        "Preprocessing läuft... Bei hochaufgelösten stl-Dateien und großen Flächengrößenunterschieden kann es zu "
        "längeren Berechnungszeiten kommen")

    # Read in the file:
    calc_trendline_of_geometry_from_stl_file(input_file) # saved in global var!

    # Interpolation:
    grid_resolution_j = grid_resolution*1j
    grid_x, max_x, max_y, min_x, min_y, z_grid_values_linear = interpolate_geometrie_in_trendline_KOS_and_draw_start_end_point( # start end saved in global var!
         grid_resolution_j, interpolate_with_corner_and_centerpoints)

    # Calculation of the bending parameters:
    calc_bending_parameters(grid_resolution, grid_x, max_distance, width_for_edge_detection, max_x, max_y, min_x, min_y, z_grid_values_linear) # saved in global var!

    ###### Startparameter extrahieren #####
    startpoint_on_surface = bend_pts_xyz_global_3D[0][:]
    endpoint_on_surface = bend_pts_xyz_global_3D[-1][:]

    Start_r_3d_atstart = norm_vector(x_direction_list_global_KOS[0])
    Start_n_3d_atstart = norm_vector(normal_patch_global_KOS[0])

    # l-list for Chromo
    l_list = np.asarray(length_list)
    # Beta_liste aus Interpolator
    beta_list = np.asarray(beta_angle_between_planes_list)#Comment_DB: THIS IS THE SAME AS BETA_LIST FOR PREPROCESSED CHROMO IN LoP, EXCEPT HERE IT IS IN DEGREES
    # Alpha_liste aus Interpolator
    alpha_list = np.asarray(alpha_angle_list)
    # L_aim #Comment_DKu_Wenzel: Kein Savgol mehr! Nimmt Kurve von Interpolierter Fläche!!
    L_aim = calc_L_aim(new_bending_direction_points_tilted_KOS_stacked)

    start_parameter = [l_list, L_aim, beta_list, alpha_list, startpoint_on_surface,
                       endpoint_on_surface, Start_r_3d_atstart, Start_n_3d_atstart]
    return start_parameter

#######################################################################################################################
# Load stl file and analyse the geometry. Calc centerpoint and trendline of the geometry
def calc_trendline_of_geometry_from_stl_file(input_file):
    patch_vectors_of_stl_input = mesh.Mesh.from_file(input_file)  # Comment_DB: stl mesh

    # Triangel vectors needed for the visualization of the patch.
    global triangle_vectors_of_stl,num_of_triangles
    triangle_vectors_of_stl = patch_vectors_of_stl_input.vectors  # Comment_DB: triangle edges (wireframe)
    num_of_triangles = len(triangle_vectors_of_stl)

    # Calc areas of triangel for compare/weight
    tri_areas = calc_tri_areas(triangle_vectors_of_stl)


    global center_point_of_cloud_weighted, trendline_global_KOS, avg_tri_normal_weighted
    avg_tri_normal_weighted = calc_avg_tri_norm_weighted_by_area(tri_areas, patch_vectors_of_stl_input)

    # Calculate approximate center of geometry
    tri_centerpoints = calc_tri_centerpoints(triangle_vectors_of_stl)  # Comment_DKu_Wenzel: basically the unweighted point cloud
    point_cloud_tri_centerpoints_weighted = calc_patch_pointcloud_weighted_by_area(tri_areas,tri_centerpoints)
    center_point_of_cloud_weighted = point_cloud_tri_centerpoints_weighted.mean(axis=0) # Mean of x,y,z-Values

    # SVD for approximate orientation of geometry (trendline KOS)
    trendline_global_KOS = calc_trendline_global(center_point_of_cloud_weighted, point_cloud_tri_centerpoints_weighted)
# Functions in calc trendline
def calc_trendline_global(center_point_of_cloud_weighted, point_cloud_tri_centerpoints_weighted):
    trendline_x_axis, trendline_y_axis, trendline_z_axis = calc_trendline_axis_with_svd(
        point_cloud_tri_centerpoints_weighted, center_point_of_cloud_weighted)
    # Wenn trendline x-Achse in negative x-Richtung definiert wurde von svd
    if trendline_x_axis[0] < 0:  # Rotation von 180° um y-Achse
        trendline_x_axis = -(trendline_x_axis)
        trendline_z_axis = -(trendline_z_axis)
    trendline_global_KOS = np.vstack((trendline_x_axis, trendline_y_axis, trendline_z_axis))
    return trendline_global_KOS
def calc_tri_normals_from_stl(stl_normals,triangle_vectors_of_stl):
    normals=[]
    #We generate our own normals with the id_list. Notice that the triangle order of stl_normals and the vertices
    #(triangles) is not the same

    #Comment_DKu_Wenzel: We are just taking the length of tri_id_sorted. So we iterate normal and not the different id
    for i in range(num_of_triangles):

        v1= triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][1]
        v2= triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][2]
        n=np.cross(v1,v2)
        n=norm_vector(n)
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
def calc_tri_corner_points(triangle_vectors_of_stl):
    tri_corner_points = []
    for i in range(num_of_triangles):
        triangle = triangle_vectors_of_stl[i]
        for j in range(3):
            tri_corner_points.append(triangle[j])

    # Doppelte Ecken löschen
    tri_corner_points = np.unique(tri_corner_points, axis=0)

    return tri_corner_points
def calc_tri_areas(triangle_vectors_of_stl):
    #Berechnung der Dreiecksflächen und Speichern in einer Liste
    tri_surface_area = []
    for i in range(num_of_triangles):
        tri_surface_area.append(0.5 * (
            np.linalg.norm((triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][1]) - (triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][2]))))
    tri_surface_area = np.asarray(tri_surface_area)
    return tri_surface_area
def calc_avg_tri_norm_weighted_by_area(tri_areas,patch_vectors_of_stl_input):
    stl_normals = (patch_vectors_of_stl_input.normals)
    tri_normals = calc_tri_normals_from_stl(stl_normals, triangle_vectors_of_stl)

    weighted_norms = []
    for i in range(num_of_triangles):
        weighted_norms.append((tri_areas[i] / sum(tri_areas)) * tri_normals[i])
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
def calc_trendline_axis_with_svd(patch_pc_weighted, center_point_of_cloud_weighted):
    # Do Principal Component Analysis(PCA) on the mean-centered data. AKA SVD
    # The first principal component contains [uu, dd, vv] , where vv[0] is the direction
    first_principal_components_pc_weighted = scipy.linalg.svd(patch_pc_weighted - center_point_of_cloud_weighted)
    # Definition der Hauptachsen
    trendline_x_axis = first_principal_components_pc_weighted[2][0] # first_principal_components_pc_weighted[2][0]: is direction of trendline
    trendline_x_axis = norm_vector(trendline_x_axis)
    # avg_tri_norm ist nicht senkrecht zur x-Achse
    # von pcc + avg_tri_norm und zurück auf x-Achse projizieren
    trendline_avg_norm_point = center_point_of_cloud_weighted + np.dot(avg_tri_normal_weighted,
                                                                       trendline_x_axis) / np.dot(trendline_x_axis,
                                                                                                  trendline_x_axis) * trendline_x_axis
    # y-Achse ist verbindung von pcc+avg_tri_norm mit dem projizierten Punkt
    trendline_z_axis = (center_point_of_cloud_weighted + avg_tri_normal_weighted) - trendline_avg_norm_point
    trendline_z_axis = norm_vector(trendline_z_axis)
    trendline_y_axis = np.cross(trendline_z_axis, trendline_x_axis)
    return trendline_x_axis, trendline_y_axis, trendline_z_axis
def find_nearest(array, value):
# find the index of the closest value in an array to a given value
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
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
def norm_vector(vector):
    vector = 1 / np.linalg.norm(vector) * vector
    return vector


#######################################################################################################################
# Interpolte geometry in rotated and translated trendline KOS
def interpolate_geometrie_in_trendline_KOS_and_draw_start_end_point(grid_resolution, interpolate_with_corner_and_centerpoints):
    grid_x, max_x, max_y, min_x, min_y, z_grid_values_linear = interpolate_geometrie(grid_resolution,interpolate_with_corner_and_centerpoints)

    show_interpolation_and_draw_start_end_points(max_x, max_y, min_x, min_y, z_grid_values_linear) # In this function x_data_drawn and y_data_drawn are calculated

    return grid_x, max_x, max_y, min_x, min_y, z_grid_values_linear
# Functions in Interpolate start_geometrie
def interpolate_geometrie(grid_resolution,interpolate_with_corner_and_centerpoints):
    # Comment_DKu_Wenzel: Interpolation mit Centerpoints teils ungenauer

    # Taking corner points of the triangles as interpolation points.
    # First thing we do is rotating the cornerpoints into the trendline KOS
    points = calc_points_in_trendline_KOS_for_interpolation(interpolate_with_corner_and_centerpoints)

    points_x_y_trendline_KOS = points[:, 0:2]
    points_z_trendline_KOS = points[:, 2]

    # Creat grid
    max_x = max(points_x_y_trendline_KOS[:, 0])
    min_x = min(points_x_y_trendline_KOS[:, 0])
    max_y = max(points_x_y_trendline_KOS[:, 1])
    min_y = min(points_x_y_trendline_KOS[:, 1])
    grid_x, grid_y = np.mgrid[min_x:max_x:grid_resolution, min_y:max_y:grid_resolution]


    # Interpolating  the Surface Geometry
    z_grid_values_linear = griddata(np.asarray(points_x_y_trendline_KOS, dtype=np.float32), np.asarray(points_z_trendline_KOS, dtype=np.float32), (grid_x, grid_y),
                                    method='linear')
    return grid_x, max_x, max_y, min_x, min_y, z_grid_values_linear
def calc_points_in_trendline_KOS_for_interpolation(interpolate_with_corner_and_centerpoints):
    tri_corner_points = calc_tri_corner_points(triangle_vectors_of_stl)
    tri_corner__points_rotatet_and_translated = translate_and_rotate_points_from_OLD_to_trendline_KOS(tri_corner_points,
                                                                                                      trendline_global_KOS,
                                                                                                      center_point_of_cloud_weighted)
    if interpolate_with_corner_and_centerpoints:
        tri_centerpoints = calc_tri_centerpoints(
            triangle_vectors_of_stl)
        tri_centerpoints_rotatet_and_translated = translate_and_rotate_points_from_OLD_to_trendline_KOS(tri_centerpoints,
                                                                                                        trendline_global_KOS,
                                                                                                        center_point_of_cloud_weighted)
        points = np.concatenate((tri_centerpoints_rotatet_and_translated, tri_corner__points_rotatet_and_translated))
    else:
        points = tri_corner__points_rotatet_and_translated
    return points
def show_interpolation_and_draw_start_end_points(max_x, max_y, min_x, min_y, z_grid_values_linear):
    # Show interpolation
    figure = pyplot.figure()  # Comment_DB: create a new figure
    plt.imshow(z_grid_values_linear.T, extent=(min_x, max_x, min_y, max_y), origin='lower')
    # plt.plot(x_values_trim, y_values, 'bo', linewidth=2.0, label='Schnitt')
    plt.title('Please select Start- and Endpoint')

    # Draw Start and Endpoint
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


#######################################################################################################################
# Translation and Rotation from Points to new trendline_axis_KOS
def translate_and_rotate_points_from_OLD_to_trendline_KOS(points_in_old_KOS, new_trendline_axis_in_old_KOS,
                                                          new_zero_point_in_old_KOS, reverse=False):
    # Gesamtidee: Erst wird die neue trendline_axis_x zu (1,0,0) rotiert, anschließend werden die Punkte in den center point weighted(cpw) verschoben

    # Basic Coordinate System
    x_axis = np.asarray((1, 0, 0), dtype=np.float32)
    y_axis = np.asarray((0, 1, 0), dtype=np.float32)
    z_axis = np.asarray((0, 0, 1), dtype=np.float32)

    # Rotationswinkel
    anglez, angley = calc_angle_for_coordinate_rotation_x_trendline(x_axis, z_axis, new_trendline_axis_in_old_KOS) # x-Axis now: (1,0,0)
    anglex = correction_angle_x_axis_for_y_z_orientation(angley, anglez, new_trendline_axis_in_old_KOS, reverse, y_axis,z_axis) # for y(0,1,0) and z(0,0,1)
    if reverse:
        anglez=-anglez
        angley=-angley
        anglex=-anglex

    # Rotation of points
    points_in_trendline_KOS = rotated_points(anglex, angley, anglez, points_in_old_KOS, reverse, x_axis, y_axis, z_axis)

    # Rotation of translation vector/new zero
    if reverse: translation_vector = -(new_zero_point_in_old_KOS)  # At the reverse case we have to shift in the old unrotated KOS
    else: translation_vector = rotate_point_around_z_y_and_x_axis_with_given_angle(anglez, angley, anglex, z_axis,
                                                                                            y_axis, x_axis,
                                                                                            new_zero_point_in_old_KOS,
                                                                                            reverse)
    #  translate rotated points to new zero
    translate_points_to_new_zero(translation_vector, points_in_trendline_KOS)

    return points_in_trendline_KOS
def calc_angle_for_coordinate_rotation_x_trendline(x_axis, z_axis, new_trendline_axis_in_old_KOS):

    # Calculation of angel z
    new_x_trendline_projected_to_x_y = project_pointtoplane((new_trendline_axis_in_old_KOS[0][:]), z_axis, np.zeros(3))
    new_x_trendline_projected_to_x_y = norm_vector(new_x_trendline_projected_to_x_y)

    anglez = math.acos(np.dot(x_axis, new_x_trendline_projected_to_x_y))
    # Wenn y negativ, in die x_Rotation in die andere Richtung korrigieren
    if new_trendline_axis_in_old_KOS[0][1] <= -0.0001: anglez = -anglez


    # Calculation of angel y
    rotated_x_trend_around_z = Quaternion(axis=z_axis, angle=-anglez).rotate(new_trendline_axis_in_old_KOS[0][:])
    rotated_x_trend_around_z = norm_vector(rotated_x_trend_around_z)

    angley = math.acos(np.dot(x_axis, rotated_x_trend_around_z))
    # Wenn z negativ, in die x_Rotation in die andere Richtung korrigieren
    if rotated_x_trend_around_z[2] <= -0.0001: angley = -angley

    return anglez, angley
def correction_angle_x_axis_for_y_z_orientation(angley, anglez, new_trendline_axis_in_old_KOS, reverse, y_axis, z_axis):
    # Rotate trendline_axis around y and z
    new_trendline_axis_in_trendline_KOS = []
    for i in range(len(new_trendline_axis_in_old_KOS[:, 0])):
        new_trendline_axis_points_rotatet_i = rotate_point_around_z_y_and_x_axis_with_given_angle(anglez, angley, 0, z_axis,
                                                                                                  y_axis, y_axis,
                                                                                                  new_trendline_axis_in_old_KOS[
                                                                                                    i], reverse)
        new_trendline_axis_in_trendline_KOS.append(new_trendline_axis_points_rotatet_i)

    # Check orientations of rotated y and z
    if new_trendline_axis_in_trendline_KOS[1][1] < 0.9:  #If new_trendline_axis_in_trendline_KOS[1][1] > 0.9 -> then y = [0,1,0], no rotation needed
        if new_trendline_axis_in_trendline_KOS[1][1] < -0.9:
            anglex = math.pi
        elif new_trendline_axis_in_trendline_KOS[1][2] > 0.9:
            anglex = -math.pi / 2
        else:  # new_trendline_axis_in_trendline_KOS[1][2] < -0.9:
            anglex = math.pi / 2
    else:
        anglex = 0

    return anglex
def rotate_point_around_z_y_and_x_axis_with_given_angle(angle1, angle2, angle3, axis1, axis2, axis3, point_to_rotate, reverse):
    if reverse:
        point_rotated_around_axis3 = Quaternion(axis=axis3, angle=angle3).rotate(point_to_rotate)
        point_rotated_around_3_and_2 = Quaternion(axis=axis2, angle=angle2).rotate(point_rotated_around_axis3)
        point_rotated_around_1_2and_3 = Quaternion(axis=axis1, angle=-angle1).rotate(point_rotated_around_3_and_2)
    else:
        rotated_point_around_axis1 = Quaternion(axis=axis1, angle=-angle1).rotate(point_to_rotate)
        point_rotated_around_1_and_2 = Quaternion(axis=axis2, angle=angle2).rotate(rotated_point_around_axis1)
        point_rotated_around_1_2and_3 = Quaternion(axis=axis3, angle=angle3).rotate(point_rotated_around_1_and_2)

    return point_rotated_around_1_2and_3
def rotated_points(anglex, angley, anglez, points_in_old_KOS, reverse, x_axis, y_axis, z_axis):
    points_in_trendline_KOS = []
    for i in range(len(points_in_old_KOS[:, 0])):
        old_points_rotatet_i = rotate_point_around_z_y_and_x_axis_with_given_angle(anglez, angley, anglex, z_axis,
                                                                                   y_axis, x_axis,
                                                                                   points_in_old_KOS[i], reverse)
        points_in_trendline_KOS.append(old_points_rotatet_i)
    points_in_trendline_KOS = np.asarray(points_in_trendline_KOS, dtype=np.float32)
    return points_in_trendline_KOS
def translate_points_to_new_zero(new_zero_point_in_old_KOS_rotated, points_in_trendline_KOS):
    points_in_trendline_KOS[:, 0] = np.subtract(points_in_trendline_KOS[:, 0], new_zero_point_in_old_KOS_rotated[0])
    points_in_trendline_KOS[:, 1] = np.subtract(points_in_trendline_KOS[:, 1], new_zero_point_in_old_KOS_rotated[1])
    points_in_trendline_KOS[:, 2] = np.subtract(points_in_trendline_KOS[:, 2], new_zero_point_in_old_KOS_rotated[2])


#######################################################################################################################
# Calculation of the bending parameters
def calc_bending_parameters(grid_ressolution_int, grid_x, max_distance, width_for_edge_detection, max_x, max_y, min_x, min_y,
                            z_grid_values_linear):

    y_0_grid_point_index = np.asarray(np.round(max_y / (max_y - min_y) * grid_ressolution_int), dtype=np.int32)
    x_0_grid_point_index = np.asarray(np.round(max_x / (max_x - min_x) * grid_ressolution_int), dtype=np.int32)


    initialize_global_lists_of_bending_and_plot_parameter()

    # Start calculating bendingpoints
    calculate_iteratively_the_tilted_bendingpoints(alpha_angle_list, grid_ressolution_int, grid_x, max_distance, width_for_edge_detection,max_x,
                                                   max_y, min_x, min_y, x_0_grid_point_index, x_values_trim_stacked,
                                                   xdata_list, y_0_grid_point_index, y_values_trim_stacked, ydata_list,
                                                   z_grid_values_linear)
    # Show results
    show_results_2D_Plots_and_Colormap(max_x, max_y, min_x, min_y, new_bending_direction_points_tilted_KOS_left_stacked,
                                       new_bending_direction_points_tilted_KOS_right_stacked,
                                       new_bending_direction_points_tilted_KOS_stacked, x_values_trim_stacked,
                                       y_values_trim_stacked, z_grid_values_linear)
# Functions in calc_bending_parameters
def initialize_global_lists_of_bending_and_plot_parameter():
    #  Startpoint and one additional Point for the direction (Endpoint). First values are drawn into the colormap.
    global xdata_list, ydata_list
    xdata_list = xdata
    ydata_list = ydata

    # List of Bendingparameters
    global x_direction_list_global_KOS, x_direction_rotated_list_global_KOS, normal_patch_global_KOS, beta_angle_between_planes_list, alpha_angle_list, length_list, edge_line_global, bend_pts_xyz_global_3D
    x_direction_list_global_KOS, x_direction_rotated_list_global_KOS, normal_patch_global_KOS, beta_angle_between_planes_list, alpha_angle_list, edge_line_global, length_list, bend_pts_xyz_global_3D = [], [], [], [], [], [], [], []

    # Blue Points for 3D-global Plot
    global new_bending_direction_points_on_surface_global_KOS_stacked, new_bending_direction_points_on_surface_global_KOS_left_stacked, new_bending_direction_points_on_surface_global_KOS_right_stacked
    new_bending_direction_points_on_surface_global_KOS_stacked, new_bending_direction_points_on_surface_global_KOS_left_stacked, new_bending_direction_points_on_surface_global_KOS_right_stacked = [], [], []

    # 2D Plot tildet local KOS
    global new_bending_direction_points_tilted_KOS_stacked, new_bending_direction_points_tilted_KOS_left_stacked, new_bending_direction_points_tilted_KOS_right_stacked, \
        bend_pts_xz_local_stacked, bend_pts_xz_local_right_stacked, bend_pts_xz_local_left_stacked
    new_bending_direction_points_tilted_KOS_stacked, new_bending_direction_points_tilted_KOS_left_stacked, new_bending_direction_points_tilted_KOS_right_stacked, \
    bend_pts_xz_local_stacked, bend_pts_xz_local_right_stacked, bend_pts_xz_local_left_stacked = [], [], [], [], [], []

    # Line in 2D Colormap
    global x_values_trim_stacked, y_values_trim_stacked
    x_values_trim_stacked, y_values_trim_stacked = [], []

    # Orientation of Tape
    global x_y_direction_tape_in_trendline_KOS_stacked
    x_y_direction_tape_in_trendline_KOS_stacked = []

    global counter_failed_matches_of_edges
    counter_failed_matches_of_edges = 0

def calculate_iteratively_the_tilted_bendingpoints(alpha_angle_list, grid_ressolution_int, grid_x, max_distance, width_for_edge_detection, max_x,
                                                   max_y, min_x, min_y, x_0_grid_point_index, x_values_trim_stacked,
                                                   xdata_list, y_0_grid_point_index, y_values_trim_stacked, ydata_list,
                                                   z_grid_values_linear):
    # calc_tilted_bending_points has no return value except num_bendpoints. All the bending parameters are saved in the global lists
    num_bendpoints = calc_tilted_bending_points(
        grid_ressolution_int, grid_x, max_y, min_y, max_x, min_x, y_0_grid_point_index, x_0_grid_point_index,
        z_grid_values_linear, [xdata_list[-2], xdata_list[-1]], [ydata_list[-2], ydata_list[-1]], max_distance, width_for_edge_detection,
        alpha_start=0, alpha_end=0)

    pyplot.figure()
    plt.imshow(z_grid_values_linear.T, extent=(min_x, max_x, min_y, max_y), origin='lower')
    plt.plot(x_values_trim_stacked[0][:], y_values_trim_stacked[0][:], 'bo', linewidth=1.0, label='Schnitt')
    plt.show(block=False)
    pyplot.pause(0.3) #Give system time to print the plot

    while num_bendpoints > 2:
        num_bendpoints = calc_tilted_bending_points(
            grid_ressolution_int, grid_x, max_y, min_y, max_x, min_x, y_0_grid_point_index, x_0_grid_point_index,
            z_grid_values_linear, [xdata_list[-2], xdata_list[-1]], [ydata_list[-2], ydata_list[-1]], max_distance, width_for_edge_detection,
            alpha_angle_list[-1], alpha_end=0)

        for i in range(len(x_values_trim_stacked)):
            plt.plot(x_values_trim_stacked[i][:], y_values_trim_stacked[i][:], 'bo', linewidth=1.0, label='Schnitt')
        plt.draw()
        pyplot.pause(0.01)
def show_results_2D_Plots_and_Colormap(max_x, max_y, min_x, min_y, new_bending_direction_points_tilted_KOS_left_stacked,
                                       new_bending_direction_points_tilted_KOS_right_stacked,
                                       new_bending_direction_points_tilted_KOS_stacked, x_values_trim_stacked,
                                       y_values_trim_stacked, z_grid_values_linear):

    # In the different tilted KOS the x-Values don´t allign. Correcting it here for nicer 2Dplot
    connect_points_in_tilted_KOS(new_bending_direction_points_tilted_KOS_stacked)
    connect_points_in_tilted_KOS(new_bending_direction_points_tilted_KOS_left_stacked)
    connect_points_in_tilted_KOS(new_bending_direction_points_tilted_KOS_right_stacked)
    connect_points_in_tilted_KOS(bend_pts_xz_local_stacked)
    connect_points_in_tilted_KOS(bend_pts_xz_local_left_stacked)
    connect_points_in_tilted_KOS(bend_pts_xz_local_right_stacked)

    #Colormap
    pyplot.figure()
    plt.subplot(221)
    plt.imshow(z_grid_values_linear.T, extent=(min_x, max_x, min_y, max_y), origin='lower')
    for i in range(len(x_values_trim_stacked)):
        plt.plot(x_values_trim_stacked[i][:], y_values_trim_stacked[i][:], 'bo', linewidth=1.0, label='Schnitt')
    plt.title('Linear')

    #2D-Sideview
    plt.subplot(222)
    for i in range(len(x_values_trim_stacked)):
        plt.plot(new_bending_direction_points_tilted_KOS_stacked[i][:, 0],
                 new_bending_direction_points_tilted_KOS_stacked[i][:, 2], 'bo', linewidth=1.0, label='Schnitt')
        plt.plot(bend_pts_xz_local_stacked[i][:, 0], bend_pts_xz_local_stacked[i][:, 1], color='green', linewidth=3.0,
             label='lineare Angleichung')  # Streckenweise linear (nur Eckpunkte)
    plt.subplot(223)
    plt.title('left')
    for i in range(len(x_values_trim_stacked)):
        plt.plot(new_bending_direction_points_tilted_KOS_left_stacked[i][:, 0],
                 new_bending_direction_points_tilted_KOS_left_stacked[i][:, 2], 'bo', linewidth=1.0, label='Schnitt')
        plt.plot(bend_pts_xz_local_left_stacked[i][:, 0], bend_pts_xz_local_left_stacked[i][:, 1], color='green', linewidth=3.0,
             label='lineare Angleichung')  # Streckenweise linear (nur Eckpunkte)
    plt.subplot(224)
    plt.title('right')
    for i in range(len(x_values_trim_stacked)):
        plt.plot(new_bending_direction_points_tilted_KOS_right_stacked[i][:, 0],
                 new_bending_direction_points_tilted_KOS_right_stacked[i][:, 2], 'bo', linewidth=1.0, label='Schnitt')
        plt.plot(bend_pts_xz_local_right_stacked[i][:, 0], bend_pts_xz_local_right_stacked[i][:, 1], color='green', linewidth=3.0,
             label='lineare Angleichung')  # Streckenweise linear (nur Eckpunkte)
    plt.show()
def connect_points_in_tilted_KOS(new_bending_direction_points_tilted_KOS_stacked):
    for i in range(1, len(new_bending_direction_points_tilted_KOS_stacked)):
        startx_i = new_bending_direction_points_tilted_KOS_stacked[i][0, 0]
        endx_i_1 = new_bending_direction_points_tilted_KOS_stacked[i - 1][-1, 0]
        difx = startx_i - endx_i_1
        new_bending_direction_points_tilted_KOS_stacked[i][:, 0] = np.subtract(
            new_bending_direction_points_tilted_KOS_stacked[i][:, 0], difx)

#######################################################################################################################
# For every iteration in calc_bending_parameters new bendpoints and parameters have to be calculated.
def calc_tilted_bending_points(grid_ressolution_int, grid_x, max_y, min_y, max_x, min_x, y_0_grid_point_index,
                               x_0_grid_point_index, z_grid_values_linear, xdata, ydata, max_distance, width_for_edge_detection, alpha_start=0 ,alpha_end=0):
    ###### Schräge Gerade mit y(x)
    # Wir brauchen 2 Geraden:   • Einmal die mathematische beschreibung mit Coordinaten für Plot.
    #                           • Einmal die mit Indizies für die Extraction aus dem Grid

    # Schrittweite
    dy = (max_y - min_y) / grid_ressolution_int
    dx = (max_x - min_x) / grid_ressolution_int
    x_values = grid_x[:, 0]

    # Steigung(x_slope) und y-Achsenabschnitt(y_intercept) mit Leastsquare, y = x_slope*x + y_intercept
    A = np.vstack([xdata, np.ones(len(xdata))]).T
    x_slope, y_intercept = np.linalg.lstsq(A, ydata,rcond=None)[0]

    # x_y_z_Axis of new KOS
    trendline_new_direction_current_KOS = calc_local_trendline_KOS(x_slope)

    # Calc bendpoints on surface in new trendline direction, including left and right for edge directions
    bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right, bend_pts_xyz_trendline, bend_pts_xz_local = calc_bend_pts_in_new_directions(
        alpha_end, alpha_start, dx, dy, grid_ressolution_int, max_distance, trendline_new_direction_current_KOS,
        width_for_edge_detection, x_0_grid_point_index, x_slope, x_values, xdata, y_0_grid_point_index, ydata,
        z_grid_values_linear)

    # Calc Tapeparameters
    edge_directions = calc_edge_directions(bend_pts_xyz_global_left, bend_pts_xyz_global_right)
    lenght_between_first_two_bends = np.linalg.norm(bend_pts_xz_local[1] - bend_pts_xz_local[0])
    x_direction_list_current_direction, \
    normal_patch_current_direction, \
    rotated_x_direction_around_edge_current_direction, \
    beta_angle_between_planes_list_current_direction,\
    alpha_angle_list_current_direction = calc_bending_parameters_with_bendpoints(bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right, bend_pts_xyz_trendline, edge_directions)

    # Add Bend/Tapeparameters from second bendpoint to global list(first bendpoint = startpoint)
    append_bend_parameters_at_second_bendpoint_to_global_list(alpha_angle_list_current_direction, bend_pts_xyz_global,
                                                              bend_pts_xyz_trendline,
                                                              beta_angle_between_planes_list_current_direction,
                                                              edge_directions, lenght_between_first_two_bends,
                                                              normal_patch_current_direction,
                                                              rotated_x_direction_around_edge_current_direction,
                                                              x_direction_list_current_direction)
    # todo: from the first iteration we can calc already the simple version without alpha
    return len(bend_pts_xz_local)

def calc_bend_pts_in_new_directions(alpha_end, alpha_start, dx, dy, grid_ressolution_int, max_distance,
                                    trendline_new_direction_current_KOS, width_for_edge_detection, x_0_grid_point_index,
                                    x_slope, x_values, xdata, y_0_grid_point_index, ydata, z_grid_values_linear):

    # Start und Endpunkt des Tapeabschnitts
    end_point_xyz_trendline_data, \
    start_point_xyz_trendline_data, \
    x_start_index, \
    x_end_index = calc_Start_End_in_trendline_KOS_from_xdata_ydata(dx, dy, grid_ressolution_int, x_0_grid_point_index,
                                                                   xdata,
                                                                   y_0_grid_point_index, ydata, z_grid_values_linear)
    # Start und Endpunkt der Seitenlinien zur Abschätzung der Biegekanten
    # Eckpunkte
    delta_length_start_bend = calc_delta_length_at_bend(width_for_edge_detection, alpha_start)
    delta_length_end_bend = calc_delta_length_at_bend(width_for_edge_detection, alpha_end)

    end_point_xyz_trendline_data_right, \
    start_point_xyz_trendline_data_right, \
    x_end_index_right, \
    x_start_index_right = calc_start_end_point_side_in_trendline_KOS(False, delta_length_end_bend,
                                                                     delta_length_start_bend,
                                                                     dx, dy, end_point_xyz_trendline_data,
                                                                     grid_ressolution_int,
                                                                     start_point_xyz_trendline_data,
                                                                     width_for_edge_detection,
                                                                     x_0_grid_point_index,
                                                                     trendline_new_direction_current_KOS[0],
                                                                     y_0_grid_point_index,
                                                                     z_grid_values_linear,
                                                                     trendline_new_direction_current_KOS[1])
    end_point_xyz_trendline_data_left, \
    start_point_xyz_trendline_data_left, \
    x_end_index_left, \
    x_start_index_left = calc_start_end_point_side_in_trendline_KOS(True, delta_length_end_bend,
                                                                    delta_length_start_bend,
                                                                    dx, dy, end_point_xyz_trendline_data,
                                                                    grid_ressolution_int,
                                                                    start_point_xyz_trendline_data,
                                                                    width_for_edge_detection,
                                                                    x_0_grid_point_index,
                                                                    trendline_new_direction_current_KOS[0],
                                                                    y_0_grid_point_index,
                                                                    z_grid_values_linear,
                                                                    trendline_new_direction_current_KOS[1])
    bend_pts_xyz_global_left, bend_pts_xyz_trendline_left, bend_pts_xz_local_left, new_bending_direction_points_on_surface_global_KOS_left, new_bending_direction_points_tilted_KOS_left, x_values_trim_left, y_values_trim_left = calc_points_on_surface_and_extract_bendline(
        dy, grid_ressolution_int, max_distance, start_point_xyz_trendline_data_left,
        trendline_new_direction_current_KOS, x_end_index_left, x_slope, x_start_index_left, x_values,
        y_0_grid_point_index, z_grid_values_linear)
    bend_pts_xyz_global_right, bend_pts_xyz_trendline_right, bend_pts_xz_local_right, new_bending_direction_points_on_surface_global_KOS_right, new_bending_direction_points_tilted_KOS_right, x_values_trim_right, y_values_trim_right = calc_points_on_surface_and_extract_bendline(
        dy, grid_ressolution_int, max_distance, start_point_xyz_trendline_data_right,
        trendline_new_direction_current_KOS,
        x_end_index_right, x_slope, x_start_index_right, x_values, y_0_grid_point_index, z_grid_values_linear)
    bend_pts_xyz_global, bend_pts_xyz_trendline, bend_pts_xz_local, new_bending_direction_points_on_surface_global_KOS, new_bending_direction_points_tilted_KOS, x_values_trim, y_values_trim = calc_points_on_surface_and_extract_bendline(
        dy, grid_ressolution_int, max_distance, start_point_xyz_trendline_data, trendline_new_direction_current_KOS,
        x_end_index, x_slope, x_start_index, x_values, y_0_grid_point_index, z_grid_values_linear)
    # Trim plot points to second bendpoint and add to global list
    append_plot_points_till_second_bendpoint_to_global_list(bend_pts_xz_local, bend_pts_xz_local_left,
                                                            bend_pts_xz_local_right,
                                                            new_bending_direction_points_on_surface_global_KOS,
                                                            new_bending_direction_points_on_surface_global_KOS_left,
                                                            new_bending_direction_points_on_surface_global_KOS_right,
                                                            new_bending_direction_points_tilted_KOS,
                                                            new_bending_direction_points_tilted_KOS_left,
                                                            new_bending_direction_points_tilted_KOS_right,
                                                            x_values_trim, y_values_trim)
    return bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right, bend_pts_xyz_trendline, bend_pts_xz_local

def calc_points_on_surface_and_extract_bendline(dy, grid_ressolution_int, max_distance,
                                                start_point_xyz_trendline_data, trendline_new_direction_current_KOS,
                                                x_end_index, x_slope, x_start_index, x_values, y_0_grid_point_index,
                                                z_grid_values_linear):
    new_bending_direction_points_on_surface_global_KOS, y_intercept, x_values_indizes_trim, \
    x_values_trim, y_values_indizes_trim, y_values_trim = calc_new_direction_points_on_surface(
        start_point_xyz_trendline_data, dy, grid_ressolution_int, x_slope, x_values, y_0_grid_point_index,
        z_grid_values_linear,
        x_start_index, x_end_index)
    new_bending_direction_points_tilted_KOS, new_bending_direction_points_on_surface_global_KOS, \
    trendline_new_direction_global_KOS = new_direction_points_in_global_and_tilted_KOS(y_intercept,
                                                                                       new_bending_direction_points_on_surface_global_KOS,
                                                                                       trendline_new_direction_current_KOS)
    bend_pts_xz_local, bend_pts_xyz_trendline, bend_pts_xyz_global = calc_local_and_global_bendpoints(max_distance,
                                                                                                      new_bending_direction_points_tilted_KOS,
                                                                                                      trendline_new_direction_current_KOS,
                                                                                                      y_intercept)
    return bend_pts_xyz_global, bend_pts_xyz_trendline, bend_pts_xz_local, new_bending_direction_points_on_surface_global_KOS, new_bending_direction_points_tilted_KOS, x_values_trim, y_values_trim

def append_bend_parameters_at_second_bendpoint_to_global_list(alpha_angle_list_current_direction, bend_pts_xyz_global,
                                                              bend_pts_xyz_trendline,
                                                              beta_angle_between_planes_list_current_direction,
                                                              edge_directions, lenght_between_first_two_bends,
                                                              normal_patch_current_direction,
                                                              rotated_x_direction_around_edge_current_direction,
                                                              x_direction_list_current_direction):

    global x_direction_list_global_KOS, x_direction_rotated_list_global_KOS, normal_patch_global_KOS, beta_angle_between_planes_list, alpha_angle_list, length_list, edge_line_global, xdata_list, ydata_list, bend_pts_xyz_global_3D

    # Tape parameter for EA
    length_list.append(lenght_between_first_two_bends)

    # Bendpoint and edge
    bend_pts_xyz_global_3D.append(bend_pts_xyz_global[0])

    edge_line_global_current_direction = [
        (bend_pts_xyz_global[:][:2] + np.multiply(edge_directions[:][:2], 50)),
        (bend_pts_xyz_global[:][:2] - np.multiply(edge_directions[:][:2], 50))]
    edge_line_global.append([edge_line_global_current_direction[0][1], edge_line_global_current_direction[1][1]])

    if len(edge_directions) == 2:  # Last part
        bend_pts_xyz_global_3D.append(bend_pts_xyz_global[1])
    if len(edge_directions) > 2:
        # xdata,ydata gibt Startpunkt und Richtung der nächsten Biegerichtung an
        rotated_x_direction_point_in_trendline_KOS = np.asarray(
            [bend_pts_xyz_global[1] + np.multiply(rotated_x_direction_around_edge_current_direction[0], 5),
             bend_pts_xyz_global[1]])
        rotated_x_direction_around_edge_trendline_KOS_current_direction = translate_and_rotate_points_from_OLD_to_trendline_KOS(
            rotated_x_direction_point_in_trendline_KOS, trendline_global_KOS, center_point_of_cloud_weighted)

        # Startpoint
        xdata_list.append(bend_pts_xyz_trendline[1][0])
        ydata_list.append(bend_pts_xyz_trendline[1][1])
        # Direction
        xdata_list.append(rotated_x_direction_around_edge_trendline_KOS_current_direction[0][0])
        ydata_list.append(rotated_x_direction_around_edge_trendline_KOS_current_direction[0][1])

        # Tape parameter for EA
        beta_angle_between_planes_list.append(beta_angle_between_planes_list_current_direction[0])
        alpha_angle_list.append(alpha_angle_list_current_direction[0])

        # Directions
        x_direction_list_global_KOS.append(x_direction_list_current_direction[0])
        x_direction_rotated_list_global_KOS.append(rotated_x_direction_around_edge_current_direction[0])
        normal_patch_global_KOS.append(normal_patch_current_direction[0])

def append_plot_points_till_second_bendpoint_to_global_list(bend_pts_xz_local, bend_pts_xz_local_left,
                                                            bend_pts_xz_local_right,
                                                            new_bending_direction_points_on_surface_global_KOS,
                                                            new_bending_direction_points_on_surface_global_KOS_left,
                                                            new_bending_direction_points_on_surface_global_KOS_right,
                                                            new_bending_direction_points_tilted_KOS,
                                                            new_bending_direction_points_tilted_KOS_left,
                                                            new_bending_direction_points_tilted_KOS_right,
                                                            x_values_trim, y_values_trim):
    # Calc trim index
    end_index_pts_global, end_index_pts_global_left, end_index_pts_global_right = calc_trim_index_at_second_bendpoint(
        bend_pts_xz_local, bend_pts_xz_local_left, bend_pts_xz_local_right, new_bending_direction_points_tilted_KOS,
        new_bending_direction_points_tilted_KOS_left, new_bending_direction_points_tilted_KOS_right)


    # Global 3D-Plot
    global new_bending_direction_points_on_surface_global_KOS_stacked, new_bending_direction_points_on_surface_global_KOS_left_stacked, new_bending_direction_points_on_surface_global_KOS_right_stacked
    new_bending_direction_points_on_surface_global_KOS_stacked.append(
        new_bending_direction_points_on_surface_global_KOS[:][:end_index_pts_global])
    new_bending_direction_points_on_surface_global_KOS_right_stacked.append(
        new_bending_direction_points_on_surface_global_KOS_right[:][:end_index_pts_global_right])
    new_bending_direction_points_on_surface_global_KOS_left_stacked.append(
        new_bending_direction_points_on_surface_global_KOS_left[:][:end_index_pts_global_left])
    # local 2D-Plot
    global new_bending_direction_points_tilted_KOS_stacked, new_bending_direction_points_tilted_KOS_left_stacked, new_bending_direction_points_tilted_KOS_right_stacked, \
        bend_pts_xz_local_stacked, bend_pts_xz_local_right_stacked, bend_pts_xz_local_left_stacked
    new_bending_direction_points_tilted_KOS_stacked.append(
        new_bending_direction_points_tilted_KOS[:end_index_pts_global])
    new_bending_direction_points_tilted_KOS_left_stacked.append(
        new_bending_direction_points_tilted_KOS_left[:end_index_pts_global])
    new_bending_direction_points_tilted_KOS_right_stacked.append(
        new_bending_direction_points_tilted_KOS_right[:end_index_pts_global])
    bend_pts_xz_local_stacked.append(bend_pts_xz_local[:2])
    bend_pts_xz_local_right_stacked.append(bend_pts_xz_local_right[:2])
    bend_pts_xz_local_left_stacked.append(bend_pts_xz_local_left[:2])
    # for trendline KOS, 2D-Plot-Colormap
    global x_values_trim_stacked, y_values_trim_stacked
    x_values_trim_stacked.append(x_values_trim[:end_index_pts_global])
    y_values_trim_stacked.append(y_values_trim[:end_index_pts_global])

def calc_trim_index_at_second_bendpoint(bend_pts_xz_local, bend_pts_xz_local_left, bend_pts_xz_local_right,
                                        new_bending_direction_points_tilted_KOS,
                                        new_bending_direction_points_tilted_KOS_left,
                                        new_bending_direction_points_tilted_KOS_right):
    try:
        end_index_pts_global = \
        [i for i, x in enumerate(new_bending_direction_points_tilted_KOS[:, 0]) if x >= bend_pts_xz_local[1][0]][0]
    except:
        end_index_pts_global = -1
    try:
        end_index_pts_global_left = [i for i, x in enumerate(new_bending_direction_points_tilted_KOS_left[:, 0]) if
                                     x >= bend_pts_xz_local_left[1][0]][0]
    except:
        end_index_pts_global_left = -1
    try:
        end_index_pts_global_right = [i for i, x in enumerate(new_bending_direction_points_tilted_KOS_right[:, 0]) if
                                      x >= bend_pts_xz_local_right[1][0]][0]
    except:
        end_index_pts_global_right = -1
    return end_index_pts_global, end_index_pts_global_left, end_index_pts_global_right

def calc_bend_pts(max_distance, x_z_surface_point, y_smooth):

    bend_pts_xz = []
    bend_pts_xz.append([x_z_surface_point[0][0], y_smooth[0]])  # Comment_DB: start point 2D (x coord, y coord)
    bend_pts_xz.append([x_z_surface_point[-1][0], y_smooth[-1]])  # Comment_DB: end point 2D (x coord, y coord)
    bend_pts_xz = np.asarray(bend_pts_xz)

    # Einfügen von weiteren Knickpunkten durch Finden von großen Abweichungen zur Kurve:
    insert_pts = True
    while insert_pts:
        points_on_line_between_bends_filled_up = []
        points_on_line_between_bends_filled_up.append([bend_pts_xz[0][0], bend_pts_xz[0][1]])  # Comment_DB: only the first bend point (starting point at edge) appended to bend points curve list

        points_on_line_between_bends_filled_up = calc_points_on_line_between_bends_filled_up(bend_pts_xz,
                                                                                             points_on_line_between_bends_filled_up, x_z_surface_point)

        max_divergence = calc_point_of_max_divergence_between_smooth_and_lin_curve(
            points_on_line_between_bends_filled_up, x_z_surface_point, y_smooth)

        # Comment_DB: We know at which x-coord of points_on_line_between_bends_filled_up the max_divergence happens --> counter i



        # no further points, if the chosen maximum distance is not surpassed
        if max_divergence[0] < max_distance:  # Comment_DB: This implies that there will be one extra bend, as the above code will have executed already, max_distance: User Input
            break

        bend_pts_xz = np.insert(bend_pts_xz, -1,
                            np.array([points_on_line_between_bends_filled_up[max_divergence[1]][0],
                                      y_smooth[max_divergence[1]]]),
                            axis=0)  # Comment_DB: insert a corner at x coord (counter i) and y coord (counter i) of max divergence
        bend_pts_xz = bend_pts_xz[bend_pts_xz[:, 0].argsort()]  # Comment_DB: Bend points sorted in an array


    return bend_pts_xz
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
def calc_bending_parameters_with_bendpoints(bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right,bend_pts_xyz_trendline,
                                            edge_directions):
    x_direction_start = norm_vector(bend_pts_xyz_global[1] - bend_pts_xyz_global[0])
    x_direction_list = [x_direction_start]

    rotated_x_direction_around_edge = []
    rotated_y_direction_around_edge = []

    normal_at_start = -calc_tape_normal(bend_pts_xyz_global[0], bend_pts_xyz_global_left[1], bend_pts_xyz_global_right[1])
    normal_patch = [normal_at_start]

    beta_angle_between_planes_list = []
    alpha_angle_between_planes_list = []

    for i in range(1, (len(edge_directions) - 1)):
        # Calc x_y_and_normal direction of the Tape at each bending point
        x_direction_before_bend = norm_vector(bend_pts_xyz_global[i] - bend_pts_xyz_global[i - 1])

        x_direction_i_without_rotation = norm_vector(bend_pts_xyz_global[i + 1] - bend_pts_xyz_global[i])

        # We get the normal of the plane defined by the right and left bendpoint of the bend and the middlepoint on the next/previous bend.
        normal_at_bendpoint_0_tape = -calc_tape_normal(bend_pts_xyz_global[i-1], bend_pts_xyz_global_left[i],
                                                      bend_pts_xyz_global_right[i])
        normal_at_bendpoint_1_tape = calc_tape_normal(bend_pts_xyz_global[i+1], bend_pts_xyz_global_left[i],
                                                      bend_pts_xyz_global_right[i])


        y_direction_tape = np.cross(normal_at_bendpoint_0_tape, x_direction_before_bend)

        trendline_patch = np.stack([x_direction_before_bend, y_direction_tape, normal_at_bendpoint_0_tape])

        #Calc angle between two consecutiv tape normals: beta
        beta_angle_between_planes = math.acos(np.dot(normal_at_bendpoint_0_tape, normal_at_bendpoint_1_tape))

        # Compare slope of old and new x_direction
        x_direction_before_bend_local = norm_vector(bend_pts_xyz_trendline[i] - bend_pts_xyz_trendline[i - 1])

        x_direction_i_without_rotation_local = norm_vector(bend_pts_xyz_trendline[i + 1] - bend_pts_xyz_trendline[i])
        if x_direction_before_bend_local[2] < x_direction_i_without_rotation_local[2]:
            beta_angle_between_planes = -beta_angle_between_planes

        # Edgedirection always in positiv y-Direction



        # Calc new x_direction after bending around edge
        rotated_x_direction_around_edge_i = Quaternion(axis=edge_directions[i], angle=beta_angle_between_planes).rotate(x_direction_before_bend)
        rotated_x_direction_around_edge.append(rotated_x_direction_around_edge_i)

        edge_directions_trendlineKOS = translate_and_rotate_points_from_OLD_to_trendline_KOS(np.asarray(edge_directions),
                                                                                             trendline_patch,
                                                                                             np.asarray([0, 0, 0]))
        if edge_directions_trendlineKOS[i][1] < 0:
            edge_directions[i] = - edge_directions[i]

        rotated_y_direction_around_edge.append(np.cross(normal_at_bendpoint_1_tape,rotated_x_direction_around_edge[-1]))


        # Calc alpha. Angle between y_tape and direction of edge. Transformation in local tape KOS to get direction of angle

        side_directions = np.stack([edge_directions[i], y_direction_tape])
        side_directions_rotated = translate_and_rotate_points_from_OLD_to_trendline_KOS(side_directions, trendline_patch,
                                                                                        np.asarray([0, 0, 0]))

        if side_directions_rotated[0][1] < 0:
            side_directions_rotated[0][:] = -side_directions_rotated[0][:]
        if side_directions_rotated[1][1] < 0:
            side_directions_rotated[1][:] = -side_directions_rotated[1][:]

        # Edge_direction auf y_direction tape projezieren, sollte eigentlich direkt drauf sein. Punkte rechts links können minimal abweichen.
        edge_direction_on_tapeplane = project_pointtoplane(edge_directions[i],normal_at_bendpoint_0_tape,np.asarray([0,0,0]))
        edge_direction_on_tapeplane = norm_vector(edge_direction_on_tapeplane)
        y_direction_tape = norm_vector(y_direction_tape)
        alpha_angle = math.acos(np.dot(edge_direction_on_tapeplane, y_direction_tape))

        if side_directions_rotated[1][0] > side_directions_rotated[0][0]:
            alpha_angle = (math.pi-alpha_angle)

        # safe x_directions, alpha, beta, and normals in a List.
        x_direction_list.append(x_direction_before_bend)    #x_direction_start is two times in the list. If there are 1+ bending points. For special case there is no bending point, the x_direction is put in the list also outside the loop
        alpha_angle_between_planes_list.append(alpha_angle)
        beta_angle_between_planes_list.append(beta_angle_between_planes)
        normal_patch.append(normal_at_bendpoint_0_tape)

    # Rotate new Tape-Orientation to Trendline KOS as referenz for the next bending direction
    if len(edge_directions)>2:
        x_y_stacked = np.stack([rotated_x_direction_around_edge[0],rotated_y_direction_around_edge[0]])
        global x_y_direction_tape_in_trendline_KOS_stacked
        x_y_direction_tape_in_trendline_KOS_stacked.append(
            translate_and_rotate_points_from_OLD_to_trendline_KOS(x_y_stacked, trendline_global_KOS, np.asarray([0, 0, 0])))


    rotated_x_direction_around_edge = np.asarray(rotated_x_direction_around_edge)
    x_direction_list = np.asarray(x_direction_list)


    return x_direction_list, normal_patch, rotated_x_direction_around_edge,beta_angle_between_planes_list,alpha_angle_between_planes_list
def calc_edge_directions(bend_pts_xyz_global_left, bend_pts_xyz_global_right):

    # Just calculate for the first 3 bends.
    global counter_failed_matches_of_edges



    while len(bend_pts_xyz_global_right) != len(bend_pts_xyz_global_left):  # Comment_DKu_Wenzel: Just looking for the first Bend would optimize stability and would be faster
        if len(bend_pts_xyz_global_right) < len(bend_pts_xyz_global_left): bend_pts_xyz_global_left = np.delete(bend_pts_xyz_global_left,-1,axis=0)
        else:  bend_pts_xyz_global_right = np.delete(bend_pts_xyz_global_right,-1,axis=0) # right < left
        print("left and right not same amount of bendingpoints")
        counter_failed_matches_of_edges += 1
        if counter_failed_matches_of_edges > 20:
            print("Please restart, edges have not been found. Maybe try other width.")
            exit()

    edge_directions = []
    for i in range(len(bend_pts_xyz_global_left)):
        edge_direction = bend_pts_xyz_global_right[i] - bend_pts_xyz_global_left[i]
        edge_direction = norm_vector(edge_direction)
        edge_directions.append(edge_direction)

    return edge_directions
def calc_tape_normal(bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right):
    v1 = bend_pts_xyz_global_left - bend_pts_xyz_global
    v2 = bend_pts_xyz_global_right - bend_pts_xyz_global
    normal_at_bendpoint_index_tape = np.cross(v2, v1)
    normal_at_bendpoint_index_tape = norm_vector(normal_at_bendpoint_index_tape)


    return normal_at_bendpoint_index_tape

# Functions in calc tilted bending points
def calc_delta_length_at_bend(width,alpha):
    if alpha > math.pi / 2:
        delta_length_bend = (width / 2) * math.tan(math.pi - alpha)

    else:  # alpha_list[0] < math.pi / 2:
        delta_length_bend = - (width / 2) * math.tan(alpha)
    return delta_length_bend
def calc_local_trendline_KOS(x_slope):
    x_trendline_new_direction = np.asarray((1, x_slope, 0), dtype=np.float32)
    x_trendline_new_direction = norm_vector(x_trendline_new_direction)
    y_trendline_new_direction = np.asarray((-x_slope, 1, 0), dtype=np.float32)
    y_trendline_new_direction = norm_vector(y_trendline_new_direction)
    z_trendline_new_direction = np.asarray((0, 0, 1), dtype=np.float32)
    trendline_new_direction_current_KOS = np.vstack(
        (x_trendline_new_direction, y_trendline_new_direction, z_trendline_new_direction))
    return trendline_new_direction_current_KOS
def calc_start_end_point_side_in_trendline_KOS(calc_left_side, delta_length_end_bend, delta_length_start_bend, dx, dy, end_point_drawn,
                                               grid_ressolution_int, start_point_drawn, width, x_0_grid_point_index,
                                               x_trendline_new_direction, y_0_grid_point_index, z_grid_values_linear,
                                               y_trendline_new_direction):

    # We need to consider how the Tape is placed 3D

    if calc_left_side:
        if delta_length_start_bend != 0:
            Start_point_side = start_point_drawn - x_y_direction_tape_in_trendline_KOS_stacked[-1][1] * width / 2 + delta_length_start_bend * x_y_direction_tape_in_trendline_KOS_stacked[-1][0]
            End_point_side = end_point_drawn - x_y_direction_tape_in_trendline_KOS_stacked[-1][1] * width / 2 + delta_length_end_bend * x_y_direction_tape_in_trendline_KOS_stacked[-1][0]
        else:
            Start_point_side = start_point_drawn - y_trendline_new_direction * width / 2 + delta_length_start_bend * x_trendline_new_direction
            End_point_side = end_point_drawn - y_trendline_new_direction * width / 2 + delta_length_end_bend * x_trendline_new_direction

    else:
        if delta_length_start_bend != 0:
            Start_point_side = start_point_drawn + x_y_direction_tape_in_trendline_KOS_stacked[-1][1] * width / 2 - delta_length_start_bend * \
                               x_y_direction_tape_in_trendline_KOS_stacked[-1][0]
            End_point_side = end_point_drawn + x_y_direction_tape_in_trendline_KOS_stacked[-1][1] * width / 2 - delta_length_end_bend * \
                             x_y_direction_tape_in_trendline_KOS_stacked[-1][0]
        else:
            Start_point_side = start_point_drawn + y_trendline_new_direction * width / 2 - delta_length_start_bend * x_trendline_new_direction
            End_point_side = end_point_drawn + y_trendline_new_direction * width / 2 - delta_length_end_bend * x_trendline_new_direction



    x_data_side_start_end = (Start_point_side[0], End_point_side[0])
    y_data_side_start_end = (Start_point_side[1], End_point_side[1])

    end_point_drawn_side, start_point_drawn_side, x_start_index_side, x_end_index_side = calc_Start_End_in_trendline_KOS_from_xdata_ydata(
        dx, dy, grid_ressolution_int,
        x_0_grid_point_index, x_data_side_start_end,
        y_0_grid_point_index, y_data_side_start_end,
        z_grid_values_linear)
    return end_point_drawn_side, start_point_drawn_side, x_end_index_side, x_start_index_side
def calc_local_and_global_bendpoints(max_distance, new_bending_direction_points_tilted_KOS,
                                     trendline_new_direction_current_KOS, y_intercept):

    x_z_points = np.stack([new_bending_direction_points_tilted_KOS[:, 0],new_bending_direction_points_tilted_KOS[:, 2]],axis=1)
    #y_smooth_left = smooth_savgol(x_z_points_filled_up, poly_order, savgol_window_quotient)
    bend_pts_xz_local = calc_bend_pts(max_distance, x_z_points,
                                           x_z_points[:,1])
    # in tilted KOS y=0. Insert this to bendpoints to get (x,y,z)
    bend_pts_xyz = np.insert(bend_pts_xz_local, 1, np.zeros(len(bend_pts_xz_local)), axis=1)
    bend_pts_xyz_trendline, bend_pts_xyz_global = new_bending_points_tilted_to_global(y_intercept, bend_pts_xyz,
                                                                  trendline_new_direction_current_KOS)
    return bend_pts_xz_local, bend_pts_xyz_trendline, bend_pts_xyz_global
def calc_Start_End_in_trendline_KOS_from_xdata_ydata(dx, dy, grid_ressolution_int, x_0_grid_point_index, xdata,
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


    start_point_xyz_data = (np.vstack((xdata[0], ydata[0], z_start_data)).T)[0][:]
    end_point_xyz_data = (np.vstack((xdata[1], ydata[1], z_end_data)).T)[0][:]
    return end_point_xyz_data, start_point_xyz_data, x_start_index, x_end_index
def trim_x_y_values_to_geometry(grid_ressolution_int, x_slope, x_values, x_values_indizes, y_values, y_values_indizes):
    # -1 default wert
    y_start_index = -1
    y_end_index = -1
    for k in range(len(x_values)-1):

        if x_slope >= 0:
            if (y_values_indizes[k] >= 0) & (y_start_index < 0): y_start_index = k
            if (y_values_indizes[k] >= grid_ressolution_int) & (y_end_index < 0): y_end_index = k - 1

        if x_slope < 0:
            if (y_values_indizes[k] < grid_ressolution_int) & (y_start_index < 0): y_start_index = k
            if (y_values_indizes[k] <= 0) & (y_end_index < 0): y_end_index = k - 1
    if y_end_index <= 0: y_end_index = grid_ressolution_int - 2
    # Trim indizes and coordinates to grid size
    x_values_indizes_trim = x_values_indizes[y_start_index:y_end_index]
    y_values_indizes_trim = y_values_indizes[y_start_index:y_end_index]
    x_values_trim = x_values[y_start_index:y_end_index]
    y_values_trim = y_values[y_start_index:y_end_index]
    return x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim
def trim_x_y_values_to_start_end_point(x_start_index, x_end_index, x_values, x_values_indizes, y_values,
                                       y_values_indizes,just_start_point_cut=True):
   # Trim indizes and coordinates to grid size
    if just_start_point_cut:
        if x_start_index < x_end_index:
            x_values_indizes_trim = x_values_indizes[x_start_index:]
            y_values_indizes_trim = y_values_indizes[x_start_index:]
            x_values_trim = x_values[x_start_index:]
            y_values_trim = y_values[x_start_index:]
        else:
            x_values_indizes_trim = x_values_indizes[x_end_index:]
            y_values_indizes_trim = y_values_indizes[x_end_index:]
            x_values_trim = x_values[x_end_index:]
            y_values_trim = y_values[x_end_index:]
    else:
        if x_start_index < x_end_index:
            x_values_indizes_trim = x_values_indizes[x_start_index:x_end_index]
            y_values_indizes_trim = y_values_indizes[x_start_index:x_end_index]
            x_values_trim = x_values[x_start_index:x_end_index]
            y_values_trim = y_values[x_start_index:x_end_index]
        else:
            x_values_indizes_trim = x_values_indizes[x_end_index:x_start_index]
            y_values_indizes_trim = y_values_indizes[x_end_index:x_start_index]
            x_values_trim = x_values[x_end_index:x_start_index]
            y_values_trim = y_values[x_end_index:x_start_index]

    return x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim
def calc_new_direction_points_on_surface(Start_point, dy, grid_ressolution_int, x_slope, x_values, y_0_grid_point_index,
                                         z_grid_values_linear, x_start_index, x_end_index):


    # Neuer y-Achsen Abschnitt für Eckpunkte: y = x_slope*x + y_intercept
    y_intercept = Start_point[1] - Start_point[0] * x_slope
    # x-y-Achsenabschnitt in Coordinaten
    y_values = np.add(np.multiply(x_values, x_slope), y_intercept)
    # x-y-Achsenabschnitt in Indizies
    x_values_indizes = np.asarray(list(range(grid_ressolution_int)), dtype=np.int32)
    y_values_indizes = np.add(np.divide(y_values, dy), (grid_ressolution_int - y_0_grid_point_index))
    y_values_indizes = np.asarray(np.round(y_values_indizes), dtype=np.int32)


    # Trim x and y valules to start/end point and if nec
    x_indizes_trim, x_values_trim, y_indizes_trim, y_values_trim = trim_x_y_values_to_start_end_point(
            x_start_index, x_end_index, x_values, x_values_indizes, y_values, y_values_indizes)

    if min(y_values_indizes) < 0 or  max(y_values_indizes) > grid_ressolution_int-1:
      x_indizes_trim, x_values_trim, y_indizes_trim, y_values_trim = trim_x_y_values_to_geometry(
                grid_ressolution_int, x_slope, x_values_trim, x_indizes_trim, y_values_trim, y_indizes_trim)

    # z Values from grid
    z_values_new_bending_direction = []
    for i in range(len(x_indizes_trim)):
        z_values_new_bending_direction.append(
            z_grid_values_linear[x_indizes_trim[i], y_indizes_trim[i]])
    z_values_new_bending_direction = np.asarray(z_values_new_bending_direction, dtype=np.float32)
    # x, y and z stacked together to 3D-Points

    # In z_grid_values_linear can be NaN. Trim those of.
    i = 0
    while z_values_new_bending_direction[i] != z_values_new_bending_direction[i]: i += 1
    j = -1
    while z_values_new_bending_direction[j] != z_values_new_bending_direction[j] :  j -= 1

    z_values_new_bending_direction = z_values_new_bending_direction[i:j]
    x_values_trim = x_values_trim[i:j]
    y_values_trim = y_values_trim[i:j]
    x_indizes_trim = x_indizes_trim[i:j]
    y_indizes_trim = y_indizes_trim[i:j]

    new_bending_direction_points_current_KOS = np.vstack(
        (x_values_trim, y_values_trim, z_values_new_bending_direction)).T
    return new_bending_direction_points_current_KOS, y_intercept, x_indizes_trim, x_values_trim, y_indizes_trim, y_values_trim
def new_direction_points_in_global_and_tilted_KOS(y_intercept, new_bending_direction_points_current_KOS,
                                                  trendline_new_direction_current_KOS):
    # new zero
    new_zero = np.asarray((0, y_intercept, 0), dtype=np.float32)

    trendline_new_direction_global_KOS = translate_and_rotate_points_from_OLD_to_trendline_KOS(
        trendline_new_direction_current_KOS, trendline_global_KOS, np.asarray([0, 0, 0]), True)

    new_bending_direction_points_global_KOS = translate_and_rotate_points_from_OLD_to_trendline_KOS(
        new_bending_direction_points_current_KOS, trendline_global_KOS, center_point_of_cloud_weighted, True)

    new_bending_direction_points_tilted_KOS = translate_and_rotate_points_from_OLD_to_trendline_KOS(
        new_bending_direction_points_current_KOS, trendline_new_direction_current_KOS, new_zero)

    return new_bending_direction_points_tilted_KOS, new_bending_direction_points_global_KOS, trendline_new_direction_global_KOS
def new_bending_points_tilted_to_global(y_intercept,
                                                new_bending_direction_points_current_KOS,
                                                trendline_new_direction_current_KOS):
    # new zero
    new_zero = np.asarray((0, y_intercept, 0), dtype=np.float32)

    new_bending_direction_points_trendline_KOS = translate_and_rotate_points_from_OLD_to_trendline_KOS(
        new_bending_direction_points_current_KOS, trendline_new_direction_current_KOS, new_zero, True)

    new_bending_direction_points_global_KOS = translate_and_rotate_points_from_OLD_to_trendline_KOS(
        new_bending_direction_points_trendline_KOS, trendline_global_KOS, center_point_of_cloud_weighted, True)

    return new_bending_direction_points_trendline_KOS, new_bending_direction_points_global_KOS

#######################################################################################################################
def calc_L_aim(new_bending_direction_points_tilted_KOS_stacked):
    L_aim = 0
    for i in range(len(new_bending_direction_points_tilted_KOS_stacked)):
        for j in range(1, len(new_bending_direction_points_tilted_KOS_stacked[i])):
            distance = calc_distance_between_two_points(new_bending_direction_points_tilted_KOS_stacked[i][j - 1],
                                                        new_bending_direction_points_tilted_KOS_stacked[i][j])

            L_aim = L_aim + distance
    return L_aim
#######################################################################################################################
def show_startstrip(bestPatch_patternpoints,patch_start,patch_end):
    ###############3D-PLOTTING################
    figure = pyplot.figure() #Comment_DB: 3D plot of objective shape
    axes = mplot3d.Axes3D(figure)

    patch_visual = mplot3d.art3d.Poly3DCollection(triangle_vectors_of_stl, linewidths=1, alpha=0.5, edgecolor=[1, 1, 1], label ='Geometriebereich') #Comment_DB: added edgecolor to make the edges visible

    axes.scatter([999999990],[9999999900],[9999999900],linewidths=0.0001, alpha = 0.5, label = "Geometriebereich") #Comment_DB: label in legend
    axes.scatter(center_point_of_cloud_weighted[0],center_point_of_cloud_weighted[1],center_point_of_cloud_weighted[2],c='g')

    # TRENDLINE
    axes.scatter(patch_start[0], patch_start[1], patch_start[2], c="black")
    axes.scatter(patch_end[0],patch_end[1],patch_end[2],c='black')

    # von PCC gemittelter Normalenvektor
    x3, y3, z3 = [center_point_of_cloud_weighted[0], center_point_of_cloud_weighted[0] + 500 * avg_tri_normal_weighted[0]], \
                 [center_point_of_cloud_weighted[1],center_point_of_cloud_weighted[1] + 500 * avg_tri_normal_weighted[1]], \
                 [center_point_of_cloud_weighted[2], center_point_of_cloud_weighted[2] + 500 * avg_tri_normal_weighted[2]]

    plt.plot(x3,y3,z3,marker='o',c='green')
    for i in range(len(new_bending_direction_points_on_surface_global_KOS_stacked)):
        plt.plot((new_bending_direction_points_on_surface_global_KOS_stacked[i][:, 0]), (new_bending_direction_points_on_surface_global_KOS_stacked[i][:, 1]), (new_bending_direction_points_on_surface_global_KOS_stacked[i][:, 2]), marker='o', c='blue')
        plt.plot((new_bending_direction_points_on_surface_global_KOS_right_stacked[i][:, 0]), (new_bending_direction_points_on_surface_global_KOS_right_stacked[i][:, 1]), (new_bending_direction_points_on_surface_global_KOS_right_stacked[i][:, 2]), marker='o', c='blue')
        plt.plot((new_bending_direction_points_on_surface_global_KOS_left_stacked[i][:, 0]), (new_bending_direction_points_on_surface_global_KOS_left_stacked[i][:, 1]), (new_bending_direction_points_on_surface_global_KOS_left_stacked[i][:, 2]), marker='o', c='blue')

    for k in range(len(x_direction_list_global_KOS)):
            axes.plot([bend_pts_xyz_global_3D[k+1][0], bend_pts_xyz_global_3D[k+1][0] + 50 * x_direction_list_global_KOS[k][0]],
                      [bend_pts_xyz_global_3D[k+1][1], bend_pts_xyz_global_3D[k+1][1] + 50 * x_direction_list_global_KOS[k][1]],
                      [bend_pts_xyz_global_3D[k+1][2], bend_pts_xyz_global_3D[k+1][2] + 50 * x_direction_list_global_KOS[k][2]], c='red')
    for k in range(len(x_direction_rotated_list_global_KOS)):
            axes.plot([bend_pts_xyz_global_3D[k+1][0], bend_pts_xyz_global_3D[k+1][0] + 50 * x_direction_rotated_list_global_KOS[k][0]],
                      [bend_pts_xyz_global_3D[k+1][1], bend_pts_xyz_global_3D[k+1][1] + 50 * x_direction_rotated_list_global_KOS[k][1]],
                      [bend_pts_xyz_global_3D[k+1][2], bend_pts_xyz_global_3D[k+1][2] + 50 * x_direction_rotated_list_global_KOS[k][2]],
                       c='red')

    for i in range(len(edge_line_global)):
        axes.plot([edge_line_global[i][0][0], edge_line_global[i][1][0]],
                [edge_line_global[i][0][1], edge_line_global[i][1][1]],
                [edge_line_global[i][0][2], edge_line_global[i][1][2]],
                 c='green')  # Comment_DB: *pc_axes is *args, and .T is np.transpose

    #"""
    for i in range(len(bestPatch_patternpoints) - 2):
        verts = [list(
            zip([bestPatch_patternpoints[i][0], bestPatch_patternpoints[i + 1][0], bestPatch_patternpoints[i + 2][0]], \
                [bestPatch_patternpoints[i][1], bestPatch_patternpoints[i + 1][1], bestPatch_patternpoints[i + 2][1]], \
                [bestPatch_patternpoints[i][2], bestPatch_patternpoints[i + 1][2],
                 bestPatch_patternpoints[i + 2][2]]))]  # Comment_DB: DARK BLUE LoP PATCH
        axes.add_collection3d(Poly3DCollection(verts), zs='z')  # Comment_DB: INSERT LoP PATCH IN GRAPH
        # patch_meshpoints.append(verts) #Comment_DB: is not used
    axes.scatter(bestPatch_patternpoints[:, 0], bestPatch_patternpoints[:, 1], bestPatch_patternpoints[:, 2], c='r')
    #"""

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
