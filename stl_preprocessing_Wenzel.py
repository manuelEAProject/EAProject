import numpy as np
import math
import scipy
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from pyquaternion import Quaternion


equidistant_step_size = 3
percentile_pc = 5
max_points_in_pc = 50000



####MAIN
def startparam(input_file, max_distance, width_for_edge_detection, grid_resolution, calc_2D_Solution, calc_2D_with_edge_detection,
               calc_3D_Solution):
    print(
        "Preprocessing in progress... With high-resolution stl files and large area size differences it can lead to longer calculation times")

    # Read in the file:
    calc_trendline_KOS_of_geometry_from_stl_file(input_file) # saved in global var!

    # Interpolation:

    z_grid_values_linear_trendline_KOS = interpolate_geometrie_in_trendline_KOS_and_draw_start_end_point( # start end saved in global var!
         grid_resolution)

    # Calculation of the bending parameters:
    calc_bending_parameters(grid_resolution, max_distance, width_for_edge_detection, z_grid_values_linear_trendline_KOS,
                            calc_2D_with_edge_detection, calc_3D_Solution)  # saved in global var!

    ##### Extract start parameters ####
    L_aim, normal_at_start_global_KOS, direction_at_start_global_KOS, alpha_list, amount_of_bends, beta_list, endpoint_on_surface, l_list, startpoint_on_surface = extract_start_parameters(calc_2D_Solution,calc_3D_Solution)

    start_parameter = [l_list, L_aim, beta_list, alpha_list, startpoint_on_surface,
                       endpoint_on_surface, direction_at_start_global_KOS, normal_at_start_global_KOS,amount_of_bends,z_grid_values_linear_trendline_KOS, max_x, max_y, min_x, min_y,gamma_start]
    return start_parameter
def extract_start_parameters(calc_2D_Solution,calc_3D_Solution):
    startpoint_on_surface = bend_pts_xyz_global_2D[0][:]
    endpoint_on_surface = bend_pts_xyz_global_2D[-1][:]
    if x_direction_list_global_KOS == []:
        print("The calculated Tape has no bendingpoints")
        exit()

    direction_at_start_global_KOS = norm_vector(x_direction_list_global_KOS[0])
    if not calc_2D_Solution and calc_3D_Solution:
        normal_at_start_global_KOS = norm_vector(normal_direction_start_3D)
    else:
        normal_at_start_global_KOS = norm_vector(normal_direction_start)  # normal_patch_global_KOS[0]
    # l-list for Chromo
    l_list = [np.asarray(length_list_3D), lenght_list_2D, length_list_2DE]
    # Beta_liste from calc_bending_parameters
    beta_list = [np.asarray(beta_angle_between_planes_list_3D), beta_angle_list_2D,
                 beta_angle_between_planes_list_2DE]
    # Alpha_liste from calc_bending_parameters
    alpha_list = [np.asarray(alpha_angle_list_3D), alpha_angle_list_2D, alpha_angle_list_2DE]
    # L_aim:  Curve of interpolated surface. 2D_solution
    L_aim = calc_L_aim(surfacepoints_between_Start_and_End)  # summed up distances
    if len(l_list[0]) > len(l_list[1]):
        amount_of_bends = len(l_list[0]) - 1
    else:
        amount_of_bends = len(l_list[1]) - 1
    return L_aim, normal_at_start_global_KOS, direction_at_start_global_KOS, alpha_list, amount_of_bends, beta_list, endpoint_on_surface, l_list, startpoint_on_surface

#######################################################################################################################
# Load stl file and analyse the geometry. Calc centerpoint and trendline of the geometry
def calc_trendline_KOS_of_geometry_from_stl_file(input_file):
    patch_vectors_of_stl_input = mesh.Mesh.from_file(input_file)  # Comment_DB: stl mesh

    # Triangel vectors needed for the visualization of the patch.
    global triangle_vectors_of_stl,num_of_triangles
    triangle_vectors_of_stl = patch_vectors_of_stl_input.vectors  # Comment_DB: triangle edges (wireframe)
    num_of_triangles = len(triangle_vectors_of_stl)

    # Calc areas of triangel for compare/weight
    tri_areas = calc_tri_areas(triangle_vectors_of_stl)

    global center_of_pointcloud_weighted_in_global_KOS, trendline_KOS_in_global_KOS, avg_tri_normal_weighted
    avg_tri_normal_weighted = calc_avg_tri_norm_weighted_by_area(tri_areas, patch_vectors_of_stl_input)

    # Calculate approximate center of geometry
    tri_centerpoints = calc_tri_centerpoints(triangle_vectors_of_stl)  # Comment_DKu_Wenzel: basically the unweighted point cloud
    point_cloud_tri_centerpoints_weighted = calc_patch_pointcloud_weighted_by_area(tri_areas,tri_centerpoints)
    center_of_pointcloud_weighted_in_global_KOS = point_cloud_tri_centerpoints_weighted.mean(axis=0) # Mean of x,y,z-Values

    # SVD for approximate orientation of geometry (trendline KOS)
    trendline_KOS_in_global_KOS = calc_trendline_global(center_of_pointcloud_weighted_in_global_KOS, point_cloud_tri_centerpoints_weighted)
# Functions in calc trendline
def calc_trendline_global(center_point_of_cloud_weighted, point_cloud_tri_centerpoints_weighted):
    trendline_x_axis, trendline_y_axis, trendline_z_axis = calc_trendline_axis_with_svd(
        point_cloud_tri_centerpoints_weighted, center_point_of_cloud_weighted)
    # If trendline x axis was defined in negative x direction by svd
    if trendline_x_axis[0] < 0:  # Rotation of 180° around y-Axis
        trendline_x_axis = -(trendline_x_axis)
        trendline_z_axis = -(trendline_z_axis)
    trendline_global_KOS = np.vstack((trendline_x_axis, trendline_y_axis, trendline_z_axis))
    return trendline_global_KOS
def calc_tri_normals_from_stl(stl_normals,triangle_vectors_of_stl):
    normals=[]
    #We generate our own normals with the id_list. Notice that the triangle order of stl_normals and the vertices
    #(triangles) is not the same


    for i in range(num_of_triangles):

        v1= triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][1]
        v2= triangle_vectors_of_stl[i][0] - triangle_vectors_of_stl[i][2]
        n=np.cross(v1,v2)
        n=norm_vector(n)
        normals.append(n)

    normals=np.asarray(normals)


    # Normal vectors are always aligned in positive z direction

    # the following average stl_normal always point at the outside of the object:
    avg_stl_normal = sum(stl_normals) / num_of_triangles
    # average of the created normals:
    avg_sorted_normal = sum(normals) / num_of_triangles


    true_when_stl_and_tri_normal_not_same_direction = avg_sorted_normal[0] * avg_stl_normal[0] < 0
    true_when_z_from_tri_normal_neg = avg_sorted_normal[2] < 0


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

    # Delete double corners
    tri_corner_points = np.unique(tri_corner_points, axis=0)

    return tri_corner_points
def calc_tri_areas(triangle_vectors_of_stl):
    #Calculation of triangle areas and saving in a list
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
    centerpoints_weights_area_tri = calc_weights_for_center_points_by_percentil_area(tri_areas)

    # for each triangle the midpoints are evaluated (put into the point cloud) as often as it is written in centerpoints_weights_area_tri
    pointcloud_weighted=[]

    for i in range(num_of_triangles):
        for j in range(centerpoints_weights_area_tri[i]):
            pointcloud_weighted.append(tri_centerpoints[i])

    pointcloud_weighted=np.asarray(pointcloud_weighted)

    return pointcloud_weighted
def calc_weights_for_center_points_by_percentil_area(tri_areas):
    # This function generates a point cloud from the triangle centers of the stl file. Since with this
    # point cloud later the main value decomposition takes place, large triangles must be stronger than small triangles
    # be weighted # This is done by calculating the ratio of the triangular areas to each other and the
    # centers of large triangles are added to the point cloud more often than the centers of small triangles


    ###Weighting large triangles
    # 1) Calculation of the triangle areas and saving in a list
    area_whole_patch = sum(tri_areas)
    # 2) The triangle surfaces are compared with each other, in order to later weight large triangles more than small ones
    # The 10th quantile of the areas is taken as the "smallest" reference area (90% of all other triangles
    # are bigger). In centerpoints_weights_area_tri, for each triangle the factor by which it is larger than
    # its reference triangle (the factor is rounded up).The center of each triangle is
    # added to the centerpoints_weights_area_tri (at least once) "factor"-times
    # To avoid too much computing time, the number of points in the point cloud is estimated in advance
    # and if a limit value( max_points_in_pc) is exceeded, the program is aborted.
    lower_percentil_area = np.percentile(tri_areas, percentile_pc)
    estimated_number_points_in_pc = math.ceil(area_whole_patch / lower_percentil_area)

    # Termination condition: maximum "max_points_in_pc" should be calculated in the point cloud.
    if max_points_in_pc < estimated_number_points_in_pc:
        print("ERROR: Please use a .stl-object with reduced resolution ")
        print("Number of triangles: ", num_of_triangles)
        print("Estimated number of points in pointcloud:", estimated_number_points_in_pc)
        print("Allowed number of points in pointcloud:", max_points_in_pc)
        exit(1)

    # In the following, each triangle is compared with the smallest triangle and recorded in centerpoints_weights_area_tri, as
    # often the smallest triangle fits into the respective triangle
    centerpoints_weights_area_tri = []
    for i in range(num_of_triangles):
        centerpoints_weights_area_tri.append(math.ceil(tri_areas[i] / lower_percentil_area))

    return centerpoints_weights_area_tri
def calc_trendline_axis_with_svd(patch_pc_weighted, center_point_of_cloud_weighted):
    # Do Principal Component Analysis(PCA) on the mean-centered data. AKA SVD
    # The first principal component contains [uu, dd, vv] , where vv[0] is the direction
    first_principal_components_pc_weighted = scipy.linalg.svd(patch_pc_weighted - center_point_of_cloud_weighted) # scipy lib is faster then numpy
    # Definition of the trendline axes
    trendline_x_axis = first_principal_components_pc_weighted[2][0] # first_principal_components_pc_weighted[2][0]: is direction of trendline
    trendline_x_axis = norm_vector(trendline_x_axis)
    # avg_tri_norm is not perpendicular to the x-axis
    # project from pcc + avg_tri_norm and back to x-axis
    trendline_avg_norm_point = center_point_of_cloud_weighted + np.dot(avg_tri_normal_weighted,
                                                                       trendline_x_axis) / np.dot(trendline_x_axis,
                                                                                                  trendline_x_axis) * trendline_x_axis
    # y-axis is connection of pcc+avg_tri_norm with the projected point
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
def interpolate_geometrie_in_trendline_KOS_and_draw_start_end_point(grid_resolution):
    z_grid_values_linear_trendline_KOS = interpolate_geometrie(grid_resolution)

    show_interpolation_and_draw_start_end_points(z_grid_values_linear_trendline_KOS) # In this function x_data_drawn and y_data_drawn are calculated

    return z_grid_values_linear_trendline_KOS
# Functions in Interpolate start_geometrie
def interpolate_geometrie(grid_resolution):
    # Comment_DKu_Wenzel: Interpolation mit Centerpoints teils ungenauer

    # Taking corner points of the triangles as interpolation points.
    # First thing we do is rotating the cornerpoints into the trendline KOS
    global tri_corner_points_trendline_KOS

    tri_corner_points_global_KOS = calc_tri_corner_points(triangle_vectors_of_stl)
    tri_corner_points_trendline_KOS = translate_and_rotate_points_from_OLD_to_NEW_KOS(tri_corner_points_global_KOS, trendline_KOS_in_global_KOS,
                                                             center_of_pointcloud_weighted_in_global_KOS)
    points_x_y_trendline_KOS = tri_corner_points_trendline_KOS[:, 0:2]
    points_z_trendline_KOS = tri_corner_points_trendline_KOS[:, 2]

    # Creat grid
    global max_x,min_x,max_y,min_y,grid_x,dx,dy,y_0_grid_point_index,x_0_grid_point_index
    grid_resolution_j = grid_resolution * 1j
    max_x = max(points_x_y_trendline_KOS[:, 0])
    min_x = min(points_x_y_trendline_KOS[:, 0])
    max_y = max(points_x_y_trendline_KOS[:, 1])
    min_y = min(points_x_y_trendline_KOS[:, 1])
    dy = (max_y - min_y) / grid_resolution
    dx = (max_x - min_x) / grid_resolution
    y_0_grid_point_index = np.asarray(np.round(max_y / (max_y - min_y) * grid_resolution), dtype=np.int32)
    x_0_grid_point_index = np.asarray(np.round(max_x / (max_x - min_x) * grid_resolution), dtype=np.int32)
    
    grid_x, grid_y = np.mgrid[min_x:max_x:grid_resolution_j, min_y:max_y:grid_resolution_j]

    # Interpolating  the Surface Geometry
    z_grid_values_linear_trendline_KOS = griddata(np.asarray(points_x_y_trendline_KOS, dtype=np.float32), np.asarray(points_z_trendline_KOS, dtype=np.float32), (grid_x, grid_y),
                                    method='linear')

    return z_grid_values_linear_trendline_KOS
def calc_points_in_trendline_KOS_for_interpolation():
    tri_corner_points = calc_tri_corner_points(triangle_vectors_of_stl)
    points = translate_and_rotate_points_from_OLD_to_NEW_KOS(tri_corner_points, trendline_KOS_in_global_KOS, center_of_pointcloud_weighted_in_global_KOS)

    return points
def show_interpolation_and_draw_start_end_points(z_grid_values_linear_trendline_KOS):
    # Show interpolation
    print("Please close window after selecting Start- and Endpoint")
    figure = pyplot.figure()  # Comment_DB: create a new figure
    plt.imshow(z_grid_values_linear_trendline_KOS.T, extent=(min_x, max_x, min_y, max_y), origin='lower')
    plt.colorbar()
    # plt.plot(x_values_trim, y_values, 'bo', linewidth=2.0, label='Schnitt')
    plt.title('Please select Start- and Endpoint')

    # Draw Start and Endpoint
    global x_start_end_point, y_start_end_point
    x_start_end_point = []
    y_start_end_point = []

    def onclick(event):
        print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              (event.button, event.x, event.y, event.xdata, event.ydata))
        plt.plot(event.xdata, event.ydata, ',')
        figure.canvas.draw()
        global x_start_end_point, y_start_end_point

        if len(x_start_end_point)<2:
            x_start_end_point.append(event.xdata)
            y_start_end_point.append(event.ydata)
        else:
            print("That was one point to much. Please select again Start and Endpoint")
            x_start_end_point.clear()
            y_start_end_point.clear()
        if len(x_start_end_point) == 2:
            if x_start_end_point[0]>x_start_end_point[1]:
                xdata_1 = x_start_end_point[1]
                ydata_1 = y_start_end_point[1]

                x_start_end_point[1] = x_start_end_point[0]
                y_start_end_point[1] = y_start_end_point[0]

                x_start_end_point[0] = xdata_1
                y_start_end_point[0] = ydata_1

    figure.canvas.mpl_connect('button_press_event', onclick)
    global continue_bool
    continue_bool = False

    def handle_close(event):
        print('Closed Figure!')
        global continue_bool
        continue_bool = True

    figure.canvas.mpl_connect('close_event', handle_close)
    plt.show()

    # Wait until points are drawn in and the window is closed
    while continue_bool is False:
        pyplot.pause(2)

# Not used in stl_preprocessing, but needed if we load settings (see load settingssheet)
def calc_trendline_and_interpolate(grid_resolution, input_file):
    # Read in the file:
    calc_trendline_KOS_of_geometry_from_stl_file(input_file)  # saved in global var!
    # Interpolation:

    z_grid_values_linear_trendline_KOS = interpolate_geometrie(grid_resolution)

    return grid_x, max_x, max_y, min_x, min_y, z_grid_values_linear_trendline_KOS,x_0_grid_point_index, y_0_grid_point_index

#######################################################################################################################
# Translation and rotation from Points to new trendline_axis_KOS
def translate_and_rotate_points_from_OLD_to_NEW_KOS(points_in_old_KOS, new_KOS_in_old_KOS,
                                                    new_zero_point_in_old_KOS, reverse=False):
    # Idea: First rotating trendline_axis_x to (1,0,0), then shifting points to new center of coordinat system(center point weighted)

    # Basic Coordinate System
    x_axis_old_KOS = np.asarray((1, 0, 0), dtype=np.float32)
    y_axis_old_KOS = np.asarray((0, 1, 0), dtype=np.float32)
    z_axis_old_KOS = np.asarray((0, 0, 1), dtype=np.float32)

    # rotation angel
    anglez, angley = calc_angle_for_coordinate_rotation_x_trendline(x_axis_old_KOS, z_axis_old_KOS, new_KOS_in_old_KOS) # x-Axis now: (1,0,0)
    anglex = correction_angle_x_axis_for_y_z_orientation(angley, anglez, new_KOS_in_old_KOS, y_axis_old_KOS, z_axis_old_KOS)  # for y(0,1,0) and z(0,0,1)
    if reverse:
        anglez=-anglez
        angley=-angley
        anglex=-anglex

    # Rotation of translation vector/new zero
    translation_vector = calc_translation_vector(anglex, angley, anglez, new_zero_point_in_old_KOS, reverse,
                                                 x_axis_old_KOS, y_axis_old_KOS, z_axis_old_KOS)

    # Rotation of points
    points_in_new_KOS = rotate_points(anglex, angley, anglez, points_in_old_KOS, reverse, x_axis_old_KOS, y_axis_old_KOS, z_axis_old_KOS)

    # Translate rotated points to new zero
    points_in_new_KOS = translate_points_to_new_zero(translation_vector, points_in_new_KOS)

    return points_in_new_KOS


def calc_translation_vector(anglex, angley, anglez, new_zero_point_in_old_KOS, reverse, x_axis_old_KOS, y_axis_old_KOS, z_axis_old_KOS):
    if reverse:
        translation_vector = -(
            new_zero_point_in_old_KOS)  # At the reverse case we have to shift in the old unrotated KOS
    else:
        translation_vector = rotate_point_around_z_y_and_x_axis_with_given_angle(anglez, angley, anglex, z_axis_old_KOS,
                                                                                 y_axis_old_KOS, x_axis_old_KOS,
                                                                                 new_zero_point_in_old_KOS,
                                                                                 reverse)
    return translation_vector
def calc_angle_for_coordinate_rotation_x_trendline(x_axis_old_KOS, z_axis_old_KOS, new_KOS_in_old_KOS):

    # Calculation of angel z
    new_x_trendline_projected_to_x_y = project_pointtoplane((new_KOS_in_old_KOS[0][:]), z_axis_old_KOS, np.zeros(3))
    new_x_trendline_projected_to_x_y = norm_vector(new_x_trendline_projected_to_x_y)

    anglez = math.acos(np.dot(x_axis_old_KOS, new_x_trendline_projected_to_x_y))
    # if y negativ, x_rotation in other direction
    if new_KOS_in_old_KOS[0][1] <= -0.0001: anglez = -anglez


    # Calculation of angel y
    rotated_x_trend_around_z = Quaternion(axis=z_axis_old_KOS, angle=-anglez).rotate(new_KOS_in_old_KOS[0][:])
    rotated_x_trend_around_z = norm_vector(rotated_x_trend_around_z)

    angley = math.acos(np.dot(x_axis_old_KOS, rotated_x_trend_around_z))
    # if z negativ, x_rotation in other direction
    if rotated_x_trend_around_z[2] <= -0.0001: angley = -angley

    return anglez, angley
def correction_angle_x_axis_for_y_z_orientation(angley, anglez, new_trendline_axis_in_old_KOS, y_axis_old_KOS, z_axis_old_KOS):
    # Rotate trendline_axis around y and z
    new_trendline_axis_in_global_KOS = []
    for i in range(len(new_trendline_axis_in_old_KOS[:, 0])):
        new_trendline_axis_points_rotatet_i = rotate_point_around_z_y_and_x_axis_with_given_angle(anglez, angley, 0, z_axis_old_KOS,
                                                                                                  y_axis_old_KOS, y_axis_old_KOS,
                                                                                                  new_trendline_axis_in_old_KOS[
                                                                                                    i], False)
        new_trendline_axis_in_global_KOS.append(new_trendline_axis_points_rotatet_i)
    angel_x = 0
    if new_trendline_axis_in_global_KOS[2][1]>0 and new_trendline_axis_in_global_KOS[1][1]>0:
        try: angel_x = -math.acos(np.dot(z_axis_old_KOS, new_trendline_axis_in_global_KOS[2]))
        except: angel_x = 0
    elif new_trendline_axis_in_global_KOS[2][1] < 0 and new_trendline_axis_in_global_KOS[1][1] > 0:
        try:angel_x = math.acos(np.dot(z_axis_old_KOS, new_trendline_axis_in_global_KOS[2]))
        except: angel_x = 0

    elif new_trendline_axis_in_global_KOS[2][1] < 0 and new_trendline_axis_in_global_KOS[1][1] < 0:
        try:angel_x = math.pi - math.acos(np.dot(z_axis_old_KOS, new_trendline_axis_in_global_KOS[2]))
        except: angel_x = 0

    elif new_trendline_axis_in_global_KOS[2][1] > 0 and new_trendline_axis_in_global_KOS[1][1] < 0:
        try:angel_x = math.pi + math.acos(np.dot(z_axis_old_KOS, new_trendline_axis_in_global_KOS[2]))
        except:   angel_x = 0


    return angel_x
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
def rotate_points(anglex, angley, anglez, points_in_old_KOS, reverse, x_axis_old_KOS, y_axis_old_KOS, z_axis_old_KOS):
    points_in_new_KOS = []
    for i in range(len(points_in_old_KOS[:, 0])):
        old_points_rotatet_i = rotate_point_around_z_y_and_x_axis_with_given_angle(anglez, angley, anglex, z_axis_old_KOS,
                                                                                   y_axis_old_KOS, x_axis_old_KOS,
                                                                                   points_in_old_KOS[i], reverse)
        points_in_new_KOS.append(old_points_rotatet_i)
    points_in_new_KOS = np.asarray(points_in_new_KOS, dtype=np.float32)
    return points_in_new_KOS
def translate_points_to_new_zero(new_zero_point_in_old_KOS_rotated, points_in_new_KOS):
    points_in_new_KOS[:, 0] = np.subtract(points_in_new_KOS[:, 0], new_zero_point_in_old_KOS_rotated[0])
    points_in_new_KOS[:, 1] = np.subtract(points_in_new_KOS[:, 1], new_zero_point_in_old_KOS_rotated[1])
    points_in_new_KOS[:, 2] = np.subtract(points_in_new_KOS[:, 2], new_zero_point_in_old_KOS_rotated[2])
    return points_in_new_KOS

#######################################################################################################################
# Calculation of the bending parameters
def calc_bending_parameters(grid_resolution_int, max_distance, width_for_edge_detection, z_grid_values_linear_trendline_KOS,
                            calc_2D_with_edge_detection, calc_3D_Solution):

    initialize_global_lists_of_3D_bending_and_plot_parameter()

    # Calc 2D-Bendpoints
    calc_bending_points(grid_resolution_int, z_grid_values_linear_trendline_KOS,
                        [x_start_end_point_list[-2], x_start_end_point_list[-1]],
                        [y_start_end_point_list[-2], y_start_end_point_list[-1]], max_distance,
                        width_for_edge_detection, 0, 0, True, calc_2D_with_edge_detection, calc_3D_Solution)

    if calc_3D_Solution:
        # Start calculating bendingpoints
        calculate_iteratively_local_bendingpoints(alpha_angle_list_3D, grid_resolution_int, max_distance,
                                                  width_for_edge_detection, x_values_trim_trendline_KOS_stacked,
                                                  x_start_end_point_list, y_values_trim_trendline_KOS_stacked,
                                                  y_start_end_point_list, z_grid_values_linear_trendline_KOS)
        # Show results
        show_results_2D_Plots_and_Colormap(surface_points_local_KOS_left_stacked,
                                           surface_points_local_KOS_right_stacked,
                                           surface_points_local_KOS_stacked, x_values_trim_trendline_KOS_stacked,
                                           y_values_trim_trendline_KOS_stacked, z_grid_values_linear_trendline_KOS)
# Functions in calc_bending_parameters
def initialize_global_lists_of_3D_bending_and_plot_parameter():
    #  Startpoint and one additional Point for the direction (Endpoint). First values are drawn into the colormap.
    global x_start_end_point_list, y_start_end_point_list, new_start_left, new_start_right
    x_start_end_point_list = x_start_end_point
    y_start_end_point_list = y_start_end_point
    new_start_left, new_start_right = [],[]
    # List of Bendingparameters
    global x_direction_list_global_KOS, x_direction_rotated_list_global_KOS, normal_patch_global_KOS, beta_angle_between_planes_list_3D, alpha_angle_list_3D, length_list_3D, edge_line_global, bend_pts_xyz_global_3D, normal_direction_start_3D, gamma_start
    x_direction_list_global_KOS, x_direction_rotated_list_global_KOS, normal_patch_global_KOS, beta_angle_between_planes_list_3D, alpha_angle_list_3D, edge_line_global, length_list_3D, bend_pts_xyz_global_3D, normal_direction_start_3D,gamma_start = [], [], [], [], [], [], [], [], [], 0

    global length_list_2DE, beta_angle_between_planes_list_2DE, alpha_angle_list_2DE
    length_list_2DE, beta_angle_between_planes_list_2DE, alpha_angle_list_2DE = [], [], []

    # Blue Points for 3D-global Plot
    global surfacepoints_between_Start_and_End, surface_points_global_KOS_stacked, surface_points_global_KOS_left_stacked, surface_points_global_KOS_right_stacked
    surfacepoints_between_Start_and_End, surface_points_global_KOS_stacked, surface_points_global_KOS_left_stacked, surface_points_global_KOS_right_stacked = [], [], [], []

    # 2D Plot tildet local KOS
    global surface_points_local_KOS_stacked, surface_points_local_KOS_left_stacked, surface_points_local_KOS_right_stacked, \
        bend_pts_xz_local_stacked, bend_pts_xz_local_right_stacked, bend_pts_xz_local_left_stacked
    surface_points_local_KOS_stacked, surface_points_local_KOS_left_stacked, surface_points_local_KOS_right_stacked, \
    bend_pts_xz_local_stacked, bend_pts_xz_local_right_stacked, bend_pts_xz_local_left_stacked = [], [], [], [], [], []

    # Line in 2D Colormap
    global x_values_trim_trendline_KOS_stacked, y_values_trim_trendline_KOS_stacked
    x_values_trim_trendline_KOS_stacked, y_values_trim_trendline_KOS_stacked = [], []


    global counter_failed_matches_of_edges
    counter_failed_matches_of_edges = 0

def calculate_iteratively_local_bendingpoints(alpha_angle_list, grid_resolution_int, max_distance,
                                              width_for_edge_detection, x_values_trim_stacked, xdata_list,
                                              y_values_trim_stacked, ydata_list, z_grid_values_linear_trendline_KOS):
    # calc_local_bending_points has no return value except num_bendpoints. All the bending parameters are saved in the global lists
    num_bendpoints = calc_bending_points(grid_resolution_int, z_grid_values_linear_trendline_KOS, [xdata_list[-2], xdata_list[-1]],
                                         [ydata_list[-2], ydata_list[-1]], max_distance, width_for_edge_detection, 0, 0)

    # Show calculated direction in colormap
    pyplot.figure()
    plt.imshow(z_grid_values_linear_trendline_KOS.T, extent=(min_x, max_x, min_y, max_y), origin='lower')
    plt.plot(x_values_trim_stacked[0][:], y_values_trim_stacked[0][:], 'bo', linewidth=1.0, label='Schnitt')
    plt.show(block=False)
    pyplot.pause(0.3) #Give system time to print the plot

    while num_bendpoints > 2:
        num_bendpoints = calc_bending_points(grid_resolution_int, z_grid_values_linear_trendline_KOS,
                                             [xdata_list[-2], xdata_list[-1]], [ydata_list[-2], ydata_list[-1]],
                                             max_distance, width_for_edge_detection, alpha_angle_list[-1], alpha_end=0)

        # Update colormap
        for i in range(len(x_values_trim_stacked)):
            plt.plot(x_values_trim_stacked[i][:], y_values_trim_stacked[i][:], 'bo', linewidth=1.0, label='Schnitt')
        plt.draw()
        pyplot.pause(0.01)
def show_results_2D_Plots_and_Colormap(surface_points_local_KOS_left_stacked,
                                       surface_points_local_KOS_right_stacked,
                                       surface_points_local_KOS_stacked, x_values_trim_stacked,
                                       y_values_trim_stacked, z_grid_values_linear_trendline_KOS):

    # In the different local KOS the x-Values don´t allign. Correcting it here for nicer 2Dplots
    connect_points_in_local_KOS(surface_points_local_KOS_stacked)
    connect_points_in_local_KOS(surface_points_local_KOS_left_stacked)
    connect_points_in_local_KOS(surface_points_local_KOS_right_stacked)
    connect_points_in_local_KOS(bend_pts_xz_local_stacked)
    connect_points_in_local_KOS(bend_pts_xz_local_left_stacked)
    connect_points_in_local_KOS(bend_pts_xz_local_right_stacked)

    #Colormap
    pyplot.figure()
    plt.title('Interpolation with start- end- connection')
    plt.imshow(z_grid_values_linear_trendline_KOS.T, extent=(min_x, max_x, min_y, max_y), origin='lower')
    plt.plot(x_values_trim_stacked[0][:], y_values_trim_stacked[0][:], 'bo', linewidth=1.0,
             label='Tape direction')
    plt.legend()
    for i in range(1,len(x_values_trim_stacked)):
        plt.plot(x_values_trim_stacked[i][:], y_values_trim_stacked[i][:], 'bo', linewidth=1.0)


    #2D-Sideview
    pyplot.figure()
    plt.title('Sideview height profil')
    plt.plot(surface_points_local_KOS_stacked[0][:, 0],
             surface_points_local_KOS_stacked[0][:, 2], 'bo', linewidth=1.0, label='surface points')
    plt.plot(bend_pts_xz_local_stacked[0][:, 0], bend_pts_xz_local_stacked[0][:, 1], color='green', linewidth=3.0,label='linear Approximation')
    plt.legend()
    for i in range(1,len(x_values_trim_stacked)):
        plt.plot(surface_points_local_KOS_stacked[i][:, 0],
                 surface_points_local_KOS_stacked[i][:, 2], 'bo', linewidth=1.0)
        plt.plot(bend_pts_xz_local_stacked[i][:, 0], bend_pts_xz_local_stacked[i][:, 1], color='green', linewidth=3.0)


    pyplot.figure()
    plt.subplot(211)
    plt.title('right')
    for i in range(len(x_values_trim_stacked)):
        plt.plot(surface_points_local_KOS_left_stacked[i][:, 0],
                 surface_points_local_KOS_left_stacked[i][:, 2], 'bo', linewidth=1.0, label='cross section')
        plt.plot(bend_pts_xz_local_left_stacked[i][:, 0], bend_pts_xz_local_left_stacked[i][:, 1], color='green', linewidth=3.0,
             label='linear Approximation')
    plt.subplot(212)
    plt.title('left')
    for i in range(len(x_values_trim_stacked)):
        plt.plot(surface_points_local_KOS_right_stacked[i][:, 0],
                 surface_points_local_KOS_right_stacked[i][:, 2], 'bo', linewidth=1.0, label='cross section')
        plt.plot(bend_pts_xz_local_right_stacked[i][:, 0], bend_pts_xz_local_right_stacked[i][:, 1], color='green', linewidth=3.0,
             label='linear Approximation')
    plt.show()
def connect_points_in_local_KOS(surface_points_local_KOS_stacked):
    for i in range(1, len(surface_points_local_KOS_stacked)):
        startx_i = surface_points_local_KOS_stacked[i][0, 0]
        endx_i_1 = surface_points_local_KOS_stacked[i - 1][-1, 0]
        difx = startx_i - endx_i_1
        surface_points_local_KOS_stacked[i][:, 0] = np.subtract(
            surface_points_local_KOS_stacked[i][:, 0], difx)

#######################################################################################################################
# Calculation of the 2D-Solution.
def calc_tape_parameter_for_2D_solution(bend_pts_xyz_global,bend_pts_xz_local):
    global bend_pts_xyz_global_2D, beta_angle_list_2D, lenght_list_2D, alpha_angle_list_2D
    bend_pts_xyz_global_2D = bend_pts_xyz_global
    beta_angle_list_2D = calc_2D_betas(bend_pts_xz_local)
    lenght_list_2D = calc_2D_lengths(bend_pts_xz_local)
    alpha_angle_list_2D = np.zeros(len(bend_pts_xz_local))

# beta-list for Chromo
def calc_2D_betas(bend_pts_xz_local):
    beta_list = []
    for i in range(1, len(bend_pts_xz_local) - 1):
        r0 = bend_pts_xz_local[i] - bend_pts_xz_local[i - 1]
        r1 = bend_pts_xz_local[i + 1] - bend_pts_xz_local[i]
        try: angle = math.acos(np.dot(r0, r1) / (np.linalg.norm(r0) * np.linalg.norm(r1)))
        except: angle = 0
        steepness_r0 = r0[1] / r0[0]
        steepness_r1 = r1[1] / r1[0]
        if steepness_r1 < steepness_r0:
            angle = -angle
        if np.abs(angle) < 0.01:  # Comment_DB: Small angles approx. to 0 degrees
            angle = 0
        beta_list.append(-angle)
    beta_list = np.asarray(beta_list)
    return beta_list
def calc_2D_lengths(bend_pts_xz_local):
    # l-list for Chromo
    l_list = []
    for i in range(1, len(bend_pts_xz_local)):
        l = np.linalg.norm(bend_pts_xz_local[i] - bend_pts_xz_local[i - 1])
        l_list.append(l)
    l_list = np.asarray(l_list)
    return l_list

#######################################################################################################################
# For every iteration in calc_bending_parameters new bendpoints and parameters have to be calculated.
def calc_bending_points(grid_resolution_int, z_grid_values_linear_trendline_KOS, x_start_end_point_trendline_KOS,
                        y_start_end_point_trendline_KOS, max_distance, width_for_edge_detection, alpha_start=0,
                        alpha_end=0, calc_tape_para_2D=False, calc_2D_with_edge_detection=False, calc_3D_Solution=True):
    


    # x_y_z_Axis of new KOS
    local_KOS_in_trendline_KOS, x_slope = calc_local_KOS_in_trendline_KOS(x_start_end_point_trendline_KOS,y_start_end_point_trendline_KOS)


    # Calc bendpoints on surface in new trendline direction, including left and right for edge directions
    bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right, bend_pts_xyz_trendline, \
    bend_pts_xz_local,local_KOS_in_global_KOS = calc_bend_pts_between_start_and_end(alpha_end, alpha_start,
                                                                                    grid_resolution_int, max_distance,
                                                                                    local_KOS_in_trendline_KOS,
                                                                                    width_for_edge_detection, x_slope,
                                                                                    x_start_end_point_trendline_KOS,
                                                                                    y_start_end_point_trendline_KOS,
                                                                                    z_grid_values_linear_trendline_KOS,
                                                                                    calc_tape_para_2D)

    normal_at_start, x_direction_start = calc_start_vectors_from_bendpoints(bend_pts_xyz_global,
                                                                            bend_pts_xyz_global_left,
                                                                            bend_pts_xyz_global_right,
                                                                            bend_pts_xz_local, local_KOS_in_global_KOS)

    # Calc Tapeparameters 2D
    if calc_tape_para_2D:
        calc_tape_parameter_for_2D_solution(bend_pts_xyz_global,bend_pts_xz_local)

    # Calc Tapeparameters 2DE/3D
    edge_directions = calc_edge_directions(bend_pts_xyz_global_left, bend_pts_xyz_global_right)
    if calc_tape_para_2D:
        global edge_directions_global_2DE
        edge_directions_global_2DE = edge_directions

    x_direction_list_current_direction_global_KOS, \
    normal_patch_current_direction_global_KOS, \
    rotated_x_direction_around_edge_current_direction_global_KOS, \
    beta_angle_between_planes_list_current_direction,\
    alpha_angle_list_current_direction,\
    lengths_between_planes_list = calc_tape_parameter_2DE_3D(bend_pts_xyz_global, bend_pts_xyz_global_left,
                                                             bend_pts_xyz_global_right,
                                                             bend_pts_xyz_trendline, edge_directions,
                                                             x_direction_start, normal_at_start,
                                                             bend_pts_xz_local, calc_2D_with_edge_detection)

    # Add Bend/Tapeparameters from first found bendpoint to global list
    append_bend_parameters_at_first_bendpoint_to_global_list(alpha_angle_list_current_direction, bend_pts_xyz_global,
                                                             bend_pts_xyz_trendline,
                                                             beta_angle_between_planes_list_current_direction,
                                                             edge_directions, lengths_between_planes_list,
                                                             normal_patch_current_direction_global_KOS,
                                                             rotated_x_direction_around_edge_current_direction_global_KOS,
                                                             x_direction_list_current_direction_global_KOS, calc_tape_para_2D, calc_3D_Solution)

    return len(bend_pts_xyz_global_left)


def calc_start_vectors_from_bendpoints(bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right,
                                       bend_pts_xz_local, local_KOS_in_global_KOS):
    # Lokal x_direction and normal
    x_direction_start = norm_vector(bend_pts_xyz_global[1] - bend_pts_xyz_global[0])  # Start_direction
    # L_aim, Start Direction and Start Normal
    global surfacepoints_between_Start_and_End, normal_direction_start, normal_direction_start_3D, gamma_start
    if surfacepoints_between_Start_and_End == []:  # Just once from the first calculation cycle
        # Start_End_connection for L_aim
        surfacepoints_between_Start_and_End = bend_pts_xz_local
        # Start normal
        normal_direction_start = np.cross(x_direction_start, local_KOS_in_global_KOS[1])
        normal_direction_start = norm_vector(normal_direction_start)


        normal_at_start = calc_tape_normal(bend_pts_xyz_global[1], bend_pts_xyz_global_left[0],
                                           bend_pts_xyz_global_right[0])  # Start_normal #normal_patch_global_KOS[-1]

        if normal_direction_start_3D == []:
            normal_direction_start_3D = normal_at_start

            gamma_start = math.acos(np.dot(normal_direction_start_3D, normal_direction_start))

            normals_rotated = translate_and_rotate_points_from_OLD_to_NEW_KOS(
                np.stack([x_direction_start, normal_direction_start_3D, normal_direction_start]),
                local_KOS_in_global_KOS, np.asarray([0, 0, 0]))

            if normals_rotated[1][1] > 0:
                gamma_start = -gamma_start


    else:
        normal_at_start = calc_tape_normal(bend_pts_xyz_global[1], bend_pts_xyz_global_left[0],
                                           bend_pts_xyz_global_right[0])  # Start_normal #normal_patch_global_KOS[-1]

    return normal_at_start, x_direction_start


def calc_bend_pts_between_start_and_end(alpha_end, alpha_start, grid_resolution_int, max_distance,
                                        local_KOS_in_trendline_KOS, width_for_edge_detection, x_slope,
                                        x_start_end_point_trendline_KOS, y_start_end_point_trendline_KOS,
                                        z_grid_values_linear_trendline_KOS, calc_tape_para_2D):

    start_point_xyz_trendline_KOS, start_point_xyz_trendline_KOS_left, start_point_xyz_trendline_KOS_right, \
    x_end_index, x_end_index_left, x_end_index_right, \
    x_start_index, x_start_index_left, x_start_index_right = calc_start_end_point_trendline_KOS(
        alpha_end, alpha_start, grid_resolution_int, local_KOS_in_trendline_KOS, width_for_edge_detection,
        x_start_end_point_trendline_KOS, y_start_end_point_trendline_KOS, z_grid_values_linear_trendline_KOS)

    # Calculation of the surfacepoints in the local, trendline and global direction and extracting bendpoints from them.

    # Left
    bend_pts_xyz_global_left, bend_pts_xyz_trendline_left, bend_pts_xz_local_left, surface_points_global_KOS_left, \
    surface_points_local_KOS_left, x_values_trim_left_trendline_KOS, y_values_trim_left_trendline_KOS, local_KOS_in_global_KOS_left = calc_points_on_surface_and_extract_bendline(
        grid_resolution_int, max_distance, start_point_xyz_trendline_KOS_left, local_KOS_in_trendline_KOS,
        x_end_index_left, x_slope, x_start_index_left, z_grid_values_linear_trendline_KOS)
    # Right
    bend_pts_xyz_global_right, bend_pts_xyz_trendline_right, bend_pts_xz_local_right, surface_points_global_KOS_right, \
    surface_points_local_KOS_right, x_values_trim_right_trendline_KOS, y_values_trim_right_trendline_KOS,local_KOS_in_global_KOS_right = calc_points_on_surface_and_extract_bendline(
        grid_resolution_int, max_distance, start_point_xyz_trendline_KOS_right, local_KOS_in_trendline_KOS,
        x_end_index_right, x_slope, x_start_index_right, z_grid_values_linear_trendline_KOS)
    # Center
    bend_pts_xyz_global, bend_pts_xyz_trendline, bend_pts_xz_local, surface_points_global_KOS, \
    surface_points_local_KOS, x_values_trim_trendline_KOS, y_values_trim_trendline_KOS,local_KOS_in_global_KOS = calc_points_on_surface_and_extract_bendline(
        grid_resolution_int, max_distance, start_point_xyz_trendline_KOS, local_KOS_in_trendline_KOS,
        x_end_index, x_slope, x_start_index, z_grid_values_linear_trendline_KOS)

    if calc_tape_para_2D:
        global surface_points_start_end_global_KOS, surface_points_start_end_global_KOS_left, surface_points_start_end_global_KOS_right
        surface_points_start_end_global_KOS = surface_points_global_KOS
        surface_points_start_end_global_KOS_left = surface_points_global_KOS_left
        surface_points_start_end_global_KOS_right = surface_points_global_KOS_right

        global bend_points_start_end_global_KOS, bend_points_start_end_global_KOS_left, bend_points_start_end_global_KOS_right
        bend_points_start_end_global_KOS = bend_pts_xyz_global
        bend_points_start_end_global_KOS_left = bend_pts_xyz_global_left
        bend_points_start_end_global_KOS_right = bend_pts_xyz_global_right
    # Trim plot points to second bendpoint and add to global list
    if not calc_tape_para_2D:
        append_plot_points_till_second_bendpoint_to_global_list(bend_pts_xz_local, bend_pts_xz_local_left,
                                                            bend_pts_xz_local_right,
                                                            surface_points_global_KOS,
                                                            surface_points_global_KOS_left,
                                                            surface_points_global_KOS_right,
                                                            surface_points_local_KOS,
                                                            surface_points_local_KOS_left,
                                                            surface_points_local_KOS_right,
                                                            x_values_trim_trendline_KOS, y_values_trim_trendline_KOS)

    global new_start_left, new_start_right
    if not calc_tape_para_2D:
        new_start_right = bend_pts_xyz_trendline_right[1]
        new_start_left = bend_pts_xyz_trendline_left[1]

    return bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right, bend_pts_xyz_trendline, bend_pts_xz_local,local_KOS_in_global_KOS


def calc_start_end_point_trendline_KOS(alpha_end, alpha_start, grid_resolution_int, local_KOS_in_trendline_KOS,
                                       width_for_edge_detection, x_start_end_point_trendline_KOS, y_start_end_point_trendline_KOS,
                                       z_grid_values_linear_trendline_KOS):
    # Start and endpoint for tape section
    end_point_xyz_trendline_KOS, \
    start_point_xyz_trendline_KOS, \
    x_start_index, \
    x_end_index = calc_Start_End_in_trendline_KOS_from_xdata_ydata(grid_resolution_int, x_start_end_point_trendline_KOS,
                                                                   y_start_end_point_trendline_KOS,
                                                                   z_grid_values_linear_trendline_KOS)
    # Start and endpoint from side line for estimating bending angels
    # Cornerpoints
    delta_length_start_bend = calc_delta_length_at_bend(width_for_edge_detection, alpha_start)
    delta_length_end_bend = calc_delta_length_at_bend(width_for_edge_detection, alpha_end)
    end_point_xyz_trendline_KOS_right, \
    start_point_xyz_trendline_KOS_right, \
    x_end_index_right, \
    x_start_index_right = calc_start_end_point_side_in_trendline_KOS(False, delta_length_end_bend,
                                                                     delta_length_start_bend,
                                                                     end_point_xyz_trendline_KOS, grid_resolution_int,
                                                                     start_point_xyz_trendline_KOS,
                                                                     width_for_edge_detection,
                                                                     local_KOS_in_trendline_KOS,
                                                                     z_grid_values_linear_trendline_KOS)
    end_point_xyz_trendline_KOS_left, \
    start_point_xyz_trendline_KOS_left, \
    x_end_index_left, \
    x_start_index_left = calc_start_end_point_side_in_trendline_KOS(True, delta_length_end_bend,
                                                                    delta_length_start_bend,
                                                                    end_point_xyz_trendline_KOS, grid_resolution_int,
                                                                    start_point_xyz_trendline_KOS,
                                                                    width_for_edge_detection,
                                                                    local_KOS_in_trendline_KOS,
                                                                    z_grid_values_linear_trendline_KOS)
    return start_point_xyz_trendline_KOS, start_point_xyz_trendline_KOS_left, start_point_xyz_trendline_KOS_right, x_end_index, x_end_index_left, x_end_index_right, x_start_index, x_start_index_left, x_start_index_right


def calc_points_on_surface_and_extract_bendline(grid_resolution_int, max_distance, start_point_xyz_trendline_KOS,
                                                local_KOS_in_trendline_KOS, x_end_index, x_slope, x_start_index,
                                                z_grid_values_linear_trendline_KOS):

    new_bending_direction_points_on_surface_trendline_KOS, y_intercept_trendline_KOS, x_values_indizes_trim, \
    x_values_trim_trendline_KOS, y_values_indizes_trim, y_values_trim_trendline_KOS = extract_points_from_interpolated_surface(
        start_point_xyz_trendline_KOS, grid_resolution_int, x_slope, z_grid_values_linear_trendline_KOS, x_start_index, x_end_index)


    surface_points_local_KOS, surface_points_global_KOS, \
    local_KOS_in_global_KOS = transform_surface_points_to_global_and_local_KOS(y_intercept_trendline_KOS,
                                                                               new_bending_direction_points_on_surface_trendline_KOS,
                                                                               local_KOS_in_trendline_KOS)
    bend_pts_xz_local, bend_pts_xyz_trendline, bend_pts_xyz_global = calc_local_and_global_bendpoints(max_distance,
                                                                                                      surface_points_local_KOS,
                                                                                                      local_KOS_in_trendline_KOS,
                                                                                                      y_intercept_trendline_KOS)
    return bend_pts_xyz_global, bend_pts_xyz_trendline, bend_pts_xz_local, surface_points_global_KOS, surface_points_local_KOS, x_values_trim_trendline_KOS, y_values_trim_trendline_KOS,local_KOS_in_global_KOS
def append_bend_parameters_at_first_bendpoint_to_global_list(alpha_angle_list_current_direction, bend_pts_xyz_global,
                                                             bend_pts_xyz_trendline,
                                                             beta_angle_between_planes_list_current_direction,
                                                             edge_directions, lengths_between_planes_list,
                                                             normal_patch_current_direction_global_KOS,
                                                             rotated_x_direction_around_edge_current_direction_global_KOS,
                                                             x_direction_list_current_direction_global_KOS, calc_tape_para_2D, calc_3D_Solution):

    global x_direction_list_global_KOS, x_direction_rotated_list_global_KOS, normal_patch_global_KOS, \
        length_list_3D,beta_angle_between_planes_list_3D, alpha_angle_list_3D,  \
        length_list_2DE, beta_angle_between_planes_list_2DE, alpha_angle_list_2DE, \
        edge_line_global, x_start_end_point_list, y_start_end_point_list, bend_pts_xyz_global_3D

    # Tape parameter for EA
    if calc_tape_para_2D:
        length_list_2DE = np.asarray(lengths_between_planes_list)
        beta_angle_between_planes_list_2DE = np.asarray(
            beta_angle_between_planes_list_current_direction)
        alpha_angle_list_2DE = np.asarray(alpha_angle_list_current_direction)

        if calc_3D_Solution == False: # If no 3D Solution is calculated, we need to save the first x_direction and normal: start_d/start_n Parameter
            x_direction_list_global_KOS.append(x_direction_list_current_direction_global_KOS[0])
            x_direction_rotated_list_global_KOS.append(rotated_x_direction_around_edge_current_direction_global_KOS[0])
            normal_patch_global_KOS.append(normal_patch_current_direction_global_KOS[0])

    else:
        length_list_3D.append(lengths_between_planes_list[0])

        # Bendpoint and edge
        bend_pts_xyz_global_3D.append(bend_pts_xyz_global[0])

        edge_line_global_current_direction = [
            (bend_pts_xyz_global[:][:2] + np.multiply(edge_directions[:][:2], 50)),
            (bend_pts_xyz_global[:][:2] - np.multiply(edge_directions[:][:2], 50))]
        edge_line_global.append([edge_line_global_current_direction[0][1], edge_line_global_current_direction[1][1]])

        if len(edge_directions) == 2:  # Last part
            bend_pts_xyz_global_3D.append(bend_pts_xyz_global[1])
        if len(edge_directions) > 2:
            # startpoint and direction of next bending direction.
            x_start_end_point_list.append(bend_pts_xyz_trendline[1][0])
            y_start_end_point_list.append(bend_pts_xyz_trendline[1][1])

            # New Direction and Endpoint in trendline KOS(same then interpolation)
            rotated_x_direction_point_in_trendline_KOS = np.asarray(
                [bend_pts_xyz_global[1] + np.multiply(rotated_x_direction_around_edge_current_direction_global_KOS[0], 1),
                 bend_pts_xyz_global[1]])
            rotated_x_direction_around_edge_trendline_KOS_current_direction = translate_and_rotate_points_from_OLD_to_NEW_KOS(
                rotated_x_direction_point_in_trendline_KOS, trendline_KOS_in_global_KOS, center_of_pointcloud_weighted_in_global_KOS)


            x_intersect, y_intersect = calc_new_endpoint(bend_pts_xyz_trendline,
                                                         rotated_x_direction_around_edge_trendline_KOS_current_direction,
                                                         x_start_end_point_list, y_start_end_point_list)

            x_start_end_point_list.append(x_intersect)
            y_start_end_point_list.append(y_intersect)

            # Tape parameter for EA
            beta_angle_between_planes_list_3D.append(beta_angle_between_planes_list_current_direction[0])
            alpha_angle_list_3D.append(alpha_angle_list_current_direction[0])


            # Directions
            x_direction_list_global_KOS.append(x_direction_list_current_direction_global_KOS[0])
            x_direction_rotated_list_global_KOS.append(rotated_x_direction_around_edge_current_direction_global_KOS[0])
            normal_patch_global_KOS.append(normal_patch_current_direction_global_KOS[1])

def calc_new_endpoint(bend_pts_xyz_trendline, rotated_x_direction_around_edge_trendline_KOS_current_direction,
                      xdata_list, ydata_list):
    # new endpoint is the intersection between the line representing the the new direction and the end line
    # Line are represented by two points

    # End line: Line perpendicular to first start-end-direction, starting at endpoint
    drawn_endpoint = np.asarray([xdata_list[1], ydata_list[1]])

    direction_perpendicular_to_first_start_end = np.asarray([-(ydata_list[1] - ydata_list[0]), (xdata_list[1] - xdata_list[0])])
    point_in_direction_perpendicular_to_start_end_drawn = drawn_endpoint + direction_perpendicular_to_first_start_end

    # new startpoint
    new_startpoint = np.asarray([bend_pts_xyz_trendline[1][0], bend_pts_xyz_trendline[1][1]])
    other_point_in_new_direction = np.asarray(
        [rotated_x_direction_around_edge_trendline_KOS_current_direction[0][0],
         rotated_x_direction_around_edge_trendline_KOS_current_direction[0][1]])

    # Calculating the two lines
    stacked_points = np.vstack([new_startpoint, other_point_in_new_direction, drawn_endpoint,
                   point_in_direction_perpendicular_to_start_end_drawn])
    homogen_coordinates = np.hstack((stacked_points, np.ones((4, 1))))
    line_in_new_direction = np.cross(homogen_coordinates[0], homogen_coordinates[1])
    end_line = np.cross(homogen_coordinates[2], homogen_coordinates[3])
    x, y, z = np.cross(line_in_new_direction, end_line)  # point of intersection
    # if z == 0:  # lines are parallel
    x_intersect = x / z
    y_intersect = y / z
    return x_intersect, y_intersect
def append_plot_points_till_second_bendpoint_to_global_list(bend_pts_xz_local, bend_pts_xz_local_left,
                                                            bend_pts_xz_local_right,
                                                            surface_points_global_KOS,
                                                            surface_points_global_KOS_left,
                                                            surface_points_global_KOS_right,
                                                            surface_points_local_KOS,
                                                            surface_points_local_KOS_left,
                                                            surface_points_local_KOS_right,
                                                            x_values_trim_trendline_KOS, y_values_trim_trendline_KOS):
    # Calc trim index
    end_index_pts_global, end_index_pts_global_left, end_index_pts_global_right = calc_trim_index_at_second_bendpoint(
        bend_pts_xz_local, bend_pts_xz_local_left, bend_pts_xz_local_right, surface_points_local_KOS,
        surface_points_local_KOS_left, surface_points_local_KOS_right)



    # Global 3D-Plot
    global surface_points_global_KOS_stacked, surface_points_global_KOS_left_stacked, surface_points_global_KOS_right_stacked
    surface_points_global_KOS_stacked.append(
        surface_points_global_KOS[:][:end_index_pts_global])
    surface_points_global_KOS_right_stacked.append(
        surface_points_global_KOS_right[:][:end_index_pts_global_right])
    surface_points_global_KOS_left_stacked.append(
        surface_points_global_KOS_left[:][:end_index_pts_global_left])
    # local 2D-Plot
    global surface_points_local_KOS_stacked, surface_points_local_KOS_left_stacked, surface_points_local_KOS_right_stacked, \
        bend_pts_xz_local_stacked, bend_pts_xz_local_right_stacked, bend_pts_xz_local_left_stacked
    surface_points_local_KOS_stacked.append(
        surface_points_local_KOS[:end_index_pts_global])
    surface_points_local_KOS_left_stacked.append(
        surface_points_local_KOS_left[:end_index_pts_global_left])
    surface_points_local_KOS_right_stacked.append(
        surface_points_local_KOS_right[:end_index_pts_global_right])
    bend_pts_xz_local_stacked.append(bend_pts_xz_local[:2])
    bend_pts_xz_local_right_stacked.append(bend_pts_xz_local_right[:2])
    bend_pts_xz_local_left_stacked.append(bend_pts_xz_local_left[:2])
    # for trendline KOS, 2D-Plot-Colormap
    global x_values_trim_trendline_KOS_stacked, y_values_trim_trendline_KOS_stacked
    x_values_trim_trendline_KOS_stacked.append(x_values_trim_trendline_KOS[:end_index_pts_global])
    y_values_trim_trendline_KOS_stacked.append(y_values_trim_trendline_KOS[:end_index_pts_global])
def calc_trim_index_at_second_bendpoint(bend_pts_xz_local, bend_pts_xz_local_left, bend_pts_xz_local_right,
                                        surface_points_local_KOS,
                                        surface_points_local_KOS_left,
                                        surface_points_local_KOS_right):
    try:
        end_index_pts_global = \
        [i for i, x in enumerate(surface_points_local_KOS[:, 0]) if x >= bend_pts_xz_local[1][0]][0]
    except:
        end_index_pts_global = -1
    try:
        end_index_pts_global_left = [i for i, x in enumerate(surface_points_local_KOS_left[:, 0]) if
                                     x >= bend_pts_xz_local_left[1][0]][0]
    except:
        end_index_pts_global_left = -1
    try:
        end_index_pts_global_right = [i for i, x in enumerate(surface_points_local_KOS_right[:, 0]) if
                                      x >= bend_pts_xz_local_right[1][0]][0]
    except:
        end_index_pts_global_right = -1
    return end_index_pts_global, end_index_pts_global_left, end_index_pts_global_right
def calc_bend_pts(max_distance, x_z_surface_point):

    bend_pts_xz = []
    bend_pts_xz.append([x_z_surface_point[0][0], x_z_surface_point[0][1]])  # Comment_DB: start point 2D (x coord, y coord)
    bend_pts_xz.append([x_z_surface_point[-1][0],x_z_surface_point[-1][1]])  # Comment_DB: end point 2D (x coord, y coord)
    bend_pts_xz = np.asarray(bend_pts_xz)

    # Inserting bendpoints if max distance to surface to bigg:
    insert_pts = True
    while insert_pts:
        points_on_line_between_bends_filled_up = []
        points_on_line_between_bends_filled_up.append([bend_pts_xz[0][0], bend_pts_xz[0][1]])  # Comment_DB: only the first bend point (starting point at edge) appended to bend points curve list

        points_on_line_between_bends_filled_up = calc_points_on_line_between_bends_filled_up(bend_pts_xz,
                                                                                             points_on_line_between_bends_filled_up, x_z_surface_point)

        max_divergence = calc_point_of_max_divergence_between_smooth_and_lin_curve(
            points_on_line_between_bends_filled_up, x_z_surface_point)

        # Comment_DB: We know at which x-coord of points_on_line_between_bends_filled_up the max_divergence happens --> counter i



        # no further points, if the chosen maximum distance is not surpassed
        if max_divergence[0] < max_distance:  # Comment_DB: This implies that there will be one extra bend, as the above code will have executed already, max_distance: User Input
            break

        bend_pts_xz = np.insert(bend_pts_xz, -1,
                            np.array([points_on_line_between_bends_filled_up[max_divergence[1]][0],
                                      x_z_surface_point[max_divergence[1]][1]]),
                            axis=0)  # Comment_DB: insert a corner at x coord (counter i) and y coord (counter i) of max divergence
        bend_pts_xz = bend_pts_xz[bend_pts_xz[:, 0].argsort()]  # Comment_DB: Bend points sorted in an array


    return bend_pts_xz
def calc_point_of_max_divergence_between_smooth_and_lin_curve(points_on_line_between_bends_filled_up, x_y_points_filled_up):
    # Comment_DB: curve_divergence in terms of y-distance # Largest deviation from smoothed curve: (COMMENT_DB: By now all the points in the above (linear) line have been appended)
    curve_divergence_y = []
    for i in range(len(points_on_line_between_bends_filled_up)):
        curve_divergence_y.append([points_on_line_between_bends_filled_up[i][0], (
                (points_on_line_between_bends_filled_up[i][0] - x_y_points_filled_up[i][0]) ** 2 + (
                points_on_line_between_bends_filled_up[i][1] - x_y_points_filled_up[i][1]) ** 2) ** 0.5])  # Comment_DB: (x-coord vs. change in y-coord) take the x coord and y-distance between linear curve and sav-gol curve and append
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
def calc_tape_parameter_2DE_3D(bend_pts_xyz_global, bend_pts_xyz_global_left, bend_pts_xyz_global_right,
                               bend_pts_xyz_trendline, edge_directions, x_direction_start, normal_at_start,
                               bend_pts_xz_local, calc_2D_with_edge_detection):

    x_direction_list_global = [x_direction_start]
    normal_direction_list_global = [normal_at_start]
    lenght_between_first_two_bends = np.linalg.norm(bend_pts_xz_local[1] - bend_pts_xz_local[0])
    lengths_between_planes_list = [lenght_between_first_two_bends]
    rotated_x_direction_around_edge_global = []
    beta_angle_between_planes_list, alpha_angle_between_planes_list = [],[]


    range_till_end  = (len(edge_directions) - 1) if calc_2D_with_edge_detection else 2

    for i in range(1, range_till_end):
        try: length_current_direction = np.linalg.norm(bend_pts_xz_local[i+1] - bend_pts_xz_local[i])
        except: break

        # Calc x_y_and_normal direction of the Tape at each bending point
        x_direction_before_bend_global = norm_vector(bend_pts_xyz_global[i] - bend_pts_xyz_global[i - 1])

        # We get the normal of the plane defined by the right and left bendpoint of the bend and the middlepoint on the next/previous bend.
        normal_before_bend_global = -calc_tape_normal(bend_pts_xyz_global[i-1], bend_pts_xyz_global_left[i],bend_pts_xyz_global_right[i])

        y_direction_before_bend_global = np.cross(normal_before_bend_global, x_direction_before_bend_global)

        tape_orientation_before_bend_KOS_in_global_KOS = np.stack([x_direction_before_bend_global, y_direction_before_bend_global, normal_before_bend_global])

        try: normal_after_bend_global = calc_tape_normal(bend_pts_xyz_global[i + 1], bend_pts_xyz_global_left[i],bend_pts_xyz_global_right[i])
        except:
            print("Not same amount of bending points")
            break

        beta_angle = calc_3D_beta(bend_pts_xyz_trendline, i, normal_before_bend_global, normal_after_bend_global)

        alpha_angle = calc_3D_alpha(edge_directions, i, normal_before_bend_global,
                                    tape_orientation_before_bend_KOS_in_global_KOS,
                                    y_direction_before_bend_global)
        # Calc new x_direction after bending around edge
        rotated_x_direction_around_edge_i_global = Quaternion(axis=edge_directions[i], angle=beta_angle).rotate(x_direction_before_bend_global)
        rotated_x_direction_around_edge_global.append(rotated_x_direction_around_edge_i_global)



        # safe x_directions, alpha, beta, and normals in a List.
        x_direction_list_global.append(x_direction_before_bend_global)    #x_direction_start is two times in the list. If there are 1+ bending points.
        alpha_angle_between_planes_list.append(alpha_angle)
        beta_angle_between_planes_list.append(beta_angle)
        normal_direction_list_global.append(normal_after_bend_global)
        lengths_between_planes_list.append(length_current_direction)


    rotated_x_direction_around_edge_global = np.asarray(rotated_x_direction_around_edge_global)
    x_direction_list_global = np.asarray(x_direction_list_global)

    return x_direction_list_global, normal_direction_list_global, rotated_x_direction_around_edge_global,beta_angle_between_planes_list,alpha_angle_between_planes_list,lengths_between_planes_list


def calc_3D_alpha(edge_directions, i, normal_at_bendpoint_0_tape, tape_orientation_before_bend_KOS_in_global_KOS, y_direction_tape):
    
    # Calc alpha. Angle between y_tape and direction of edge. Transformation in local tape KOS to get direction of angle
    edge_and_tape_side = np.stack([edge_directions[i], tape_orientation_before_bend_KOS_in_global_KOS[1]])
    edge_and_tape_side_TapeKOS = translate_and_rotate_points_from_OLD_to_NEW_KOS(edge_and_tape_side,tape_orientation_before_bend_KOS_in_global_KOS,np.asarray([0, 0, 0]))

    # Orienting Edge
    if edge_and_tape_side_TapeKOS[0][1] < 0:
        edge_and_tape_side_TapeKOS[0][:] = -edge_and_tape_side_TapeKOS[0][:]
        edge_directions[i] = -edge_directions[i]
    # Orienting Tape
    if edge_and_tape_side_TapeKOS[1][1] < 0:
        edge_and_tape_side_TapeKOS[1][:] = -edge_and_tape_side_TapeKOS[1][:]
    try:
        alpha_angle = math.acos(np.dot(edge_directions[i], tape_orientation_before_bend_KOS_in_global_KOS[1]))
    except:
        alpha_angle = 0
    if  edge_and_tape_side_TapeKOS[0][0] < 0:
        alpha_angle = (math.pi - alpha_angle)
    return alpha_angle
def calc_3D_beta(bend_pts_xyz_trendline, i, normal_at_bendpoint_0_tape, normal_at_bendpoint_1_tape):
    # Calc angle between two consecutiv tape normals: beta
    beta_angle_between_planes = math.acos(np.dot(normal_at_bendpoint_0_tape, normal_at_bendpoint_1_tape))
    # Compare slope of old and new x_direction
    x_direction_before_bend_local = norm_vector(bend_pts_xyz_trendline[i] - bend_pts_xyz_trendline[i - 1])
    x_direction_after_bend_without_rotation_local = norm_vector(
        bend_pts_xyz_trendline[i + 1] - bend_pts_xyz_trendline[i])
    if x_direction_before_bend_local[2] < x_direction_after_bend_without_rotation_local[2]:
        beta_angle_between_planes = -beta_angle_between_planes
    return beta_angle_between_planes
def calc_edge_directions(bend_pts_xyz_global_left, bend_pts_xyz_global_right):
    global counter_failed_matches_of_edges

    while len(bend_pts_xyz_global_right) != len(bend_pts_xyz_global_left):  # Comment_DKu_Wenzel: Just looking for the first Bend would optimize stability and would be faster
        if len(bend_pts_xyz_global_right) < len(bend_pts_xyz_global_left): bend_pts_xyz_global_left = np.delete(bend_pts_xyz_global_left,-1,axis=0)
        else:  bend_pts_xyz_global_right = np.delete(bend_pts_xyz_global_right,-1,axis=0) # right < left
        print("left and right not same amount of bendingpoints")
        counter_failed_matches_of_edges += 1
        if counter_failed_matches_of_edges > 25:
            print("Please restart, matching of left an right edges failed. Maybe try other width.")
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

# Functions in calc local bending points
def calc_delta_length_at_bend(width,alpha):
    if alpha > math.pi / 2:
        delta_length_bend = (width / 2) * math.tan(math.pi - alpha)

    else:  # alpha_list[0] < math.pi / 2:
        delta_length_bend = - (width / 2) * math.tan(alpha)
    return delta_length_bend
def calc_local_KOS_in_trendline_KOS(x_start_end_point_trendline_KOS,y_start_end_point_trendline_KOS):
    # incline/x_slope and y-intercept with Leastsquare, y = x_slope*x + y_intercept
    A = np.vstack([x_start_end_point_trendline_KOS, np.ones(len(x_start_end_point_trendline_KOS))]).T
    x_slope, y_intercept = np.linalg.lstsq(A, y_start_end_point_trendline_KOS, rcond=None)[0]

    x_trendline_new_direction = norm_vector(np.asarray((1, x_slope, 0), dtype=np.float32))
    y_trendline_new_direction = norm_vector(np.asarray((-x_slope, 1, 0), dtype=np.float32))
    z_trendline_new_direction = np.asarray((0, 0, 1), dtype=np.float32)

    local_KOS_in_trendline_KOS = np.vstack((x_trendline_new_direction, y_trendline_new_direction, z_trendline_new_direction))
    return local_KOS_in_trendline_KOS, x_slope
def calc_start_end_point_side_in_trendline_KOS(calc_left_side, delta_length_end_bend, delta_length_start_bend,
                                               end_point_drawn, grid_resolution_int, start_point_drawn, width,
                                               local_KOS_in_trendline_KOS, z_grid_values_linear_trendline_KOS):

    # We need to consider how the Tape is placed 3D

    if calc_left_side:
        Start_point_side_trendline_KOS = start_point_drawn - local_KOS_in_trendline_KOS[1] * width / 2 + delta_length_start_bend * local_KOS_in_trendline_KOS[0]
        End_point_side_trendline_KOS = end_point_drawn - local_KOS_in_trendline_KOS[1] * width / 2 + delta_length_end_bend * local_KOS_in_trendline_KOS[0]

        if new_start_left != []:
            Start_point_side_trendline_KOS = new_start_left

    else:
        Start_point_side_trendline_KOS = start_point_drawn + local_KOS_in_trendline_KOS[1] * width / 2 - delta_length_start_bend * local_KOS_in_trendline_KOS[0]
        End_point_side_trendline_KOS = end_point_drawn + local_KOS_in_trendline_KOS[1] * width / 2 - delta_length_end_bend * local_KOS_in_trendline_KOS[0]

        if new_start_right != []:
            Start_point_side_trendline_KOS = new_start_right



    x_data_side_start_end = (Start_point_side_trendline_KOS[0], End_point_side_trendline_KOS[0])
    y_data_side_start_end = (Start_point_side_trendline_KOS[1], End_point_side_trendline_KOS[1])

    end_point_drawn_side, start_point_drawn_side, x_start_index_side, x_end_index_side = calc_Start_End_in_trendline_KOS_from_xdata_ydata(
        grid_resolution_int, x_data_side_start_end, y_data_side_start_end, z_grid_values_linear_trendline_KOS)
    return end_point_drawn_side, start_point_drawn_side, x_end_index_side, x_start_index_side
def calc_local_and_global_bendpoints(max_distance, surface_points_local_KOS,
                                     local_KOS_in_trendline_KOS, y_intercept_trendline_KOS):

    x_z_points = np.stack([surface_points_local_KOS[:, 0],surface_points_local_KOS[:, 2]],axis=1)
    #y_smooth_left = smooth_savgol(x_z_points_filled_up, poly_order, savgol_window_quotient)
    bend_pts_xz_local = calc_bend_pts(max_distance, x_z_points)
    # in local KOS y=0. Insert this to bendpoints to get (x,y,z)
    bend_pts_xyz_local = np.insert(bend_pts_xz_local, 1, np.zeros(len(bend_pts_xz_local)), axis=1)
    bend_pts_xyz_trendline, bend_pts_xyz_global = new_bending_points_local_to_global(y_intercept_trendline_KOS, bend_pts_xyz_local,
                                                                                     local_KOS_in_trendline_KOS)
    return bend_pts_xz_local, bend_pts_xyz_trendline, bend_pts_xyz_global
def calc_Start_End_in_trendline_KOS_from_xdata_ydata(grid_resolution_int, xdata, ydata, z_grid_values_linear_trendline_KOS):
    # Start und Endpunkt der Eingezeichnet wurde

    y_end_index = np.asarray(np.round(np.add(np.divide(ydata[1], dy), (grid_resolution_int - y_0_grid_point_index))),
                             dtype=np.int32)
    x_end_index = np.asarray(np.round(np.add(np.divide(xdata[1], dx), (grid_resolution_int - x_0_grid_point_index))),
                             dtype=np.int32)
    y_start_index = np.asarray(np.round(np.add(np.divide(ydata[0], dy), (grid_resolution_int - y_0_grid_point_index))),
                               dtype=np.int32)
    x_start_index = np.asarray(np.round(np.add(np.divide(xdata[0], dx), (grid_resolution_int - x_0_grid_point_index))),
                               dtype=np.int32)

    # Left and right starting point can be outside the grit, when they are outside, they get a default value. The z-Data would be needed for the plot.
    if x_start_index < 0 or x_start_index >= grid_resolution_int or \
            y_start_index < 0 or y_start_index >= grid_resolution_int:
        z_start_data = z_grid_values_linear_trendline_KOS[0, 0]
    else: z_start_data = z_grid_values_linear_trendline_KOS[x_start_index, y_start_index]

    if x_end_index < 0 or x_end_index >= grid_resolution_int or \
            y_end_index < 0 or y_end_index >= grid_resolution_int:
        z_end_data = z_grid_values_linear_trendline_KOS[grid_resolution_int - 1, grid_resolution_int - 1]
    else: z_end_data = z_grid_values_linear_trendline_KOS[x_end_index, y_end_index]


    start_point_xyz_data = (np.vstack((xdata[0], ydata[0], z_start_data)).T)[0][:]
    end_point_xyz_data = (np.vstack((xdata[1], ydata[1], z_end_data)).T)[0][:]
    return end_point_xyz_data, start_point_xyz_data, x_start_index, x_end_index
def trim_x_y_values_to_grid(grid_resolution_int, x_slope, x_values_trendline_KOS, x_values_indizes, y_values_trendline_KOS, y_values_indizes):

    y_start_index = -1
    y_end_index = -1
    for k in range(len(x_values_trendline_KOS) - 1):

        if x_slope >= 0:
            if (y_values_indizes[k] >= 0) & (y_start_index < 0): y_start_index = k
            if (y_values_indizes[k] >= grid_resolution_int) & (y_end_index < 0): y_end_index = k - 1

        if x_slope < 0:
            if (y_values_indizes[k] < grid_resolution_int) & (y_start_index < 0): y_start_index = k
            if (y_values_indizes[k] <= 0) & (y_end_index < 0): y_end_index = k - 1
    if y_end_index <= 0: y_end_index = grid_resolution_int - 2 # default value

    # Trim indizes and coordinates to grid size
    x_values_indizes_trim = x_values_indizes[y_start_index:y_end_index]
    y_values_indizes_trim = y_values_indizes[y_start_index:y_end_index]
    x_values_trim = x_values_trendline_KOS[y_start_index:y_end_index]
    y_values_trim = y_values_trendline_KOS[y_start_index:y_end_index]
    return x_values_indizes_trim, x_values_trim, y_values_indizes_trim, y_values_trim
def trim_x_y_values_to_start_end_point(x_start_index, x_end_index, x_values_trendline_KOS, x_values_indizes, y_values_trendline_KOS,
                                       y_values_indizes):
     # Trim indizes and coordinates to grid size
    if x_start_index < x_end_index:
        x_values_indizes_trim = x_values_indizes[x_start_index:x_end_index]
        y_values_indizes_trim = y_values_indizes[x_start_index:x_end_index]
        x_values_trim_trendline_KOS = x_values_trendline_KOS[x_start_index:x_end_index]
        y_values_trim_trendline_KOS = y_values_trendline_KOS[x_start_index:x_end_index]
    else:
        x_values_indizes_trim = x_values_indizes[x_end_index:x_start_index]
        y_values_indizes_trim = y_values_indizes[x_end_index:x_start_index]
        x_values_trim_trendline_KOS = x_values_trendline_KOS[x_end_index:x_start_index]
        y_values_trim_trendline_KOS = y_values_trendline_KOS[x_end_index:x_start_index]

    return x_values_indizes_trim, x_values_trim_trendline_KOS, y_values_indizes_trim, y_values_trim_trendline_KOS
def extract_points_from_interpolated_surface(Start_point_trendline_KOS, grid_resolution_int, x_slope, z_grid_values_linear_trendline_KOS, x_start_index,
                                             x_end_index):
    x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim, y_intercept_trendline_KOS, y_values_trim_trendline_KOS = calc_x_y_values_in_trendline_KOS_and_with_indizes(
        Start_point_trendline_KOS, grid_resolution_int, x_end_index, x_slope, x_start_index)

    # z Values from grid
    z_values_trendline_KOS = extract_z_values_from_grid(x_indizes_trim, y_indizes_trim,
                                                                z_grid_values_linear_trendline_KOS)

    # In z_grid_values_linear_trendline_KOS can be NaN. Trim those of.
    x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim, y_values_trim_trendline_KOS, z_values_trendline_KOS = trim_x_y_z_values_to_geometry(
        x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim, y_values_trim_trendline_KOS,
        z_values_trendline_KOS)
    
    # x, y and z stacked together to 3D-Points
    new_bending_direction_points_trendline_KOS = np.vstack(
        (x_values_trim_trendline_KOS, y_values_trim_trendline_KOS, z_values_trendline_KOS)).T
    return new_bending_direction_points_trendline_KOS, y_intercept_trendline_KOS, x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim, y_values_trim_trendline_KOS
def extract_z_values_from_grid(x_indizes_trim, y_indizes_trim, z_grid_values_linear_trendline_KOS):
    z_values_trendline_KOS = []
    for i in range(len(x_indizes_trim)):
        z_values_trendline_KOS.append(
            z_grid_values_linear_trendline_KOS[x_indizes_trim[i], y_indizes_trim[i]])
    z_values_trendline_KOS = np.asarray(z_values_trendline_KOS, dtype=np.float32)
    return z_values_trendline_KOS
def trim_x_y_z_values_to_geometry(x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim,
                                  y_values_trim_trendline_KOS, z_values_trendline_KOS):
    i = 0
    while z_values_trendline_KOS[i] != z_values_trendline_KOS[i]: i += 1
    j = -1
    while z_values_trendline_KOS[j] != z_values_trendline_KOS[j]:  j -= 1
    z_values_trendline_KOS = z_values_trendline_KOS[i:j]
    x_values_trim_trendline_KOS = x_values_trim_trendline_KOS[i:j]
    y_values_trim_trendline_KOS = y_values_trim_trendline_KOS[i:j]
    x_indizes_trim = x_indizes_trim[i:j]
    y_indizes_trim = y_indizes_trim[i:j]
    return x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim, y_values_trim_trendline_KOS, z_values_trendline_KOS
def calc_x_y_values_in_trendline_KOS_and_with_indizes(Start_point_trendline_KOS, grid_resolution_int, x_end_index,
                                                      x_slope, x_start_index):
    ##### Oblique linear line y(x)
    # 2 lines needed:   • mathematical decription with coordinates for Plot.
    #                   • with Indizies for extraction from grid
    # new y_intercept_trendline_KOS, inserting startpoint into formula:y = x_slope*x + y_intercept_trendline_KOS
    y_intercept_trendline_KOS = Start_point_trendline_KOS[1] - Start_point_trendline_KOS[0] * x_slope
    # x-y-values in coordinates
    x_values_trendline_KOS = grid_x[:, 0]
    y_values_trendline_KOS = np.add(np.multiply(x_values_trendline_KOS, x_slope), y_intercept_trendline_KOS)
    # x-y-values with Indizies
    x_values_indizes = np.asarray(list(range(grid_resolution_int)), dtype=np.int32)
    y_values_indizes = np.add(np.divide(y_values_trendline_KOS, dy), (grid_resolution_int - y_0_grid_point_index))
    y_values_indizes = np.asarray(np.round(y_values_indizes), dtype=np.int32)
    # Trim x and y valules to start/end point and if nec
    x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim, y_values_trim_trendline_KOS = trim_x_y_values_to_start_end_point(
        x_start_index, x_end_index, x_values_trendline_KOS, x_values_indizes, y_values_trendline_KOS, y_values_indizes)
    if min(y_values_indizes) < 0 or max(y_values_indizes) > grid_resolution_int - 1:
        x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim, y_values_trim_trendline_KOS = trim_x_y_values_to_grid(
            grid_resolution_int, x_slope, x_values_trim_trendline_KOS, x_indizes_trim, y_values_trim_trendline_KOS,
            y_indizes_trim)
    return x_indizes_trim, x_values_trim_trendline_KOS, y_indizes_trim, y_intercept_trendline_KOS, y_values_trim_trendline_KOS
def transform_surface_points_to_global_and_local_KOS(y_intercept_trendline_KOS, new_bending_direction_points_trendline_KOS,
                                                     local_KOS_in_trendline_KOS):
    # new zero
    new_zero = np.asarray((0, y_intercept_trendline_KOS, 0), dtype=np.float32)

    local_KOS_in_global_KOS = translate_and_rotate_points_from_OLD_to_NEW_KOS(
        local_KOS_in_trendline_KOS, trendline_KOS_in_global_KOS, np.asarray([0, 0, 0]), True)

    new_bending_direction_points_global_KOS = translate_and_rotate_points_from_OLD_to_NEW_KOS(
        new_bending_direction_points_trendline_KOS, trendline_KOS_in_global_KOS, center_of_pointcloud_weighted_in_global_KOS, True)

    surface_points_local_KOS = translate_and_rotate_points_from_OLD_to_NEW_KOS(
        new_bending_direction_points_trendline_KOS, local_KOS_in_trendline_KOS, new_zero)

    return surface_points_local_KOS, new_bending_direction_points_global_KOS, local_KOS_in_global_KOS
def new_bending_points_local_to_global(y_intercept_trendline_KOS,
                                       surface_points_local_KOS,
                                       local_KOS_in_trendline_KOS):
    # new zero
    new_zero = np.asarray((0, y_intercept_trendline_KOS, 0), dtype=np.float32)

    new_bending_direction_points_trendline_KOS = translate_and_rotate_points_from_OLD_to_NEW_KOS(
        surface_points_local_KOS, local_KOS_in_trendline_KOS, new_zero, True)

    new_bending_direction_points_global_KOS = translate_and_rotate_points_from_OLD_to_NEW_KOS(
        new_bending_direction_points_trendline_KOS, trendline_KOS_in_global_KOS, center_of_pointcloud_weighted_in_global_KOS, True)

    return new_bending_direction_points_trendline_KOS, new_bending_direction_points_global_KOS

#######################################################################################################################
def calc_L_aim(x_y_points_filled_up):
    L_aim = 0
    for i in range(1, len(x_y_points_filled_up)):

        L_aim =L_aim + np.linalg.norm(x_y_points_filled_up[i]-x_y_points_filled_up[i - 1])

    """L_aim = 0
    for i in range(len(surface_points_local_KOS_stacked)):
        for j in range(1, len(surface_points_local_KOS_stacked[i])):
            distance = calc_distance_between_two_points(surface_points_local_KOS_stacked[i][j - 1],
                                                        surface_points_local_KOS_stacked[i][j])

            L_aim = L_aim + distance
    """
    return L_aim
#######################################################################################################################
def show_startstrip(bestPatch_patternpoints,patch_start,patch_end,dimension):
    ###############3D-PLOTTING################
    figure = pyplot.figure() #Comment_DB: 3D plot of objective shape
    axes = mplot3d.Axes3D(figure)
    plt.title('Start solution preprocessor '+ dimension)

    patch_visual = mplot3d.art3d.Poly3DCollection(triangle_vectors_of_stl, linewidths=1, alpha=0.5, edgecolor=[1, 1, 1], label ='Geometry') #Comment_DB: added edgecolor to make the edges visible

    axes.scatter([999999990],[9999999900],[9999999900],linewidths=0.0001, alpha = 0.5, label = "Geometry") #Comment_DB: label in legend
    """
    for corner_point in tri_corner_points_trendline_KOS:
        axes.scatter(corner_point[0], corner_point[1], corner_point[2], c='black')
    
    # Trendline KOS
    # x-Achse
    global trendline_KOS_in_global_KOS
    x1, y1, z1 = [center_of_pointcloud_weighted_in_global_KOS[0],
                  center_of_pointcloud_weighted_in_global_KOS[0] + 200 * trendline_KOS_in_global_KOS[0][0]], \
                 [center_of_pointcloud_weighted_in_global_KOS[1],
                  center_of_pointcloud_weighted_in_global_KOS[1] + 200 * trendline_KOS_in_global_KOS[0][1]], \
                 [center_of_pointcloud_weighted_in_global_KOS[2],
                  center_of_pointcloud_weighted_in_global_KOS[2] + 200 * trendline_KOS_in_global_KOS[0][2]]
    axes.plot(x1,y1,z1, c='red', label ='x-Achse')
    # y-Achse
    x2, y2, z2 = [center_of_pointcloud_weighted_in_global_KOS[0],
                  center_of_pointcloud_weighted_in_global_KOS[0] + 200 * trendline_KOS_in_global_KOS[1][0]], \
                 [center_of_pointcloud_weighted_in_global_KOS[1],
                  center_of_pointcloud_weighted_in_global_KOS[1] + 200 * trendline_KOS_in_global_KOS[1][1]], \
                 [center_of_pointcloud_weighted_in_global_KOS[2],
                  center_of_pointcloud_weighted_in_global_KOS[2] + 200 * trendline_KOS_in_global_KOS[1][2]]
    axes.plot(x2, y2, z2, c='blue', label ='y-Achse')
    # z-Achse
    x3, y3, z3 = [center_of_pointcloud_weighted_in_global_KOS[0],
                  center_of_pointcloud_weighted_in_global_KOS[0] + 200 * trendline_KOS_in_global_KOS[2][0]], \
                 [center_of_pointcloud_weighted_in_global_KOS[1],
                  center_of_pointcloud_weighted_in_global_KOS[1] + 200 * trendline_KOS_in_global_KOS[2][1]], \
                 [center_of_pointcloud_weighted_in_global_KOS[2],
                  center_of_pointcloud_weighted_in_global_KOS[2] + 200 * trendline_KOS_in_global_KOS[2][2]]
    axes.plot(x3, y3, z3, c='green', label ='z-Achse')
    """
    # Midpoint
    """
    axes.scatter(center_of_pointcloud_weighted_in_global_KOS[0],center_of_pointcloud_weighted_in_global_KOS[1],center_of_pointcloud_weighted_in_global_KOS[2],c='g')
    """

    # Start- Endpoint
    axes.scatter(patch_start[0], patch_start[1], patch_start[2], c="green")
    axes.scatter(patch_end[0],patch_end[1],patch_end[2],c='black')
    """

    axes.plot(
        [center_of_pointcloud_weighted_in_global_KOS[0], center_of_pointcloud_weighted_in_global_KOS[0] + 20 * trendline_KOS_in_global_KOS[0][0]],
        [center_of_pointcloud_weighted_in_global_KOS[1], center_of_pointcloud_weighted_in_global_KOS[1] + 20 * trendline_KOS_in_global_KOS[0][1]],
        [center_of_pointcloud_weighted_in_global_KOS[2], center_of_pointcloud_weighted_in_global_KOS[2] + 20 * trendline_KOS_in_global_KOS[0][2]],
        c='red', label='x-Axis_Trendline')
    axes.plot(
        [center_of_pointcloud_weighted_in_global_KOS[0],
         center_of_pointcloud_weighted_in_global_KOS[0] + 20 * trendline_KOS_in_global_KOS[1][0]],
        [center_of_pointcloud_weighted_in_global_KOS[1],
         center_of_pointcloud_weighted_in_global_KOS[1] + 20 * trendline_KOS_in_global_KOS[1][1]],
        [center_of_pointcloud_weighted_in_global_KOS[2],
         center_of_pointcloud_weighted_in_global_KOS[2] + 20 * trendline_KOS_in_global_KOS[1][2]],
        c='green', label='y-Axis_Trendline')
    axes.plot(
        [center_of_pointcloud_weighted_in_global_KOS[0],
         center_of_pointcloud_weighted_in_global_KOS[0] + 20 * trendline_KOS_in_global_KOS[2][0]],
        [center_of_pointcloud_weighted_in_global_KOS[1],
         center_of_pointcloud_weighted_in_global_KOS[1] + 20 * trendline_KOS_in_global_KOS[2][1]],
        [center_of_pointcloud_weighted_in_global_KOS[2],
         center_of_pointcloud_weighted_in_global_KOS[2] + 20 * trendline_KOS_in_global_KOS[2][2]],
        c='blue', label='z-Axis_Trendline')

    axes.plot(
        [0,0 + 20 *1],
        [0, 0 + 20 * 0],
        [0, 0 + 20 * 0],
        c='red', label='x-Axis_Global')
    axes.plot(
        [0, 0 + 20 * 0],
        [0, 0 + 20 * 1],
        [0, 0 + 20 * 0],
        c='green', label='y-Axis_Global')
    axes.plot(
        [0, 0 + 20 * 0],
        [0, 0 + 20 * 0],
        [0, 0 + 20 * 1],
        c='blue', label='z-Axis_Global')

    
    axes.plot(
        [bend_pts_xyz_global_2D[0][0],
         bend_pts_xyz_global_2D[0][0] + 20 * x_direction_list_global_KOS[0][0]],
        [bend_pts_xyz_global_2D[0][1],
         bend_pts_xyz_global_2D[0][1] + 20 * x_direction_list_global_KOS[0][1]],
        [bend_pts_xyz_global_2D[0][2],
         bend_pts_xyz_global_2D[0][2] + 20 * x_direction_list_global_KOS[0][2]],
        c='blue', label='z-Axis')
    
    axes.plot(
        [bend_pts_xyz_global_2D[0][0],
         bend_pts_xyz_global_2D[0][0] + 20 * normal_direction_start[0]],
        [bend_pts_xyz_global_2D[0][1],
         bend_pts_xyz_global_2D[0][1] + 20 * normal_direction_start[1]],
        [bend_pts_xyz_global_2D[0][2],
         bend_pts_xyz_global_2D[0][2] + 20 * normal_direction_start[2]],
        c='blue', label='z-Axis')
    
    """
    """
    for i in range(len(surface_points_global_KOS_stacked)):
        plt.plot((surface_points_global_KOS_stacked[i][:, 0]), (surface_points_global_KOS_stacked[i][:, 1]), (surface_points_global_KOS_stacked[i][:, 2]), marker='o', c='blue')
        plt.plot((surface_points_global_KOS_right_stacked[i][:, 0]), (surface_points_global_KOS_right_stacked[i][:, 1]), (surface_points_global_KOS_right_stacked[i][:, 2]), marker='o', c='blue')
        plt.plot((surface_points_global_KOS_left_stacked[i][:, 0]), (surface_points_global_KOS_left_stacked[i][:, 1]), (surface_points_global_KOS_left_stacked[i][:, 2]), marker='o', c='blue')
     """
    """
    plt.plot((surface_points_start_end_global_KOS[:, 0]), (surface_points_start_end_global_KOS[:, 1]),
             (surface_points_start_end_global_KOS[:, 2]), marker='o', c='blue')
    plt.plot((surface_points_start_end_global_KOS_left[:, 0]), (surface_points_start_end_global_KOS_left[:, 1]),
             (surface_points_start_end_global_KOS_left[:, 2]), marker='o', c='blue')
    plt.plot((surface_points_start_end_global_KOS_right[:, 0]), (surface_points_start_end_global_KOS_right[:, 1]),
             (surface_points_start_end_global_KOS_right[:, 2]), marker='o', c='blue')
    """
    # Bendpoints and edge direction 2DE
    """plt.plot((bend_points_start_end_global_KOS[:, 0]), (bend_points_start_end_global_KOS[:, 1]),
             (bend_points_start_end_global_KOS[:, 2]), marker='o', c='red')
    plt.plot((bend_points_start_end_global_KOS_left[:, 0]), (bend_points_start_end_global_KOS_left[:, 1]),
             (bend_points_start_end_global_KOS_left[:, 2]), marker='o', c='red')
    plt.plot((bend_points_start_end_global_KOS_right[:, 0]), (bend_points_start_end_global_KOS_right[:, 1]),
             (bend_points_start_end_global_KOS_right[:, 2]), marker='o', c='red')

    global edge_directions_global_2DE
    edge_line_global_current_direction = [
        (bend_points_start_end_global_KOS[:] + np.multiply(edge_directions_global_2DE[:], 50)),
        (bend_points_start_end_global_KOS[:] - np.multiply(edge_directions_global_2DE[:], 50))]

    for i in range(len(edge_line_global_current_direction[0])):
        axes.plot([edge_line_global_current_direction[0][i][0], edge_line_global_current_direction[1][i][0]],
                  [edge_line_global_current_direction[0][i][1], edge_line_global_current_direction[1][i][1]],
                  [edge_line_global_current_direction[0][i][2], edge_line_global_current_direction[1][i][2]],
                  c='green')  # Comment_DB: *pc_axes is *args, and .T is np.transpose"""
    #bend_points_start_end_global_KOS, bend_points_start_end_global_KOS_left, bend_points_start_end_global_KOS_right


    #surface_points_global_KOS_left, surface_points_global_KOS_right, bend_pts_xyz_global_left, bend_pts_xyz_global_right

    """
    
    for k in [0]: #range(len(x_direction_list_global_KOS)):
            axes.plot([bend_pts_xyz_global_3D[k+1][0], bend_pts_xyz_global_3D[k+1][0] + 50 * x_direction_list_global_KOS[k][0]],
                      [bend_pts_xyz_global_3D[k+1][1], bend_pts_xyz_global_3D[k+1][1] + 50 * x_direction_list_global_KOS[k][1]],
                      [bend_pts_xyz_global_3D[k+1][2], bend_pts_xyz_global_3D[k+1][2] + 50 * x_direction_list_global_KOS[k][2]], c='red', label ='Tape direction')
    for k in [0]: #range(len(x_direction_rotated_list_global_KOS)):
            axes.plot([bend_pts_xyz_global_3D[k+1][0], bend_pts_xyz_global_3D[k+1][0] + 50 * x_direction_rotated_list_global_KOS[k][0]],
                      [bend_pts_xyz_global_3D[k+1][1], bend_pts_xyz_global_3D[k+1][1] + 50 * x_direction_rotated_list_global_KOS[k][1]],
                      [bend_pts_xyz_global_3D[k+1][2], bend_pts_xyz_global_3D[k+1][2] + 50 * x_direction_rotated_list_global_KOS[k][2]],
                       c='yellow', label ='Tape direction rotated')
    
    for i in [0]: #range(len(edge_line_global)):
        axes.plot([edge_line_global[i][0][0], edge_line_global[i][1][0]],
                [edge_line_global[i][0][1], edge_line_global[i][1][1]],
                [edge_line_global[i][0][2], edge_line_global[i][1][2]],
                 c='green', label ='Edge')  # Comment_DB: *pc_axes is *args, and .T is np.transpose
    
    for k in [1]:  # range(len(x_direction_list_global_KOS)):
        axes.plot([bend_pts_xyz_global_3D[k + 1][0],
                   bend_pts_xyz_global_3D[k + 1][0] + 50 * x_direction_list_global_KOS[k][0]],
                  [bend_pts_xyz_global_3D[k + 1][1],
                   bend_pts_xyz_global_3D[k + 1][1] + 50 * x_direction_list_global_KOS[k][1]],
                  [bend_pts_xyz_global_3D[k + 1][2],
                   bend_pts_xyz_global_3D[k + 1][2] + 50 * x_direction_list_global_KOS[k][2]], c='red')
    for k in [1]:  # range(len(x_direction_rotated_list_global_KOS)):
        axes.plot([bend_pts_xyz_global_3D[k + 1][0],
                   bend_pts_xyz_global_3D[k + 1][0] + 50 * x_direction_rotated_list_global_KOS[k][0]],
                  [bend_pts_xyz_global_3D[k + 1][1],
                   bend_pts_xyz_global_3D[k + 1][1] + 50 * x_direction_rotated_list_global_KOS[k][1]],
                  [bend_pts_xyz_global_3D[k + 1][2],
                   bend_pts_xyz_global_3D[k + 1][2] + 50 * x_direction_rotated_list_global_KOS[k][2]],
                  c='yellow')
    
    for i in range(len(edge_line_global)):
        axes.plot([edge_line_global[i][0][0], edge_line_global[i][1][0]],
                  [edge_line_global[i][0][1], edge_line_global[i][1][1]],
                  [edge_line_global[i][0][2], edge_line_global[i][1][2]],
                  c='green')  # Comment_DB: *pc_axes is *args, and .T is np.transpose
    #"""
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



