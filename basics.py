# VRGS Python methods
import vrgs

# Return the number list of meshes in the project
vrgs.mesh_list()

# Return the list of point clouds in the project
vrgs.pointcloud_list()

# Return the list of vetices in a mesh_list
vrgs.meshvertices()

# Return the list of vertices in a point cloud
vrgs.pclvertices()

# Return the number of attributes of a mesh
vrgs.meshattributelist()

# Return the number of the attributes of a point cloud 
vrgs.pclattributelist()

# Return the attributes of a mesh
vrgs.meshattribute()

# Return the attributes of a point cloud 
vrgs.pclattribute()

# Add attributes to a mesh after computation
vrgs.addmeshattribute()

# Return a list of all orientations in the project
vrgs.orientations()
