import vrgs
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import pandas as pd

def classify(dip, coplanarity):
    return 1 if coplanarity > 4 else 0
 ## return 1 if dip > 80 and coplanarity > 4 else 0

def best_clusters(df, maximum_K):
    return [KMeans(n_clusters=k, n_init=10).fit(df).inertia_ for k in range(1, maximum_K)], list(range(1, maximum_K))

def elbow_plot(centres, k_values):
    plt.figure(figsize=(12, 6))
    plt.plot(k_values, centres, 'o-', color='orange')
    plt.xlabel("Number of Clusters (K)")
    plt.ylabel("Cluster Inertia")
    plt.title("Elbow Plot of KMeans")
    plt.show()

meshes = vrgs.mesh_list()
if meshes:
    print("List of meshes", meshes)
    mesh_name = meshes[0][0]
    attributes = vrgs.mesh_attribute_list(mesh_name)
    vertices = vrgs.mesh_vertices(mesh_name) 
    x,y,z = zip(*vertices)
    if attributes:
        print(mesh_name, attributes)
        attribute_names = ["Dip", "Azimuth", "CoPlanarity", "NormalX", "NormalY", "NormalZ", "ChangeOfCurvature"]
        attribute_values = [vrgs.mesh_attribute(mesh_name, attr) for attr in attribute_names]
        dip, azimuth, coplanarity, normalx, normaly, normalz, deltacurve = attribute_values
        
        classified_result = list(map(classify, dip, coplanarity))
        if classified_result:
            vrgs.mesh_add_attribute(mesh_name, "classified", classified_result)
        
        frame = pd.DataFrame({
            "dip": dip,
            "x": x,
            "y": y,
            "z": z,
            "coplanarity": classified_result,
            "curve_change": deltacurve
        })
        print(frame)

        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(frame)
##        centres, k_values = best_clusters(scaled_data, 30)
##        elbow_plot(centres, k_values)
        print(scaled_data)
        kmeans_model = KMeans(n_clusters=30, n_init=10).fit(scaled_data)
        frame["clusters"] = kmeans_model.labels_

        plt.scatter(azimuth, frame["dip"], c=frame["clusters"], s=coplanarity)
        plt.show()

        vrgs.mesh_add_attribute(mesh_name, "classified", list(frame["clusters"]))