__project__ = "KarstNSim_Public"
__filename__ = "notebook_functions"
__author__ = "Rachel Vasseur"
__date__ = "August 2026"

# Import
import numpy as np
import re
import pandas as pd
import plotly.graph_objects as go
import math

## 3d visualization
# Function to extract triangles and vertices from a file
def mesh_extraction(data):
    """
    Read a file with vertices and triangles.

    Parameters
    ----------
    data: The .txt file already read by pandas

    Returns
    -------
    The vertices and triangles extracted from the .txt file
    """
    vertices = data[data.Type == "VRTX"][['X', 'Y', 'Z']].to_numpy()
    triangles = data[data.Type == "TRGL"][['X', 'Y', 'Z']].to_numpy()

    return vertices, triangles


# Function to plot the points on each segment
def separate_segment(data):
    """
    Connect all the points on each segment.

    Parameters
    ----------
    data: The .txt file already read by pandas

    Returns
    -------
    All the coordinates of the points of the simulation grouped by segment (index)
    """
    x_all, y_all, z_all = [], [], []

    for index, grp in data.groupby("Index"):
        x_all.extend(grp.X.tolist())
        y_all.extend(grp.Y.tolist())
        z_all.extend(grp.Z.tolist())
        # Add None to separate indexes
        x_all.append(None)
        y_all.append(None)
        z_all.append(None)

    return x_all, y_all, z_all


# Function to create a sphere with a set radius
def sphere_mesh(x0, y0, z0, r, resolution=20):
    """
    Transform a point into a sphere mesh.

    Parameters
    ----------
    x0: coordinate X of the sphere
    y0: coordinate Y of the sphere
    z0: coordinate Z of the sphere
    r: radius of the sphere
    resolution: define the resolution of the sphere, default is 20

    Returns
    -------
    Returns a meshed version of the elements of a point
    """
    u = np.linspace(0, 2*np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    u, v = np.meshgrid(u, v)
    x = x0 + r * np.cos(u) * np.sin(v)
    y = y0 + r * np.sin(u) * np.sin(v)
    z = z0 + r * np.cos(v)
    # Conversion in triangles for Mesh3d
    x_flat = x.flatten()
    y_flat = y.flatten()
    z_flat = z.flatten()
    # Index generation of the triangles
    i, j, k = [], [], []
    n = resolution
    for a in range(n-1):
        for b in range(n-1):
            idx = a*n + b
            i.append(idx)
            j.append(idx+1)
            k.append(idx+n)
            i.append(idx+1)
            j.append(idx+n+1)
            k.append(idx+n)

    return x_flat, y_flat, z_flat, i, j, k


## Code to change several parameters in an .txt file
def modify_instruction_file(file, modifs):
    """
    Change a .txt file according to the input modifications.

    Parameters
    ----------
    file: The .txt file
    modifs: a dictionary with the keys to change and their new values

    Returns
    -------
    The .txt file is changed
    """
    with open(file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    new_lines = []

    for line in lines:
        line_stripped = line.strip()

        # Case 1 : commented line such as //key
        match_comment = re.match(r"^\s*//\s*([a-zA-Z0-9_]+)\s*:(.*)$", line)

        if match_comment:
            key = match_comment.group(1)

            if key in modifs:
                # Uncommente + replace value
                new_value = modifs[key]
                new_lines.append(f"{key}: {new_value}\n")
            else:
                # Keep it commented
                new_lines.append(line)
            continue

        # Case 2 : normal line key: value
        match_normal = re.match(r"^\s*([a-zA-Z0-9_]+)\s*:(.*)$", line)

        if match_normal:
            key = match_normal.group(1)

            if key in modifs:
                new_lines.append(f"{key}: {modifs[key]}\n")
            else:
                new_lines.append(line)
        else:
            # Non concerned line (comments, sections, etc.)
            new_lines.append(line)

    with open(file, "w", encoding="utf-8") as f:
        f.writelines(new_lines)


def load_box(filepath):
    """
    Create the box from a given file.

    Parameters
    ----------
    filepath: Path where the .txt file is located

    Returns
    -------
    x_lines, y_lines and z_lines, the data to create the box
    """
    vectors = {}

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            # Ignore empty lines
            if not line:
                continue

            # Stop when table starts
            if line.startswith("Index"):
                break

            parts = line.split()
            key = parts[0]
            values = parts[1:]

            # Only keep 3D vectors
            if len(values) == 3:
                vectors[key] = np.array(list(map(float, values)))

    # Extract basis vectors
    O = vectors["basis"]
    u = vectors["u"]
    v = vectors["v"]
    w = vectors["w"]

    # Build vertices
    vertices = np.array([O, O + u, O + v, O + w, O + u + v, O + u + w, O + v + w, O + u + v + w])

    # Define edges
    edges = [(0,1), (0,2), (0,3), (1,4), (1,5), (2,4), (2,6), (3,5), (3,6), (4,7), (5,7), (6,7)]

    # Build line coordinates
    x_lines, y_lines, z_lines = [], [], []

    for i, j in edges:
        x_lines += [vertices[i, 0], vertices[j, 0], None]
        y_lines += [vertices[i, 1], vertices[j, 1], None]
        z_lines += [vertices[i, 2], vertices[j, 2], None]

    return x_lines, y_lines, z_lines


def add_inputs(fig):
    """
    Import and shows all the inputs on the figure.

    Parameters
    ----------
    fig: Figure to which the inputs will be added

    Returns
    -------
    Add the inputs to the figure
    """
    # Data
    inlets = pd.read_csv('Input_files/1_base/example_sinks.txt', delimiter='\t')  # Inlets
    outlets = pd.read_csv('Input_files/1_base/example_springs.txt', delimiter='\t')  # Outlets
    WT_1 = pd.read_csv('Input_files/1_base/example_watertable_surf1.txt', sep=r'\s+')  # Water table 1
    WT_2 = pd.read_csv('Input_files/1_base/example_watertable_surf2.txt', sep=r'\s+')  # Water table 2
    topo = pd.read_csv('Input_files/1_base/example_topo_surf.txt', sep=r'\s+')  # Topographic surface
    Inc_1 = pd.read_csv('Input_files/1_base/example_inception_surf1.txt', sep=r'\s+')  # Inception surface 1
    Inc_2 = pd.read_csv('Input_files/1_base/example_inception_surf2.txt', sep=r'\s+')  # Inception surface 2
    nkspheres = pd.read_csv('Input_files/2_polygenic_karst/example_nokarstspheres.txt', sep=r'\s+')  # Karst-free spheres

    ## Data extraction
    # Mesh extraction of the surfaces (allows you to display them using Mesh3d afterwards)
    vrtx_WT1, trgl_WT1 = mesh_extraction(WT_1)  # Water table 1
    vrtx_WT2, trgl_WT2 = mesh_extraction(WT_2)  # Water table 2
    vrtx_topo, trgl_topo = mesh_extraction(topo)  # Topographic surface
    vrtx_Inc1, trgl_Inc1 = mesh_extraction(Inc_1)  # Inception surface 1
    vrtx_Inc2, trgl_Inc2 = mesh_extraction(Inc_2)  # Inception surface 2

    # Retrieving data from the boundary box
    X_box, Y_box, Z_box = load_box('Input_files/1_base/example_box.txt')  # Box / limits

    # Inputs
    fig.add_trace(go.Scatter3d(x=X_box, y=Y_box, z=Z_box, mode="lines", line=dict(color="black", width=5), name="Box"))  # Box / limits
    fig.add_trace(go.Scatter3d(x=inlets.X, y=inlets.Y, z=inlets.Z, mode='markers', marker=dict(color="red", size=8, line=dict(color='black', width=2)), name="Inlets"))  # Inlets
    fig.add_trace(go.Scatter3d(x=outlets.X, y=outlets.Y, z=outlets.Z, mode='markers', marker=dict(color="blue", size=8, line=dict(color='black', width=2)), name="Outlets"))  # Outlets
    
    waypoints = pd.read_csv('Input_files/2_polygenic_karst/example_waypoints.txt', delimiter='\t')
    fig.add_trace(go.Scatter3d(x=waypoints.X, y=waypoints.Y, z=waypoints.Z, mode='markers', marker=dict(color="purple", size=8, line=dict(color='black', width=2)), name="Waypoints"))  # Waypoints

    fig.add_trace(go.Mesh3d(x=vrtx_topo[:, 0], y=vrtx_topo[:, 1], z=vrtx_topo[:, 2], i=trgl_topo[:, 0], j=trgl_topo[:, 1], k=trgl_topo[:, 2], intensity=vrtx_topo[:, 2], colorscale='Viridis', opacity=0.3, colorbar=dict(title='Z(m)'), showscale=True, name="Topo", showlegend=True))  # Topographic surface
    fig.add_trace(go.Mesh3d(x=vrtx_WT1[:, 0], y=vrtx_WT1[:, 1], z=vrtx_WT1[:, 2], i=trgl_WT1[:, 0], j=trgl_WT1[:, 1], k=trgl_WT1[:, 2], color='lightblue', name='Water table surface 1', showlegend=True))  # Water table 1
    fig.add_trace(go.Mesh3d(x=vrtx_WT2[:, 0], y=vrtx_WT2[:, 1], z=vrtx_WT2[:, 2], i=trgl_WT2[:, 0], j=trgl_WT2[:, 1], k=trgl_WT2[:, 2], color='lightblue', name='Water table surface 2', showlegend=True))  # Water table 2
    fig.add_trace(go.Mesh3d(x=vrtx_Inc1[:, 0], y=vrtx_Inc1[:, 1], z=vrtx_Inc1[:, 2], i=trgl_Inc1[:, 0], j=trgl_Inc1[:, 1], k=trgl_Inc1[:, 2], color='orange', name='Inception surface 1', showlegend=True))  # Inception surface 1
    fig.add_trace(go.Mesh3d(x=vrtx_Inc2[:, 0], y=vrtx_Inc2[:, 1], z=vrtx_Inc2[:, 2], i=trgl_Inc2[:, 0], j=trgl_Inc2[:, 1], k=trgl_Inc2[:, 2], color='orange', name='Inception surface 2', showlegend=True))  # Inception surface 2

    for idx, row in nkspheres.iterrows():  # Karst-free spheres
        x_mesh, y_mesh, z_mesh, i, j, k = sphere_mesh(row.X.astype(float), row.Y.astype(float), row.Z.astype(float), row.radius.astype(float))
        fig.add_trace(go.Mesh3d(x=x_mesh, y=y_mesh, z=z_mesh, i=i, j=j, k=k, color='lightgreen', opacity=0.5, name=f"Karst-free sphere {int(row.Index.astype(int))}", showlegend=True))
