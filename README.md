# KarstNSim_Public
Public version of KarstNSim, a C++ code for graph-based and geologically-driven simulation of 3D karst networks.

* [2024 Publication](https://doi.org/10.1016/j.jhydrol.2024.130878)
* [2025 Thesis](https://hal.univ-lorraine.fr/tel-05114757v1)

Its inputs and outputs are ASCII files and it can be run through a single command.
It adapts the Karst simulation code proposed by <b> Paris, A., Guérin, E., Peytavie, A., Collon, P., Galin, E., 2021. Synthesizing Geologically Coherent Cave Networks. Comput. Graph. Forum 40, 277–287. https://doi.org/10.1111/cgf.14420 which is available on Github at : https://github.com/aparis69/Karst-Synthesis. </b>
This implementation includes modifications as compared to this initial independant version, in order to better suit geological data and information. 

The first version of KarstNSim was done in the frame of <b> Benoit Thebault </b> master's thesis, supervised by Pauline Collon. It was presented in the 2022 RINGMeeting in: <b> Thebault, B., Collon, P., Antoine, C., Paris, A., Galin, E., 2022. Karstic network simulation with γ -graphs, in: 2022 RING Meeting. </b>
Since 2022 summer, KarstNSim is developed in the frame of <b> Augustin Gouy</b> PhD thesis, supervised by Pauline Collon and Vincent Bailly-Comte.
This public version corresponds to the version of KarstNSim used to generate results in the 2025 PhD thesis.

It is recommended to read the methodology presented in the 2024 article and/or in the thesis to better apprehend the code.

If you use this code, please cite the associated article :

```
@article{Gouy2024,
author = {Gouy, Augustin and Collon, Pauline and Bailly-Comte, Vincent and Galin, Eric and Antoine, Christophe and Thebault, Beno{\^{i}}t and Landrein, Philippe},
doi = {10.1016/j.jhydrol.2024.130878},
issn = {00221694},
journal = {Journal of Hydrology},
month = {feb},
title = {{KarstNSim: A graph-based method for 3D geologically-driven simulation of karst networks}},
year = {2024}
}
```

## Changelog

Please find the link to the changelog [here](Changelog.md).

## Requirements

* [CMake](https://cmake.org/download/) 3.8 to 3.28 (select the Windows x64 Installer). During installation, check the option "add CMake to the system PATH".
* [Visual Studio 2017](https://my.visualstudio.com/Downloads?q=visual%20studio%202017&wt.mc_id=o~msft~vscom~older-downloads) or newer (for Windows).
* C++14 or newer (can be installed from Visual Studio).
* (Optional) [Doxygen 1.9.6](https://www.doxygen.nl/download.html) or newer for generating documentation.

## Compatibility

KarstNSim is designed to operate on Windows 10. While it hasn't been directly tested on Linux, there is an indication of compatibility based on a successful CMake test build.

## Installation

* Download the archive and unzip it somewhere (avoid spaces and special characters in the path).
* Go to the KarstNSim folder and run the batch file "build.bat", which will create a build folder and run CMake to generate build files and build the project (including compilation).
* An executable should have been generated in build/release/karstnsim.exe. To run the code, you can either double-click directly on karstnsim.exe if the instructions file is named exactly "instructions.txt", or else open a command prompt anywhere and type :

```
cd path/to/your/executable
karstnsim.exe ../../../Input_files/[name_of_instruction_file]
```

The instruction file is in the Input_files directory in the root of the archive. *Make sure to use "/" or "\\\\" but never "\\" for the paths.*

Outputs are stored in the outputs directory.

## Generate documentation files

A doxyfile is present in the archive. To automatically generate the documentation, type `doxygen path/to/YourDoxyfile` in a command prompt, or simply `doxygen doxyfile` if already
 in the root folder. Once generated, you will find the documentation in the *hmtl folder*. It is advised to open the documentation starting from the main page, which provides important
general information about the project structure. You can find it by opening the "index.html" file, or by opening any other .html file and clicking on "Main Page" in the upper left.
A complete documentation of all user input parameters is available in the KarstNSim::ParamsSource struct page.

## Testing

We provide here a synthetic dataset (inspired from the one used in the 2024 article) and an instruction file to simulate results shown in a figure of the 2024 article. <!-- three instruction files (instruction_file_step1, instruction_file_step2 and instruction_file_step3) in the correct format used for KarstNSim, which correspond to inputs used to generate examples in three steps. They are in the Input_files folder.
Since the name of these three files is not "instruction_file.txt", you have to launch the executable manually from a command prompt (see ##Installation for more information).

The first step consists in showing you a simulation, using five inlets going to the upper spring (S1). Subcosts involved are: vadose/phreatic partition, fractures (chosen at orientations N000 and N060), intrinsic karstification potential (including a ghost-rock corridor following the syncline axis), and inception surfaces.
Also, an amplification of the network is made to add some dead-end branches to increase network density, and some cycles to change network topology.

See below pictures of a) model settings (inlets in red, outlets in blue), b) the background grid Poisson-sphere radius property *r*, c) the background grid intrinsic karstification potential property *P*, and d) an example result.

XXX -->
<img src="vadose_contexts_example.png" alt="Figure 12 (Gouy et al., 2024)" width="100%" align="center">

<!-- The second step consists in presenting a polygenic karst simulation: the network generated during step 1 is reused by reducing the cost of edges already traversed. This time, the spring considered is the lower one (S2): this mimicks a base level drop. The ghost-rock volume is still present, influencing the chosen paths.
The cost function is mostly the same, but the fracture component was removed (mimicking a shift from an epikarstic fractured zone to more homogeneous rocks).
Moreover, two waypoints and two karst-free points were added on the substratum of the aquifer, to show their capability to control path position. The amplification step is also applied, with the same parameters as in the first step.

Below you can see a) Settings with waypoints in purple and karst-free points in green, and b) an example result.

XXX

[] The third and last step consists in generating equivalent conduit dimensions on the edges of the full network generated. The variogram parameters used are the ones proposed by Frantz et al. (2021), based on cave survey data. This step involves using a parameter option that allows to add the previously generated networks as
inputs, skipping all steps except the dimensions simulation one. Below is an example simulation result with the equivalent radii represented in log-scale, with values proportional to the sphere radii at each skeleton node.

XXX -->

If you want to make modifications to the input parameters, create a new file instructions.txt, and change options as you want with inspiration from the instruction files provided.

It is advised to use a 3D viewer software (e.g., [ParaView](https://www.paraview.org/download/) ) to visualize complex inputs and outputs (the viewer is not provided in the archive).

<b>For any information </b>, please contact : 
* Augustin Gouy : augustin.gouy@univ-lorraine.fr
* Pauline Collon : pauline.collon@univ-lorraine.fr 
* Christophe Antoine : christophe.antoine@univ-lorraine.fr