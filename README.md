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
* An executable should have been generated in 'build/release/karstnsim.exe'.

**Running the code:**
- **Double-click** `karstnsim.exe`: the program scans `Input_files/` for subfolders and shows a numbered menu of available examples. Type the number of the example you want and press Enter.
- **Command line**
  ```bash
  cd path/to/your/executable
  karstnsim.exe ../../../Input_files/[path/to/your/instruction_file.txt]
  ```

*Make sure to use "/" or "\\\\" but never "\\" for the paths.*

Outputs are stored in the outputs directory.

## Documentation

### Input configuration documentation

If you're not looking to code in KarstNSim and only want to use it for simulation applications, you only need information on how to configure a simulation with the (many) inputs of KarstNSim. 
We provide a file in the archive ([here](config_reference.md)) which precises the type, theoretical and typical range of values, meaning, behavior and practical effect of all input parameters.

### Generate documentation files for full code documentation

If you're aiming at coding in KarstNSim or better understand how the code works, a doxyfile is present in the archive. To automatically generate the documentation, 
type `doxygen path/to/YourDoxyfile` in a command prompt, or simply `doxygen doxyfile` if already
 in the root folder. Once generated, you will find the documentation in the *hmtl folder*. It is advised to open the documentation starting from the main page, which provides important
general information about the project structure. You can find it by opening the "index.html" file, or by opening any other .html file and clicking on "Main Page" in the upper left.
A complete documentation of all user input parameters is available in the KarstNSim::ParamsSource struct page.

## Testing

We provide **four ready-to-run examples** (see the `Input_files/` folder). Each example lives in its own subfolder under `Input_files/`. When you double-click `karstnsim.exe`,
the interactive picker lists these folders and preselects the appropriate instruction file for each. The illustrations provided below have been generated with a slightly different version of KarstNSim,
meaning the seed will provide different random numbers and networks will be slightly different by using KarstNSim_Public.

> Tip: Use a 3D viewer such as [ParaView](https://www.paraview.org/download/) to inspect inputs/outputs.

### 1) Base simulation
A minimal, single-phase run on the synthetic dataset (as in Fig. 12 of the 2024 paper, but with some improvements). It demonstrates:
- vadose/phreatic partitioning,
- fracture guidance (two families),
- intrinsic karstification potential,
- inception surfaces.

<img src="" alt="base_example.png" width="100%" align="center">

### 2) Polygenic karst (multi-phase reuse)
Reuses a previously generated network by **down-weighting already traversed edges**, mimicking multi-phase karstification (here, base-level drop with flow redirected to another outlet). The example disables the fracture term to emulate a more homogeneous medium in the bottom of the aquifer, and showcases **waypoints** and **karst-free points** close to the aquifer substratum.

<img src="" alt="polygenic_example.png" width="100%" align="center">

### 3) Amplification
Starting from the **base** network, this example runs the amplification step to control density and topology by adding:
- small-scale dead-end branches,  
- a limited number of cycles/loops.

<img src="" alt="amplification_example.png" width="100%" align="center">

### 4) Karst section generation (curvilinear SGS)
Given the **amplified** network, this example simulates **equivalent radii on nodes** using a **curvilinear SGS** algorithm (Frantz et al., 2021). Variogram parameters mirror those in the paper; radii are visualized in log-scale and mapped to node sphere size for rendering.

<img src="" alt="karst_section_generation_example.png" width="100%" align="center">

The main option of KarstNSim not presented in the above examples is the ghost-rock-driven simulation. Read the config reference file for information on correct use.


---


<b>For any information </b>, please contact : 
* Augustin Gouy : a.gouy.proaddress@gmail.com
* Pauline Collon : pauline.collon@univ-lorraine.fr 
* Christophe Antoine : christophe.antoine@univ-lorraine.fr

<b> To report a bug in the code </b>, please contact :
* Augustin Gouy : a.gouy.proaddress@gmail.com
Provide a screenshot of the log as well as complete information on your input parameters.