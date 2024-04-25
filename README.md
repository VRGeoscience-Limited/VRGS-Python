# VRGS Python Integration

![VRGS + Python](/assets/images/VRGeoscience_plus_python_v1.png)

This repository contains python code and scripts that can be integrated directly _within_ [Virtual Reality Geological Studio](https://www.vrgeoscience.com/). The aim is to leverage the extensive and expanding ecosystem of earth science specific python packages. The goal is to facilitate advanced 3D outcrop analysis by incorporating python into and complementing current **VRGS** workflows.  

This repository will be periodically updated with code and scripts. 

## Installation
- Install python on your machine --version >= 3.10
    - Recommended to install from the main [python](https://www.python.org/downloads/) build.
    - VRGS will only work with a main python install, it currently does not work with virtual environment installs (e.g. anaconda)

- Connect and point the `Python Path` in the VRGS `Projet Properties` menu to the python directory. This will ensure that the python code is correctly executed. 
    - e.g. `C:\python\python_xxx`

    - Save the VRGS project to implement the directory addition/change.
        - You may also need to restart VRGS.

    - If runtime errors are encountered, please install necessary dependencies.
        - e.g. `pip install xxx`

## Basic Usage
The following steps set up the python script development environment in **VRGS**. 
1. Navigate to the `Collections` menu. 

2. Right + Click on the `Python` option 
    - Click `"New"` to create a new python script. 
    - The built-in python interpreter will open.

3. Import necessary packages and write the awesome python script.

4. :floppy_disk: `Save script` to write python file to project directory. 

5. To run the script, click the :arrow_forward: `Run` command.
    - :exclamation: **Top Tip** - open a `Python Output` console window to view any data output from the `print()` method, or any errors that may occur.
