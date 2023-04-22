"Dynamic Traffic Simulation Package with Multi-Resolution Modelling" (DLSim-MRM) is an open source, high-fidelity multi-resolution (i.e., macroscopic, mesoscopic, and microscopic simulation) traffic simulation package which users jointly apply varying temporal and spatial resolutions to solve a single question or set of questions that mirror the physical world with complex intersections. Users can perform traffic assignments and feed results from one model to another while maintaining consistency between the model assumptions. DLSim-MRM typically takes the following steps for simulation based on [General Modeling Network Specification (GMNS)](https://github.com/zephyr-data-specs/GMNS) format:

1. Use demand forecasting models to determine overall trip patterns in a regional network, including trip generation, trip distribution, mode split, and initial O-D matrices.
2. Use mesoscopic simulation-based dynamic traffic assignment (DTA) to realistically assign traffic to the network by accounting for strategic traveler behavior.
3. Use microscopic analysis of traffic at the corridor level or subnetwork level.

DLSim-MRM uses 3 open-source packages; [OSM2GMNS](https://github.com/asu-trans-ai-lab/OSM2GMNS), [Path4GMNS](https://github.com/asu-trans-ai-lab/Path4GMNS) and [Vol2Timing](https://github.com/asu-trans-ai-lab/Vol2Timing) with the additional developments along the multi-resolution modelling and dyanmic traffic simulation.

-OSM2GMNS can help users easily convert networks from OpenStreetMap to .csv files with standard GMNS format for visualization, traffic simulation and planning purpose.

-Path4GMNS is an open-source AMS library for efficiently macroscopic and mesoscopic traffic assignment based on General Modeling Network Specification (GMNS) format.

-Vol2Timing is a python tool aims to offer a light-weight computational engine to generate optimize signal control timing data, and analyze the effectiveness of signal control strategies.

## Installation:

DLSim has been published on [PyPI](https://pypi.org/project/dlsim/), and can be installed by using package manager [pip](https://pip.pypa.io/en/stable/) to install DLSim.

```
pip install DLSim
```

If you need a specific version of DLSim, say, 0.2.1,

```
$ pip install dlsim==0.2.11
```

## Usage

Find the shortest path (based on distance) and output it in the format of a sequence of node/link IDs.

```
from DLSim import DLSim

# load the DLSim class
DL = DLSim()

# check the working directory
DL.check_working_directory()

# check all the required files exist
DL.check_DLSim_input_files()

# load and update settings
DL.DLSim_settings

# perform kernel network assignment simulation
DL.perform_kernel_network_assignment_simulation()



```

## Contributing

Pull requests are welcome. For major changes, please [open an issue](https://github.com/asu-trans-ai-lab/DLSim-MRM) first to discuss what you would like to change.

## License

[License](https://github.com/asu-trans-ai-lab/DLSim)
