![alt text](https://github.com/asu-trans-ai-lab/DLSim/blob/main/media/DLSIM_Logo_4.jpg)

"Dynamic Traffic Simulation Package with Multi-Resolution Modelling" (DLSim-MRM) is an open source, high-fidelity multi-resolution (i.e., macroscopic, mesoscopic, and microscopic simulation) traffic simulation package which users jointly apply varying temporal and spatial resolutions to solve a single question or set of questions that mirror the physical world with complex intersections. Users can perform traffic assignments and feed results from one model to another while maintaining consistency between the model assumptions. DLSim-MRM typically takes the following steps for simulation based on [General Modeling Network Specification (GMNS)](https://github.com/zephyr-data-specs/GMNS) format, as shown in Figure 1:
1.	Use demand forecasting models to determine overall trip patterns in a regional network, including trip generation, trip distribution, mode split, and initial O-D matrices.
2.	Use mesoscopic simulation-based dynamic traffic assignment (DTA) to realistically assign traffic to the network by accounting for strategic traveler behavior. 
3.	Use microscopic analysis of traffic at the corridor level or subnetwork level.

![alt text](https://github.com/asu-trans-ai-lab/DLSim/blob/main/media/Multiresolution2.jpg)

Figure 1: Different levels of sample road networks for Phoenix, AZ.

DLSim-MRM uses 3 open-source packages; [OSM2GMNS](https://github.com/asu-trans-ai-lab/OSM2GMNS), [Path4GMNS](https://github.com/asu-trans-ai-lab/Path4GMNS) and [Vol2Timing](https://github.com/asu-trans-ai-lab/Vol2Timing) with the additional developments along the multi-resolution modelling and dyanmic traffic simulation.

-OSM2GMNS can help users easily convert networks from OpenStreetMap to .csv files with standard GMNS format for visualization, traffic simulation and planning purpose.

-Path4GMNS is an open-source AMS library for efficiently macroscopic and mesoscopic traffic assignment based on General Modeling Network Specification (GMNS) format.

-Vol2Timing is a python tool aims to offer a light-weight computational engine to generate optimize signal control timing data, and analyze the effectiveness of signal control strategies.

DLSim-MRM has 4 main steps in the simulation process from the beginning as shown in Fig. 2. 

1. After generating GMNS structured node and link files for all simulation levels, first step is reading nodes, links, lanes and agents (auto, walk, bike, bus, truck, cav or ev).

2. In the second step, DLSim-MRM reads the agents demand for desired OD pairs, demand period and cumulative level as presented in input file.

3. Checking OD connectivity and accesibility in the microscopic level, column generation and column-pool based flow updating for traffic assingmnet, and finally OD estimation are the key processes for the traffic assignment in 3rd step.

4. DLSim-MRM simulate the agents behaviour along the multi-resolution network by using queue-VDF, parallel processing, re-routing based memeory management, trajectory generation focusing internal consistency of multi-resolution road network. 

![alt text](https://github.com/asu-trans-ai-lab/DLSim/blob/main/media/DLSIMFlow2.jpg)

Figure 2: Key processes of DLSim-MRM.

## Remarks:

### Developers:

[Dr. Cafer Avci](https://github.com/caferavci)

[Jiawei (Jay) LU](https://github.com/jiawlu)

[Dr. Peiheng Li](https://github.com/jdlph/Path4GMNS)

[Dr. Xuesong (Simon) Zhou](https://github.com/xzhou99)

### Introductory Videos & Slides

[Recent Advances](https://www.youtube.com/watch?v=dj6c6h4mWfI) by Dr. Cafer Avci

Tutorials in [*ASU Transportation AI Lab](https://www.youtube.com/channel/UCpwXRD0kEkR5iQ77iCXCNuQ/videos) by Dr. Xuesong (Simon) Zhou

[Mini-lessson](https://www.youtube.com/watch?v=rorZAhNNOf0&feature=youtu.be) by Dr. Cafer Avci and Dr. Xuesong (Simon) Zhou : What is the best way to learn dynamic traffic simulation and network assignment for a beginner? Do you want to integrate a powerful traffic simulator in your deep learning framework? We would like to offer a collaborative learning experience through 500 lines of python codes and real-life data sets. This is part of our mini-lessons through teaching dialog.

[In-Depth Video Discussions]()

[Slides introducing DLSIM]()

## Theory
[DLSim-MRM-Theory](https://github.com/asu-trans-ai-lab/DLSim/wiki/DLSIM-Theory)

[In-Depth Learning from References](https://github.com/asu-trans-ai-lab/DLSim/wiki/References)

[Data flow chart](https://github.com/asu-trans-ai-lab/DLSim/wiki/data-flow-chart)

[User Guide](https://github.com/asu-trans-ai-lab/DLSim/wiki/User-Guide)

## Tutorials

[DDI Simulation](https://github.com/asu-trans-ai-lab/DLSim/wiki/DDI-tutorial )

[I-10 Corridor Simulation](https://github.com/asu-trans-ai-lab/DLSim/wiki/I10-corridor)

[3 Corridor Simulation](https://github.com/asu-trans-ai-lab/DLSim/wiki/3-corridor)

[Vehicle Routing](https://github.com/asu-trans-ai-lab/DLSim/wiki/vehicle-routing)

## Levels of DLSIM:
Levels of modeling elements:


|Category | Element | GMNS File Names | Learning Goals |
| --- | --- | --- | --- |
| A | Network | node.csv, link.csv | Free-flow speed, capacity, multiresolution network  |
| B | Demand | input_path.csv | Zone structure, OD demand matrix mapping to road network  |
| C | Signal | Timing in link.csv | Micro |
| D | Scenario | setting.csv |  |
| E1 | Link output | link_performance.csv | Macro |
| E2 | Route assignment output | Route_assignment.csv | Macro|
| E3 | Agent output | agent.csv | Agent |
| E4 | Trajectory output | trajectory.csv, trace.csv | Micro |


## A. Network Generation

[Network basics](https://github.com/asu-trans-ai-lab/DLSim/wiki/network-basics)

[Network generation with osm2gmns](https://github.com/asu-trans-ai-lab/DLSim/wiki/network-generation)

[Modification of networks](https://github.com/asu-trans-ai-lab/DLSim/wiki/network-modification)

[Visualization of networks with NeXTA](https://github.com/asu-trans-ai-lab/DLSim/wiki/network-visualization)

## B. Demand Generation

[Trip generation](https://github.com/asu-trans-ai-lab/DLSim/wiki/trip-generation)

[OD demand estimation](https://github.com/asu-trans-ai-lab/DLSim/wiki/OD-demand-estimation)

## C. Traffic Signal Modelling

[Traffic-Signal-Timing-Modeling](https://github.com/asu-trans-ai-lab/DLSim/wiki/Traffic-Signal-Timing-Modeling)

## D. Scenario Management

[Settings of DLSIM](https://github.com/asu-trans-ai-lab/DLSim/wiki/Settings)

## E. Output

[E1. Link Output](https://github.com/asu-trans-ai-lab/DLSim/wiki/link-output)

[E2. Route Assignment Output](https://github.com/asu-trans-ai-lab/DLSim/wiki/route-assignment-output)

[E3. Agent Output](https://github.com/asu-trans-ai-lab/DLSim/wiki/agent-output)

[E4. Trajectory Output](https://github.com/asu-trans-ai-lab/DLSim/wiki/trajectory-output )

## Appendices

[FAQ](https://github.com/asu-trans-ai-lab/DLSim/wiki/Frequently-Asked-Questions)

[ChangeLog](https://github.com/asu-trans-ai-lab/DLSim/wiki/ChangeLog)

[Streamlined-Workflow](https://github.com/asu-trans-ai-lab/DLSim/wiki/Streamlined-Workflow)

## License
[License](https://github.com/asu-trans-ai-lab/DLSim)
