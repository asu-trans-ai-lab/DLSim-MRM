![alt text](https://github.com/asu-trans-ai-lab/DLSim/blob/main/media/DLSIM-LOGO.jpg)
# DLSim

"Dynamic Traffic Simulation Package" (DLSIM) is an open source, high-fidelity multi-resolution (i.e., macroscopic, mesoscopic, and microscopic simulation) traffic simulation package which users jointly apply varying temporal and spatial resolutions to solve a single question or set of questions that mirror the physical world with complex intersections. Users can perform traffic assignments and feed results from one model to another while maintaining consistency between the model assumptions. DLSIM typically takes the following steps for simulation based on [General Modeling Network Specification (GMNS)](https://github.com/zephyr-data-specs/GMNS) format, as shown in Figure 1:
1.	Use demand forecasting models to determine overall trip patterns in a regional network, including trip generation, trip distribution, mode split, and initial O-D matrices.
2.	Use mesoscopic simulation-based dynamic traffic assignment (DTA) to realistically assign traffic to the network by accounting for strategic traveler behavior. 
3.	Use microscopic analysis of traffic at the corridor level or subnetwork level.

![alt text](https://github.com/asu-trans-ai-lab/DLSim/blob/main/media/Multiresolution2.jpg)

Figure 1: Different levels of sample road networks for Phoenix, AZ.



## Introductory Videos & Slides

[Recent Advances](https://www.youtube.com/watch?v=dj6c6h4mWfI) by Dr. Cafer Avci

Tutorials in [*ASU Transportation AI Lab](https://www.youtube.com/channel/UCpwXRD0kEkR5iQ77iCXCNuQ/videos) by Dr. Xuesong (Simon) Zhou

[[In-Depth Video Discussions]]

[[Slides introducing DLSIM]]

## Theory
[DLSIM-Theory](https://github.com/asu-trans-ai-lab/DLSim/wiki/DLSIM-Theory)

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
| E4 | Trajectory output | trajectory.csv, trace.csv | -Micro |


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


