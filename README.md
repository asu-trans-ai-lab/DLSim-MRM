# DLSim

DLSim is an open-source, cross-platform, lightweight, and fast Python traffic assignment tool 
adopted and modified from [ASU Trans+AI Lab](https://github.com/asu-trans-ai-lab/DTALite)

## Installation

DLSim has been published on [PyPI](https://pypi.org/project/dlsim/), and can be installed by using package manager [pip](https://pip.pypa.io/en/stable/) to install DLSim.

```bash
pip install dlsim
```
If you need a specific version of DLSim, say, 0.2.1,
```
$ pip install dlsim==0.2.1
```

If you want to test the latest features of DLSim, you can build the package from sources and install it offline, where **Python 3.x** is required.
```
# from the root directory of DLSim
$ python setup.py sdist bdist_wheel
$ cd dist
$ python -m pip install dlsim-version.tar.gz
``` 
## Usage

Find the shortest path (based on distance) and output it in the format of a sequence of node/link IDs.

```python
import DLSim as ds
import time

network = ds.Network()

ds.g_ReadInputData(network.node_list, 
                   network.link_list, 
                   network.agent_list,
                   network.internal_node_seq_no_dict,
                   network.external_node_id_dict,
                   network.agent_td_list_dict,
                   network.zone_to_nodes_dict,network.node_seq_to_link_seq)
network.allocate()

ds.begin_time = time.time()

ds.g_find_shortest_path_for_agent(network)
ds.g_TrafficSimulation(network.node_list, network.link_list, 
                       network.agent_list, network.agent_td_list_dict)
ds.end_time = time.time()
ds.g_OutputFiles(network.link_list, 
                 network.agent_list, 
                 network.external_node_id_dict)
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
[License](https://github.com/asu-trans-ai-lab/DLSim)