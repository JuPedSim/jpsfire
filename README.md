**JPSfire**

Framework for the coupling between [JuPedSim](http://jupedsim.org) and fire simulation results obtained from [FDS](https://github.com/firemodels/fds-smv)

Have a look at the [Wiki](https://gitlab.version.fz-juelich.de/jupedsim/jpsfire/wikis/home)


Usage:
-----

- Create a directory in `demos` (e.g. `Test`)
- In `Test` create two directories `FDS` and `JuPedSim`
  - In `FDS`: put fds-file and all corresponding `fds`-simulation results
  - In `JuPedSim`: put the files necessary to run a `jpscore` simulation, especially an ini-file and a geometry file.
    Have a look into [0_demo/JuPedSim/jpsfire_ini.xml](0_demo/JuPedSim/jpsfire_ini.xml) to see
    how to define the `jpsfire`-relevant information.
  - Run [src/preprocess_fds.py](src/preprocess_fds.py) to generate out of the FDS-simulation the relevant `npz-files`,
    which will be read during `jpscore` simulations.
    You may want to see the available option by running

    ```python
    python preprocess_fds.py -h
    ```
    
    Assuming your files live in a directory called `test` a call could be:
    ```python
    python preprocess_fds.py -j ../demos/test/JuPedSim -f ../demos/test/FDS
    ```
  - Run `jpscore`. Note that in the ini-file, the right paths for `fds` data should be correct. 

Available components:
--------------------

- Smoke sensor

- Walking speed reduction in smoke

- Fire hazard analyses


Upcoming components:
-------------------
Detection


Roadmap
-----

- Further examples for SmokeSensor (3D?)
- Individual detection of smoke
- Use floor field for everything


