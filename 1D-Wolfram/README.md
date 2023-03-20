# Code to implement wolfram cellular automata and generate time space plots corresponding to different rules.

## Description of parameters

### User inputs:

- ***sys_size*** _(int)_: size of the system
- ***init_type*** _(str)_: type of initialization. Two possible values:
  - '***r***'-random (system cells assigned 0 or 1 values randomly)
  - '***c***'-centre (only centre cell assigned 1 while rest all zero)
- ***init_rand_state*** _(int)_: [needed only if init_type is 'r']
- ***BC_type*** _(str)_: Boundary condition to be used. Possible values:
  - '***p***'-periodic boundary condition
  - '***fix-L-R***'-fixed boundary condition whered L and R are fixed states on left and right boundaries e.g. 'fix-1-0' or 'fix-1-1'
- ***rule_number*** _(int)_: Rule to use for evolution of system. Any integer between 0-255
- ***time_steps*** _(int)_: Number of time steps



### Fixed parameters
- ***c_map*** _(str)_: 'summer'; colormap used for plotting system state. Can be changed to any colormap recognized by matplotlib.
- ***plot_type*** _(str)_: 'space-time'; type of plot to be mage. Only one option in current code.
- ***nb_size*** _(int)_: 3; neighborhood size. Can be changed to any odd integer, but computation becomes too expensive for higher values.
- ***n_states*** _(int)_: 2; no. of possible states of a cell. Can be increased, but computation becomes too expensive for higher values.
