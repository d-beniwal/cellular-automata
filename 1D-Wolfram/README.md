# Code to implement wolfram cellular automata and generate time space plots corresponding to different rules.

**Description of parameters**

**User inputs**:

- _sys_size (int)_: size of the system
- _init_type (str)_: type of initialization. Two possible values:
  - 'r'-random (system cells assigned 0 or 1 values randomly)
  - 'c'-centre (only centre cell assigned 1 while rest all zero)
- _init_rand_state (int)_: [needed only if init_type is 'r']
