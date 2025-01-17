# distributed-flexible-control

This repository holds the code base used to produce the results for the "Optimal Attitude Control of Large Flexible Space Structures with Distributed Momentum Actuators" paper submitted to the 2025 IEEE Aerospace Conference, available at [arXiv:2410.07376](https://arxiv.org/abs/2410.07376).

## Running the simulations

The main script used to generate the plots and results shown in the paper is the `src/Main.m` with MATLAB R2024a. The code has in lines 16-18 a set of three flags to define the scenario to run.

The `ActPlacementFlagList` defines the actuator placement configurations to run. It holds the name of the two actuator placement configurations that are compared throughout the paper, the "Centralized" and "Distributed" actuator configurations.

The user should set `ScenarioFlag` to "FinePointing" to reproduce the results of the Fine Pointing scenario with reaction wheel jitter perturbation rejection, and "Slewing" to reproduce the results of the slew maneuver with input shaping in coarse pointing mode.

To select the guidance profile, the user should have the `ScenarioFlag` set to "Slewing" plus the `GuidProfileFlag` set to either "TimeOptimal" or "Smoothed". Both will see an input-shaped, bang-coast-bang profile generated for the slewing maneuver. The "TimeOptimal" option will take the reference input profile and generate a dynamically consistent reference trajectory for all states for the controller to track. The "Smoothed" option gives to the controller a reference trajectory considering only rigid-body dynamics with no gyroscopic torque, with the reference for the remaining states all set to zero.

The main script will store the data results in the `Results/SimData` folder, and any generated plots in the `Results/Plots/` folder.

## Citing this work

If you wish to elaborate further on this research or use this code, we would be grateful if you would cite it as:

Cachim, P., Kraus, W., Manchester, Z., Lourenco, P., & Ventura, R. (2024). Optimal Attitude Control of Large Flexible Space Structures with Distributed Momentum Actuators. arXiv preprint arXiv:2410.07376.

```bibtex
@misc{cachim2024optimalattitudecontrollarge,
      title={Optimal Attitude Control of Large Flexible Space Structures with Distributed Momentum Actuators}, 
      author={Pedro Cachim and Will Kraus and Zachary Manchester and Pedro Lourenco and Rodrigo Ventura},
      year={2024},
      eprint={2410.07376},
      archivePrefix={arXiv},
      primaryClass={astro-ph.IM},
      url={https://arxiv.org/abs/2410.07376}, 
}
```
