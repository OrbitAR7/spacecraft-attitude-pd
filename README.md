# Spacecraft Attitude Control (MATLAB)

Attitude-hold simulation of a rigid spacecraft with a 4-wheel RWA and PD control. One full-featured script and one simple, quick-look script.

## Quick start

```matlab
% Full simulation (multi-figure; saves PNGs to results/figures/)
>> attitude_hold_full

% Simple demo (multiple plots, no file export)
>> attitude_hold_simple
```

## What’s inside

* **attitude\_hold\_full.m** — robust shapes, quaternion/Euler dynamics, min-norm wheel allocation, history/metrics, and exported figures.&#x20;
* **attitude\_hold\_simple.m** — straight-ahead demo with quaternion, rates, error, wheel speeds, and torque plots.&#x20;

## Tuning (edit at the top of `attitude_hold_full.m`)

* Inertia: `spacecraft.inertia = diag([6400 4730 8160]);`
* Gains: `controller_params.Kp = 10;  controller_params.Kd = 450;`
* Timebase: `sim_params.dt = 0.1;  sim_params.duration = 1000;`
* Initial state: `spacecraft.quaternion`, `spacecraft.angular_velocity`

> Want more dramatic plots? Start at a 90° attitude, e.g. `[sqrt(2)/2, 0, 0, sqrt(2)/2]`, with zero rates.

## Output

* Figures (full script only): `results/figures/*.png`
* Data (full script only): `results/simulation_*.mat`

## Requirements

* MATLAB (R2018b+) — no toolboxes required.

