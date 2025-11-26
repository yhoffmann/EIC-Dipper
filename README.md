# EIC-Dipper  

## Paper (to be added once it's published)

## How to install  

1. Clone this repository:  
```
$ git clone https://github.com/yhoffmann/EIC-Dipper.git
```
2. Enter the directory and clone/update the submodules:
```
$ cd EIC-Dipper/
EIC-Dipper$ git submodule update --init --recursive --remote
```
_NOTE_: If you are ever stuck with old versions of submodules, for example if the above command yields the error `fatal: Unable to find refs/remotes/origin/branch-name revision in submodule path 'path/to/submodule'`, usually just running
```
git submodule sync
```
should fix the problem. Make sure to run the `update` command again afterwards.

## How to build and run
1. Build submodules and dependencies and the main part of the code.
```
EIC-Dipper$ make all DILUTE=X
```
Chose `X = 0` for dense model and `X = 1` for dilute model. Everything else can be configured from the command line. Just run
```
$ ./eic --help
```
to see the possible flags. Typical usage could look like this:
```
EIC-Dipper$ ./eic -A 1 -H 3 -rH2 0.7 -b --g2mu02-factor 0.5 -p --threads $(nproc)
```
This would run one event with 3 hotspots (in one nucleon/proton), where the hotspot radius squared is 0.7, bottom quarks in the vector meson, half the color charge in every hotspot, with progress monitoring to stdout, and your hardware's logical processors fully utilized.

For now, do not set `--g2mu02-factor` to anything but `0.5`, `1.0`, or `2.0`.


## Job submission on high-throughput and high-performance computing clusters
* The directory `hpc-job-submission-scripts` contains a few mock scripts to easily be configured for use on `HTCondor` or `slurm`. _NOTE: The scripts will intentionally not work as provided because the necessary setup might be slightly different on different systems._
