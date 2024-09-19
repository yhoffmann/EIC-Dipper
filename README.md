# EIC-Dipper  

Documentation:  
https://www.overleaf.com/read/fxvswvktwmfd
https://www.overleaf.com/read/rtjrybdvgprj#db7d9b

## How to install  

1. Clone this repository:  
```
$ git clone https://github.com/yhoffmann/EIC-Dipper.git
```
2. Enter the directory and clone the submodules (not working right now because I have not made all submodules public):
```
$ cd EIC-Dipper/
EIC-Dipper$ git submodule update --init --recursive
```

## How to build and run
1. Build submpdules and dependencies and the main part of the code.
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