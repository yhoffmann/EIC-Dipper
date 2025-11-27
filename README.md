# EIC-Dipper (Electron-Ion Collision DIPole-scattering-Process Event generatoR)

_I will probably make the entire calculation a bit more streamlined in the
future._

## Sources
Definitely also check out [this paper](https://arxiv.org/abs/2206.05207)
>_Proton hot spots and exclusive vector meson production_,
Phys. Rev. D **106** (2022) no.7, 074025
[doi:10.1103/PhysRevD.106.074025](https://doi.org/10.1103/PhysRevD.106.074025)

as some of their results are used for verification of the calulations described
in our paper (see following section).

For the rest of the sources, also refer to our paper.

## Paper (to be added once it's published)

## Functionality

A brief explanation of the functionality of the code:
  * Each execution of the binary calculates scattering amplitudes for a single
    event (meaning a single configuration of hotspots)  
     → a seed of `0` (which is the default) results in a random seed; enter a
    seed other than `0` with `--seed` or `-s` for reproducibility, for example
    to remove event-dependence in dense-vs-dilute tests  
    → the main way to change parameters and options (for example number of
    nucleons or size of hotspot) is via command-line arguments
  * There are default values for most settings, but make sure to set the
    `--threads` option, otherwise running will be slow

1. Some parameters are set based on the provided command-line arguments (see
   below)
    * In this step, interpolator data will also be either loaded from file if a
      matching file exists, and if no such file exists, it will be generated and
      written to disk
      * since we do not want every invocation of the executable to rerun the data
        generation, it is recommended to run a single execution with the desired
        parameters and then start a batch job, where every instance can now read
        fully prepared interpolation data
      * I might make precalculated sets of interpolator data available in
      [my CERNBOX](https://cernbox.cern.ch/s/E2nfl1eqASaEcCv) at some point

2. Amplitude results are calculated
    * Calculation might be slow, especially if you forget to specify the number
      of threads to run using the `--threads` (equivalent to `-t`) option
    * You can enable logging to stdout of calculated results with the `-p` or
    `-pp` options
      * `-p`: log only finished results (you usually want this unless testing
        convergence)  
        → one amplitude result for one set of parameters (Delta, phi, etc.)  
        → the output can be understood as follows:  
         `<co/inco> <Delta> Converged after <num of sub-domains in b/bbar> <val> <err> rel err: <rel_err> after <num of integration points total>`
      * `-pp`: also log intermediate results from the integration in sub-domains  
         → intermediate results, meaning results from when an integration over
        a sub-domain finishes  
         → the output for each sub-domain integration is  
          `<co/inco> <Q> <Delta> <phi_Delta> <sub-domain count>
          <"root" if Bessel root else "">(<bmin>,<bmax>)
          <sum>(<partial result>) Error: <partial result err> after
          <num evals for partial><warning if precision not reached>`

3. Amplitude results are written to file

    * If no output file is specified with the `-o` option, the binary will take
      care of automatically placing the file in the correct directory  
      → the directory structure is created when running most `make` targets; it
      is specifically run before executing the binary when using the job
      submission scripts

4. Additional steps that are necessary for cross-section calculation
    * run a few events (128 could be a good start)
    * open `analysis/dsigma_dt.ipynb` and in `enter the folder where data for the
    process was saved (for example `data/samples/c10/de`, meaning sampled
    results, not analytical, charm with a factor of `1.0` ("10") on g²µ_0² and
    dense model ("de") using a proton target)  
    → due to too many parameters and not getting this far in the analysis, any
    nucleus runs don't specify quarkonium type or the factor

## How to install

1. Clone this repository:

```
git clone https://github.com/yhoffmann/EIC-Dipper.git
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
./eic --help
```

to see the possible flags. Typical usage could look like this:

```
EIC-Dipper$ ./eic -A 1 -H 3 -rH2 0.7 -b --g2mu02-factor 0.5 -p --threads $(nproc)
```

This would run one event with 3 hotspots (in one nucleon/proton), where the hotspot radius squared is 0.7, bottom quarks in the vector meson, half the color charge in every hotspot, with progress monitoring to stdout, and your hardware's logical processors fully utilized.

For now, do not set `--g2mu02-factor` to anything but `0.5`, `1.0`, or `2.0`.

## Job submission on high-throughput and high-performance computing clusters WIP

- The directory `job-submission-scripts` contains a few mock scripts to easily be configured for use on `HTCondor` or `slurm`. _NOTE: The scripts will intentionally not work as provided because the necessary setup might be slightly different on different systems._
  - I've never used HTCondor for this specifically, so the code is not really set up to be comfortable to run on it
    - You will need to manually adjust the dir to read from in the analysis scripts  
      → basically just set it to `data/samlpes/` and then set the parameters manually with `set_parameters_manually(...)`  
      → only files that have matching parameters (so presumably only files from the same run) are read and analyzed
    - If you're on a shared filesystem (between login and compute nodes), the above steps are probably not necessary either
  - The `slurm` submit script assumes that you are on a shared filesystem because it is basically always the case  
    → will require some manual setup if your system does not use shared filesystems
- Make sure to submit the scripts while the project root is cwd
  - Example: `condor_submit job-submission-scripts/condor.submit`
