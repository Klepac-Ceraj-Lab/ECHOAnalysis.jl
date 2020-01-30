# Echo RESONANCE Microbiome Group

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://klepac-ceraj-lab.github.io/echo_analysis/dev/)
[![Build Status](https://travis-ci.org/Klepac-Ceraj-Lab/echo_analysis.svg?branch=master)](https://travis-ci.org/Klepac-Ceraj-Lab/echo_analysis)

[![Kevin Bonham, PhD](https://img.shields.io/badge/Author-Kevin%20Bonham%2C%20PhD-blueviolet)](http://nequals.me)
[![Sophie Rowland](https://img.shields.io/badge/Author-Sophie%20Rowland-blueviolet)](http://sophierowland.com/)
[![Vanja Klepac-Ceraj](https://img.shields.io/badge/Author-Vanja%20Klepec--Ceraj%2C%20PhD-blueviolet)](https://www.vkclab.com/)


## Purpose

Computational analysis of fecal samples and cognitive data from the Resonance cohort by the Klepac-Ceraj lab. This work is funded by the [National Institutes of Health (NIH)][1]

[1]: https://www.nih.gov/echo

## Code

The code in this repository is written in [julia][2]
and supports the analysis of Microbiome data.
The easiest way to use it is to install it in your active project
using julia's `Pkg` package manager.

This package depends on packages in the [BioJulia][3] registry,
which you can also add using the `Pkg` REPL
(this only needs to be done once).

[2]: https://julialang.org/
[3]: https://biojulia.net/getting-started/

```
$ cd ~/my_project
$ julia
Starting Julia...
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.2.0 (2019-08-20)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> ] # this activates the pkg REPL prompt

(v1.3) pkg> registry add https://github.com/BioJulia/BioJuliaRegistry.git

(v1.3) pkg> activate .
Activating new environment at `~/my_project/Project.toml`

(my_project) pkg> add https://github.com/Klepac-Ceraj-Lab/echo_analysis
Updating registry at `~/.julia/registries/BioJuliaRegistry`
  Updating git-repo `https://github.com/BioJulia/BioJuliaRegistry`
  Updating registry at `~/.julia/registries/General`
  Updating git-repo `https://github.com/JuliaRegistries/General.git`
   Cloning git-repo `https://github.com/Klepac-Ceraj-Lab/echo_analysis`
  Updating git-repo `https://github.com/Klepac-Ceraj-Lab/echo_analysis`
  Updating git-repo `https://github.com/Klepac-Ceraj-Lab/echo_analysis`
 Resolving package versions...
 # lots more output ...
 # ...
(my_project) pkg> # press backspace to get back to julia prompt
```
