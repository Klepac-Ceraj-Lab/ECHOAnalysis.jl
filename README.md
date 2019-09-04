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

The code in this repository is written in [julia][2].
In order to run it, clone the repo, and activate the environment. If you've just downloaded it, you may need to instantiate it.

[2]: https://julialang.org/

```
$ git clone https://github.com/Klepac-Ceraj-Lab/echo_analysis.git
$ cd echo-analysis
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

(v1.2) pkg> activate . # to get the pkg repl prompt, just press ']'

(echo-analysis) pkg> instantiate
```
