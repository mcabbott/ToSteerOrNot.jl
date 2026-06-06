# ToSteerOrNot.jl

How should Sebastian find the island? 

![](docs/sebastian.png)

This repository has code related to "Discrete turn strategies emerge in information-limited navigation" by Jose M. Betancourt, Matthew P. Leighton, Thierry Emonet, Benjamin B. Machta, Michael C. Abbott, [arxiv:2602.23324](https://arxiv.org/abs/2602.23324).

It's organised as a Julia package. The folder `docs` contains the code used to generate the paper's figures. There's also a Pluto notebook which re-creates simpler versions of the main figures, and a Jupyter notebook with the same code. With luck you can run these in the cloud, [live on Binder](https://pluto.land/n/ibrtftuh) and [live Google Colab](https://colab.research.google.com/drive/1J4lp4oJiF0LtSYwnwdjh3IwmTYLOFB_a?usp=sharing) respectively.

If you have [Julia installed locally](https://julialang.org/downloads/) (version 1.11 or later), the following steps will install everything:

```
using Pkg; Pkg.add(url="https://github.com/mcabbott/ToSteerOrNot.jl")
using ToSteerOrNot, Plots, Pluto, DelimitedFiles
```
