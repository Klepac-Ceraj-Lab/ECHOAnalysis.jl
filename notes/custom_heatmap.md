From Daniel Schwabeneder via julialang slack:

```julia
using LinearAlgebra, Plots
values = [0, 3, 4, 7, 10]
colors = plot_color.([:white, :yellow, :green, :blue, :black])
grad = cgrad(colors, normalize(values, Inf))
heatmap(rand(10, 10), color = grad)
```

>Actually, you shoul also add `clims = extrema(values)`` in the heatmap call.
