# Meeting notes

## 2024-09-30

- Fergus did a [Gradus.jl](https://github.com/astro-group-bristol/Gradus.jl) demo today, looking a the [documentation](https://astro-group-bristol.github.io/Gradus.jl/dev/).
- We're using VS Code with the Julia extension.
- You can start the "REPL" in VS Code using SHIFT + Command + P and start typing "Julia: Start REPL" and hitting ENTER.
- You can run commands in the REPL by hitting SHIFT + ENTER.

In the REPL you can open the package manager, find the status of installed packages, add packages, and update packages as follows.
```julia
]
st
add Gradus, Plots
up
```
Then hitting backspace will get you out of the package manager.

- We looked at writing some example Julia code which I've put in the file [test.jl](test.jl).
- You can click on the plot image and step backwards and forwards with the cursor keys.
- If you click on the Julia extension icon on the left hand side of VS Code you will also see a Plot Navigator window on the left which shows the plots and when they were created.
