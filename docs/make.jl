using Documenter
using CausalTools 

makedocs(modules=[CausalTools],
         format=Documenter.HTML(),
         sitename="CausalTools Documentation",
          pages=["Home" => "index.md",],)
