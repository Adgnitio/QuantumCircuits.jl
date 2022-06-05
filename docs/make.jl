#Pkg.activate(".")
#Pkg.instantiate() 

using Documenter
using DocumenterTools: Themes

#Pkg.add(PackageSpec(path="../../QuantumCircuits.jl"))
push!(LOAD_PATH,"../src/")
using QuantumCircuits

makedocs(
    modules = [QuantumCircuits],
    sitename = "QuantumCircuits.jl",
    format = Documenter.HTML(mathengine = Documenter.MathJax(), 
                             prettyurls = get(ENV, "CI", nothing) == "true",
                             sidebar_sitename=false
                             ),
    strict = true,
    authors = "Rafał Pracht",
    pages = [
        "Introduction" => "index.md",
        "Quick Start guide" => "quickguide.md",
        "Quantum Gates Library" => [
            "Single-qubit gates s" => "0_Single-qubit-gates.md",
            "Two-qubit gates" => "1_Two-qubit-gates.md",
        ],        
        "Novel algorithm to Simulation on NISQ device" => [
            "Problem definition" => "Simulation/problem_definition.md",
            "U4 - Cartan's KAK decomposition" => "Simulation/U4_Cartan_decomposition.md",
            "Algorithm description" => "Simulation/Trotterization.md",
        ],
        "Examples" => [
            "Log-Normal state preparation" => "examples/state_preparation.md",
        ],        
        "Library References" => [
            "QCircuits" => [
                "Circuit" => "QCircuits/Circuit.md",
                "Gates" => "QCircuits/Gates.md",
                "Qiskit" => "QCircuits/Qiskit.md",
                "Registers" => "QCircuits/Registers.md",
                "Graph" => "QCircuits/Graph.md",
                "OtherGates" => "QCircuits/OtherGates.md",
                "ComplexGates" => "QCircuits/ComplexGates.md",
                "Instructions" => "QCircuits/Instructions.md",
                "Math" => "QCircuits/Math.md",
                "QBase" => "QCircuits/QBase.md",  
            ],
            "Execute" => "Execute.md",
            "QML" => [
                "QML" => "QML/QML.md",
                "Optimization" => "QML/Optimization.md",
                "Gates" => "QML/CircuitOptimization.md",
            ],
        ],
    ],
)
#QuantumCircuits.QCircuits.Gates

deploydocs(
    repo = "github.com/Adgnitio/QuantumCircuits.jl",
    push_preview = true
)
