abstract type BorderCondition end

"""
    Dirichlet(value::Union{Function, Float64})

Structure to define Dirichlet boundary conditions. The value can be a constant or a function of the space variable.
T = g
g = Dirichlet(0.0) # Constant Dirichlet boundary condition
g = Dirichlet(x -> sin(x)) # Dirichlet boundary condition that depends on the space variable
g = Dirichlet((x, t) -> sin(x) * cos(t)) # Dirichlet boundary condition that depends on the space and time variables
"""
struct Dirichlet <: BorderCondition
    value::Union{Function, Float64}
end

"""
    Neumann(value::Union{Function, Float64})

Structure to define Neumann boundary conditions. The value can be a constant or a function of the space variable.
∇T.n = g
g = Neumann(0.0) # Constant Neumann boundary condition
g = Neumann(x -> sin(x)) # Neumann boundary condition that depends on the space variable
g = Neumann((x, t) -> sin(x) * cos(t)) # Neumann boundary condition that depends on the space and time variables
"""
struct Neumann <: BorderCondition
    value::Union{Function, Float64}
end

"""
    Robin(alpha::Union{Function, Float64}, beta::Union{Function, Float64}, value::Union{Function, Float64})

Structure to define Robin boundary conditions. The value can be a constant or a function of the space variable.
αT + β∇T.n = g
g = Robin(0.0, 1.0, 0.0) # Constant Robin boundary condition
g = Robin(x -> sin(x), 1.0, 0.0) # Robin boundary condition that depends on the space variable
g = Robin((x, t) -> sin(x) * cos(t), 1.0, 0.0) # Robin boundary condition that depends on the space and time variables
"""
struct Robin <: BorderCondition
    alpha::Union{Function, Float64}
    beta::Union{Function, Float64}
    value::Union{Function, Float64}
end

"""
    Periodic()

Structure to define periodic boundary conditions.
"""
struct Periodic <: BorderCondition
end

abstract type InterfaceCondition end

"""
    ScalarJump(α₁::Union{Function,Float64}, α₂::Union{Function,Float64}, value::Union{Function,Float64})

Structure to define scalar jump conditions. The value can be a constant or a function of the space variable.
[[αT]] = α₂T2 - α₁T1 = g
g = ScalarJump(0.0, 1.0, 0.0) # Constant scalar jump condition
g = ScalarJump(x -> sin(x), 1.0, 0.0) # Scalar jump condition that depends on the space variable
g = ScalarJump((x, t) -> sin(x) * cos(t), 1.0, 0.0) # Scalar jump condition that depends on the space and time variables
"""
struct ScalarJump <: InterfaceCondition
    α₁::Union{Function,Float64}
    α₂::Union{Function,Float64}
    value::Union{Function,Float64}
end

"""
    FluxJump(β₁::Union{Function,Float64}, β₂::Union{Function,Float64}, value::Union{Function,Float64})

Structure to define flux jump conditions. The value can be a constant or a function of the space variable.
[[β∇T.n]] = β₂∇T2.n - β₁∇T1.n = g
g = FluxJump(0.0, 1.0, 0.0) # Constant flux jump condition
g = FluxJump(x -> sin(x), 1.0, 0.0) # Flux jump condition that depends on the space variable
g = FluxJump((x, t) -> sin(x) * cos(t), 1.0, 0.0) # Flux jump condition that depends on the space and time variables
"""
struct FluxJump <: InterfaceCondition
    β₁::Union{Function,Float64}
    β₂::Union{Function,Float64}
    value::Union{Function,Float64}
end