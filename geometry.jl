
abstract type Particle end;
abstract type Grid end;

mutable struct Particle2D <: Particle
    x::Vector{Float32};
    v::Vector{Float32};
    C::Matrix{Float32};
    V::Float32;
    J::Float32;
    property::Property2D;
end

mutable struct Grid2D <: Grid
    size_x::Int32;
    size_y::Int32;
    dh::Float32;
    field_v::Matrix{Vector{Float32}};
    field_m::Matrix{Float32};
end

function Particle2D(x0::Vector{Float32}, v0::Vector{Float32}, V0::Float32, p::Property2D)
    C = zeros(2, 2);
    return Particle2D(x0, v0, C, V0, 1.0, p);
end

function getMass(self::Particle2D)
    return self.property.M.rho * self.V;
end

function Grid2D(size_x::Int64, size_y::Int64, dh::Float32)
    filed_v = Matrix{Vector{Float32}}(undef, size_x, size_y);
    filed_m = Matrix{Float32}(undef, size_x, size_y);
    return Grid2D(size_x, size_y, dh, filed_v, filed_m);
end

function initGrid(self::Grid2D)
    for i = 1 : self.size_x
        for j = 1 : self.size_y
            self.field_m[i, j] = 0.0;
            self.field_v[i, j] = [0.0, 0.0];
        end
    end
end