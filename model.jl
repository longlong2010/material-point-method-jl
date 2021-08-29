using LinearAlgebra;

mutable struct Model
    particles::Array{Particle2D};
    grid::Grid;
    gravity::Float32;
end

function Model(grid::Grid2D, gravity::Float32)
    return Model(Array{Particle2D}([]), grid, gravity);
end

function addParticle(self::Model, p::Particle2D)
    if !(p in self.particles)
        push!(self.particles, p);
    end
end

function particle2Grid(self::Model, dt::Float32)
    local dh::Float32 = self.grid.dh;
    local N::Int32 = 3;
    local grid = self.grid;
    for particle in self.particles
        local X = particle.x / dh;
        local base = floor.(Int32, X .- 0.5);
        local x = X - base;

        local w = [0.5 * (1.5 .- x).^2, 0.75 .- (x .- 1).^2, 0.5 * (x .- 0.5).^2];
        local mass = getMass(particle);
        local stress = -dt * 4 * particle.property.M.E * particle.V * (particle.J - 1) / dh^2;
        local affine = [stress 0; 0 stress] + mass * particle.C;
        for i = 1 : N
            for j = 1 : N
                local offset = Vector{Int32}([i - 1, j - 1]);
                local dr = (offset - x) * dh;
                local weight = w[i][1] * w[j][2];
                grid.field_m[base[1] + i, base[2] + j] += weight * mass;
                grid.field_v[base[1] + i, base[2] + j] += weight * (mass * particle.v + affine * dr);
            end
        end
    end
end

function boundaryCondition(self::Model, dt::Float32)
    local bound = 3;
    local grid = self.grid;
    local field_v = grid.field_v;
    local field_m = grid.field_m;
    for i = 1 : grid.size_x
        for j = 1 : grid.size_y
            if field_m[i, j] > 0
                field_v[i, j] /= field_m[i, j];
            end
            field_v[i, j][2] += self.gravity * dt;
            if i < bound && field_v[i, j][1] < 0
                field_v[i, j][1] = 0;
            end
            if i > grid.size_x - bound && field_v[i, j][1] > 0
                field_v[i, j][1] = 0;
            end
            if j < bound && field_v[i, j][2] < 0
                field_v[i, j][2] = 0;
            end
            if j > grid.size_y - bound && field_v[i, j][2] > 0
                field_v[i, j][2] = 0;
            end
        end
    end
end

function grid2Particle(self::Model, dt::Float32)
    local dh::Float32 = self.grid.dh;
    local N::Int32 = 3;
    local grid = self.grid;
    local field_v = grid.field_v;
    for particle in self.particles
        local X = particle.x / dh;
        local base = floor.(Int32, X .- 0.5);
        local x = X - base;
        local w = [0.5 * (1.5 .- x).^2, 0.75 .- (x .- 1).^2, 0.5 * (x .- 0.5).^2];
        local v1 = zeros(2);
        local C1 = zeros(2, 2);
        for i = 1 : N
            for j = 1 : N
                local offset = Vector{Int32}([i - 1, j - 1]);
                local dr = (offset - x) * dh;
                local weight = w[i][1] * w[j][2];
                local grid_v = field_v[base[1] + i, base[2] + j];
                v1 += weight * grid_v;
                C1 += 4 * weight * grid_v * dr' / dh^2; 
            end
        end
        particle.v = v1;
        particle.x += dt * v1;
        particle.J *= 1 + dt * tr(C1);
        particle.C = C1;
    end
end

function solve(self::Model, dt::Float32)
    initGrid(self.grid);
    particle2Grid(self, dt);
    boundaryCondition(self, dt);
    grid2Particle(self, dt);
end

function getPositions(self::Model)
    local r = zeros(length(self.particles), 2);
    local k::Int32 = 1;
    for particle in self.particles
        r[k, :] = particle.x;
        k += 1;
    end
    return r;
end

function getGridPositions(self::Model)
    local grid = self.grid;
    local r = zeros(grid.size_x * grid.size_y, 2);
    local k::Int32 = 1;
    local dh = grid.dh;
    for i = 1 : grid.size_x
        for j = 1 : grid.size_y
            r[k, :] = [(i - 1) * dh, (j - 1) * dh];
            k += 1;
        end
    end
    return r;
end

function getGridVelocities(self::Model)
    local grid = self.grid;
    local v = zeros(grid.size_x * grid.size_y, 2);
    local k::Int32 = 1;
    for i = 1 : grid.size_x
        for j = 1 : grid.size_y
            v[k, :] = grid.field_v[i, j];
            k += 1;
        end
    end
    return v;
end