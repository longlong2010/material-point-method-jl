include("property.jl");
include("geometry.jl");
include("model.jl");

using Plots;
pyplot();
begin
    local n_particle = 8192;
    local grid_size = 128;
    local rho = 1;
    local nu = 0.3;
    local E = 400;
    local dx::Float32 = 1 / grid_size;
    local dt = 2f-4;
    local g = -9.8f0;

    local m = Material(E, nu, rho);
    local grid = Grid2D(grid_size, grid_size, dx);
    local p = Property2D(m);
    local model = Model(grid, g);
    for i  = 1 : n_particle
        local x = Vector{Float32}([rand() * 0.4 + 0.2, rand() * 0.4 + 0.2]);
        local v = Vector{Float32}([0, -1])
        local particle = Particle2D(x, v, (dx * 0.5f0)^2, p);
        addParticle(model, particle);
    end
    for i = 1 : 5000
        solve(model, dt);
        println(i);
        if i % 50 == 0
            local j::Int32 = floor(i / 50);
            local r = getPositions(model);
            local fig = plot(r[:, 1], r[:, 2], seriestype=:scatter, aspect_ratio=:equal, xlims = (0, 1), ylims = (0, 1), size = (800, 800));
            #local r = getGridPositions(model);
            #local v = getGridVelocities(model);
            #fig = quiver(r[:, 1], r[:, 2], quiver=(v[:, 1] * 0.05, v[:, 2] * 0.05), size = (600 * 4, 600 * 4), aspect_ratio=:equal, xlims = (0, 1), ylims = (0, 1));
            savefig(fig, "image/$j.png");
        end
    end
end
