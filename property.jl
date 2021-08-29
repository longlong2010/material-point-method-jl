abstract type Property end;

struct Material <: Property
	E::Float64;
	nu::Float64;
	rho::Float64;
end

struct Property3D <: Property
	M::Material;
end

struct Property2D <: Property
	M::Material;
end
