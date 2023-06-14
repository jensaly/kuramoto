
function meshgrid(x, y)
    X = [x for _ in y, x in x]
    Y = [y for y in y, _ in x]
    X, Y
end

u0 = [0.0, 1.0, 2.0]   # Initial phase values

θi, θj = meshgrid(u0, u0)

vxv = θj - θi
sins = sin.(vxv)
summed = sum(sins, dims=2)
plusu0 = u0 + summed
print(vxv)
print("\n")
print(sins)
print("\n")
print(summed)
print("\n")
print(plusu0)