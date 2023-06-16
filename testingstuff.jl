using DifferentialEquations

function vector_to_row_matrix(vector::Vector{T}, n::Int) where T
    matrix = reduce(hcat, repeat([vector], n, 1))
    return matrix
end


vector = [1, 2, 3, 4, 5]
n = 5
matrix1 = reduce(hcat, repeat([vector], n, 1))
matrix2 = reducerepeat([vector], n, 1)

display(matrix1)
display(matrix2)
#display(matrix2)

print("Exit code 0")