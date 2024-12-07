def gauss_elimination_with_pivoting(A, b):
   
    n = len(b)

    # Augment the matrix A with vector b
    Ab = []
    for i in range(n):
        Ab.append(A[i] + [b[i]])

    # Forward elimination with partial pivoting
    for i in range(n):
        # Pivoting
        max_row_index = i
        for r in range(i + 1, n):
            if abs(Ab[r][i]) > abs(Ab[max_row_index][i]):
                max_row_index = r
        
        # Swap the rows if needed
        if i != max_row_index:
            Ab[i], Ab[max_row_index] = Ab[max_row_index], Ab[i]

        # Elimination process
        for j in range(i + 1, n):
            factor = Ab[j][i] / Ab[i][i]
            for k in range(i, n + 1):  # Update row j
                Ab[j][k] -= factor * Ab[i][k]

    # Back substitution
    x = [0] * n
    for i in range(n - 1, -1, -1):
        sum_ax = 0
        for j in range(i + 1, n):
            sum_ax += Ab[i][j] * x[j]
        x[i] = (Ab[i][-1] - sum_ax) / Ab[i][i]

    return x

# Example usage
A = [
    [1, 2, 1],
    [2, 6, 1],
    [1, 1, 4]
]
b = [2, 7, 3]

solution = gauss_elimination_with_pivoting(A, b)

# Print the solution using f-string formatting
solution = ', '.join(f'x[{i}] = {value:.4f}'for i, value in enumerate(solution))
print(f"Solution: {solution}")
