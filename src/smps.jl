using JuMP, MathOptInterface
using SparseArrays

"""
A problem template defined by SMPS format.
"""
struct spCorType{T}
    problem_name::String
    template_matrix::AbstractMatrix{T}
    directions::Vector{String}
    row_names::Vector{String}
    col_names::Vector{String}
end

"""
Parse cor file into tokens for intermediate representation.
"""
function _tokenize_cor(cor_path::String)
    tokens = Dict(
        "NAME" => [],
        "ROWS" => [],
        "COLUMNS" => [],
        "RHS" => [],
        "BOUNDS" => []
    )

    # Read cor and remove empty or comment rows
    lines = readlines(cor_path)
    filter!(s -> (!isempty(s) && s[1] != '*'), lines)

    section = ""

    for line in lines
        token = split(line)

        # Look at section line first
        if line[1] != ' '
            # Split each row into tokens first
            section = token[1]

            # special case for NAME section because it is on
            # the same line
            if token[1] == "NAME"
                push!(tokens["NAME"], token[2])
            end
        else
            # Data lines
            push!(tokens[section], token)
        end
    end

    return tokens
end

# Parse row tokens into constraint
# direction and list of row names
function _parse_row_tokens(tokens)
    direction::Vector{Char} = [t[1][1] for t in tokens]
    row_names::Vector{String} = [t[2] for t in tokens]
    return direction, row_names
end

# Extract variable names from column tokens
function _parse_unique_columns(tokens)
    col_names::Vector{String} = [t[1] for t in tokens]
    return unique(col_names)
end

# Extract the coefficient matrix
function _parse_column_to_matrix(tokens, row_names, col_names)
    # Build mapping from the names to the assigned indices
    col_mapping = Dict(col => i for (i, col) in enumerate(col_names))
    row_mapping = Dict(row => i for (i, row) in enumerate(row_names))

    M = spzeros(length(row_names), length(col_names))

    for token in tokens
        col_name = token[1]
        col_index = col_mapping[col_name]
        
        # Iterate the rest of the line in pairs, and populate the entry
        for (row_name, vstring) in Iterators.partition(token[2:end], 2)
            row_index = row_mapping[row_name]
            v = parse(Float64, vstring)
            M[row_index, col_index] = v
        end
    end

    return M
end

# Extract the rhs from RHS tokens, assuming zero entry for unfound
function _parse_rhs(tokens, row_names)
    row_mapping = Dict(row => i for (i, row) in enumerate(row_names))
    rhs = zeros(length(row_names))
    for token in tokens
        for (row_name, vstring) in Iterators.partition(token[2:end], 2)
            row_index = row_mapping[row_name]
            rhs[row_index] = parse(Float64, vstring)
        end
    end
    return rhs
end