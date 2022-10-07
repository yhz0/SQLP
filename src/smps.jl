import JuMP
import MathOptInterface

function read_smps(cor_path::AbstractString, tim_path::AbstractString, sto_path::AbstractString)

end

"Read SMPS cor file cor_path either by MPS format."
function read_cor(cor_path::AbstractString)::JuMP.Model
    return JuMP.read_from_file(cor_path, format = MathOptInterface.FileFormats.FORMAT_MPS)
end

struct spTimeType

end

"Read time file."
function read_tim(tim_path::AbstractString)

end