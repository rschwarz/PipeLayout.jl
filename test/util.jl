import PipeLayout: nonzero

facts("numerical comparisons") do
    @fact nonzero(0.0) --> false
    @fact nonzero(1e-15) --> false
    @fact nonzero(-1e-15) --> false
    @fact nonzero(1e-3) --> true
    @fact nonzero(-1e-3) --> true
end
