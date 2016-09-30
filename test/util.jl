import PipeLayout: nonzero

facts("numerical comparisons") do
    @fact nonzero(0.0) --> false
    @fact nonzero(1e-15) --> false
    @fact nonzero(-1e-15) --> false
    @fact nonzero(1e-3) --> true
    @fact nonzero(-1e-3) --> true
end

facts("random subset sampling") do
    base = collect(1:6)
    sample0 = select_subset(base, 0)
    sample3 = select_subset(base, 3)
    sample6 = select_subset(base, 6)

    @fact length(sample0) --> 0
    @fact length(sample3) --> 3
    @fact unique(sample3) --> sample3
    @fact length(sample6) --> 6
    @fact unique(sample6) --> sample6
    @fact_throws select_subset(base, 9)
end
